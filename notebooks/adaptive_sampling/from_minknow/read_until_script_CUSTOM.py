from __future__ import annotations

import argparse
import ast
import json
import logging
import os
import pprint
import random
import signal
import threading
import time
from collections import defaultdict
from typing import Any, Callable, Optional, TypeAlias

import minknow
import minknow_api.manager
import numpy as np
import pandas as pd
from bream4.device_interfaces import device_wrapper
from bream4.device_interfaces.devices.base_device import BaseDeviceInterface, get_env_grpc_port
from bream4.toolkit.procedure_components.output_locations import get_run_output_path
from google.protobuf import json_format
from grpc import RpcError
from minknow_api import protocol_service
from minknow_api.statistics_pb2 import ReadEndReason
from pybasecall_client_lib.pyclient import PyBasecallClient
from read_until import AccumulatingCache, ReadUntilClient
from sequencing.read_until.read_until_utils import WatchAcquisitionStatus, basecall, wait_for_sequencing_to_start
from utility.config_argparse import boolify
from utility.config_file_utils import expand_grouped_list

"""
Basic Read-Until script.

Note that because this uses a bream device_wrapper you can do log_to_gui etc.
The logs will also appear in the relevant device position logs.

If any log output is generated, due to how log rotation works on windows vs unix,
 * a separate log file (bream-0) will be created on unix
 * the same bream log as the main process will be used on Windows.

Due to the fact there is a bream process already interacting with the device,
please alter this at your own peril.

The current implementation of bream does no flicking of strand so we should be
safe to do this.
"""

GUPPY_PARAMS_RENAME = {
    # minknow -> guppy
    "barcoding_kits": "barcode_kits",
    "min_score": "min_score_barcode_front",
    "min_score_rear": "min_score_barcode_rear",
    "min_score_mid": "min_score_barcode_mid",
    "min_score_mask": "min_score_barcode_mask",
}

ReadType: TypeAlias = "dict[str, Any]"  # The read that comes back from the basecall client


class ReadUntil(threading.Thread):
    """This class is the meat of Read Until. This is the class that hooks in to the read until api
    and basecalls reads and determines whether to flick or accept the reads based on configuration passed in.

    * run is the main method that is on a separate thread and does all the lifting.
    * _enrich/_deplete get called depending on config to give a response of what to do
        - unblock/unblock_hit_outside_bed/stop_receiving/stop_receiving_hit_outside_bed/no_decision can be returned

    """

    def __init__(
        self,
        config: dict[str, Any],
        device: BaseDeviceInterface,
        read_until_client: ReadUntilClient,
        guppy_client: PyBasecallClient,
        barcode_counts: Optional[dict[str, int]] = None,
        exception_callback: Optional[Callable] = None,
    ):
        self.logger = logging.getLogger(__name__)
        self.exception_callback = exception_callback
        self.exception = None
        self.keep_going = True
        self.guppy_client = guppy_client
        self.read_until_client = read_until_client

        self.flick_time = (
            time.monotonic() + config["barcode_balance_after_x_s"] if not barcode_counts else time.monotonic()
        )

        self.config = config
        self.barcode_balance_min_percent = config["barcode_balance_min_percent"] / 100
        if config["filter_type"] not in ["enrich", "deplete", "barcode_balance"]:
            raise ValueError(f"{config['filter_type']} not enrich, deplete, or barcode_balance")
        self.mode = config["filter_type"]

        if barcode_counts is None:
            self.barcode_counts: dict[str, int] = defaultdict(int)
            self.all_barcode_count = 0
            self.median_barcode_count = None
            self.barcode_inv_freq = {}
        else:
            self.barcode_counts = barcode_counts
            self._recalculate_barcode_balance_stats()

        self.report_read_path = config["report_read"]
        if self.report_read_path:
            self.report_read_path = get_run_output_path(device, "adaptive_sampling.csv")
            self.exists = os.path.exists(self.report_read_path)

        self.report_sam_path = config["report_sam"]
        if self.report_sam_path:
            self.report_sam_path = get_run_output_path(device, "adaptive_sampling_sam.csv")

        self.sample_rate = device.get_sample_rate()
        super().__init__()

    def _recalculate_barcode_balance_stats(self) -> None:
        """It can be expensive to do certain calculations for every read, especially if there are lots of
        barcodes present. This function recalculates stats used in _barcode_balance
        """

        self.all_barcode_count = sum(self.barcode_counts.values())
        self.median_barcode_count = np.median(list(self.barcode_counts.values()))

        self.barcode_inv_freq = {
            k: self.all_barcode_count / v
            for k, v in self.barcode_counts.items()
            if (v / self.median_barcode_count) >= self.barcode_balance_min_percent
        }

    def _barcode_balance(self, read: ReadType) -> str:
        """Attempts to balance the barcode in the read by maybe flicking if that barcode has
        have already been sequenced a lot

        :param barcode: String of the barcode to deal with
        :returns: unblock*/stop_receiving*/no_decision to indicate what to do

        :rtype: str

        """
        barcode = read["metadata"]["barcode_arrangement"]

        # Some old kits report barcode 12 as 12a. Make sure to unify them
        if barcode == "barcode12a":
            barcode = "barcode12"
        # Can remove above if we then start changing padding
        if barcode == "barcode012a":
            barcode = "barcode012"

        # If we don't know it get rid of it
        if barcode == "unclassified" or barcode is None:
            return "unblock_barcode_unclassified"

        if read["metadata"]["barcode_score"] < self.config["min_score_filter"]:
            return "stop_receiving_low_score"

        if self.config["barcode_balance_barcodes"]:
            if barcode not in self.config["barcode_balance_barcodes"]:
                return "unblock_bad_barcode"

        self.barcode_counts[barcode] += 1

        # If we have no data in our cache, blanket accept
        # Or not yet allowed to flick
        if self.all_barcode_count == 0 or time.monotonic() < self.flick_time:
            return "stop_receiving"

        assert self.median_barcode_count is not None  # nosec Type hint
        proportion = self.barcode_counts[barcode] / self.median_barcode_count

        # If we haven't accepted many of this barcode, accept it
        if proportion < self.barcode_balance_min_percent or barcode not in self.barcode_inv_freq:
            return "stop_receiving"

        # If we've seen a few then accept it inversely proportional to how many we've seen
        accept_prob = self.barcode_inv_freq[barcode] / max(self.barcode_inv_freq.values())
        if accept_prob > random.uniform(0, 1):
            return "stop_receiving"

        return "unblock"

    def _enrich(self, read: ReadType) -> str:
        """Figures out what to do with a read if trying to enrich from the ref/bed file.

        It accepts a read from guppy - This assumes it is a dict with metadata as a key
        where the relevant information can be accessed.

        Returns a string of the decision (unblock*/stop_receiving*/no_decision)
        """

        if self.config["bed_file"]:
            hits = read["metadata"]["alignment_bed_hits"]
            if hits:
                return "stop_receiving"
            elif read["metadata"]["sequence_length"] >= self.config["enrich_unblock_min_sequence_length"]:
                return "unblock"
            elif read["metadata"]["alignment_genome"] != "*":
                return "unblock_hit_outside_bed"
        else:
            hits = read["metadata"]["alignment_genome"]
            if hits and hits != "*":
                return "stop_receiving"
            elif read["metadata"]["sequence_length"] >= self.config["enrich_unblock_min_sequence_length"]:
                return "unblock"

        return "no_decision"

    def _deplete(self, read: ReadType) -> str:
        """Figures out what to do with a read if trying to deplete from the ref/bed file.

        It accepts a read from guppy - This assumes it is a dict with metadata as a key
        where the relevant information can be accessed.

        Returns a string of the decision (unblock*/stop_receiving*/no_decision)
        """
        if self.config["bed_file"]:
            hits = read["metadata"]["alignment_bed_hits"]
            if hits:
                return "unblock"
            elif read["metadata"]["sequence_length"] >= self.config["deplete_stop_receiving_min_sequence_length"]:
                return "stop_receiving"
            elif read["metadata"]["alignment_genome"] != "*":
                return "stop_receiving_hit_outside_bed"

        else:
            hits = read["metadata"]["alignment_genome"]
            if hits and hits != "*":
                return "unblock"
            elif read["metadata"]["sequence_length"] >= self.config["deplete_stop_receiving_min_sequence_length"]:
                return "stop_receiving"

        return "no_decision"

    def run(self) -> None:
        """Wrap the run method to record if we have an exception and to signal to exit"""
        try:
            self._run()
            return
        except RpcError as exc:
            # If a stream was purposefully cancelled, then not our error
            code = exc.code
            if code == code.CANCELLED or code == code.ABORTED:
                return
            self.logger.info("rpc_error", exc_info=True)
            self.exception = exc
        except Exception as exc:
            self.logger.info("error", exc_info=True)
            self.exception = exc

        if self.exception_callback:
            self.exception_callback()
        raise self.exception

    def _run(self) -> None:
        """Constantly:

        1. Pulls the latest reads from the read until api client.
        2. Basecalls them
        3. For each read, figures out whether to flick/accept/leave
        4. Potentially write reports

        To trigger a stop, call self.stop()

        """
        info_time = time.monotonic() + 60
        read_basecalled_count = 0
        read_chunk_count = 0
        sample_rate = float(self.sample_rate)

        while self.keep_going:

            batch_time = time.time()  # Used to track acquisition timings in CSV - Needs absolute time
            start_time = time.monotonic()  # Everything else tracked as monotonic to ignore clock changes

            # get the most recent read chunks from the client
            read_chunks = self.read_until_client.get_read_chunks(batch_size=self.config["batch_size"], last=True)

            # Move any "bad" chunks into a separate list so not basecalled
            bad_chunks = []
            if self.config["flick_strand_once"]:
                # Check the mux change end reason of the previous read. 2 scenarios (Processing * read chunk)
                # <strand> <flick> *<strand:previous_end_reason=mux_change>
                # This is most likely the same strand so flick failed, so just accept this read
                # <strand> <flick> <pore:previous_end_reason=mux_change> *<strand: previous_end_reason!=mux_change>
                # This is a new strand so we want to continue with AS logic
                reason = ReadEndReason.MuxChange
                bad_chunks.extend([(ch, read) for (ch, read) in read_chunks if read.previous_read_end_reason == reason])
                read_chunks = [(ch, read) for (ch, read) in read_chunks if read.previous_read_end_reason != reason]

            read_chunk_count += len(read_chunks)
            chunk_time = time.monotonic()

            called_batch = basecall(
                guppy_client=self.guppy_client,
                reads=read_chunks,
                dtype=self.read_until_client.signal_dtype,
                daq_values=self.read_until_client.calibration_values,
                basecall_timeout=self.config["basecall_timeout"],
                sample_rate=sample_rate,
            )

            # Tell progress
            if info_time <= start_time:
                self.logger.info(
                    f"Basecalled {read_basecalled_count} reads and received {read_chunk_count} reads in the last minute"
                )
                read_basecalled_count = 0
                read_chunk_count = 0
                info_time = start_time + 60

            read_until_data = []
            sam_data = []

            unblock_reads = []
            stop_receiving_reads = []

            if self.mode == "barcode_balance":
                self._recalculate_barcode_balance_stats()

            for (channel, read_number), read in called_batch:
                read_basecalled_count += 1

                if read["metadata"]["sequence_length"] <= 0:
                    continue  # Bad basecall, some metadata may not be present

                if self.mode == "enrich":
                    decision = self._enrich(read)
                elif self.mode == "deplete":
                    decision = self._deplete(read)
                elif self.mode == "barcode_balance":
                    decision = self._barcode_balance(read)
                else:
                    assert False  # nosec Type checking helper

                if decision.startswith("stop_receiving"):
                    stop_receiving_reads.append((channel, read_number))
                elif decision.startswith("unblock"):
                    unblock_reads.append((channel, read_number))

                if self.report_read_path:
                    read_until_data.append(
                        [
                            batch_time,
                            read_number,
                            channel,
                            read["metadata"]["duration"],
                            read["metadata"]["read_id"],
                            read["metadata"]["sequence_length"],
                            decision,
                            read["metadata"].get("barcode_arrangement", "N/A"),
                            read["metadata"]["mean_qscore"],
                            round(chunk_time - start_time, 3),
                            round(read["metadata"]["send_read_time"] - chunk_time, 3),
                            round(read["metadata"]["receive_read_time"] - read["metadata"]["send_read_time"], 3),
                            round(time.monotonic() - read["metadata"]["receive_read_time"], 3),
                        ]
                    )

                if self.report_sam_path:
                    sam_data.append(str(read["metadata"]["alignment_sam_record"]) + "\n")

            # Ignore the rest of the bad_chunk reads
            for (channel, read) in bad_chunks:
                stop_receiving_reads.append((channel, read.number))
                # Put the bad reads into data so that there is a record of them
                if self.report_read_path:
                    read_until_data.append(
                        [
                            batch_time,
                            read.number,
                            channel,
                            len(read.raw_data),
                            read.id,
                            0,  # not basecalled, no sequence length
                            "stop_receiving_bad_previous_flick",
                            "N/A",  # not basecalled, no barcode
                            0,  # not basecalled, no qscore
                            0,  # Not basecalled, no timings for basecalling
                            0,  # ..
                            0,  # ..
                            0,  # ..
                        ]
                    )

            time_to_make_decisions = time.monotonic()

            if unblock_reads:
                self.read_until_client.unblock_read_batch(unblock_reads, duration=self.config["unblock_duration"])
            if stop_receiving_reads:
                self.read_until_client.stop_receiving_batch(stop_receiving_reads)

            time_to_send_decisions = time.monotonic()

            if read_until_data:
                columns = [
                    "batch_time",
                    "read_number",
                    "channel",
                    "num_samples",
                    "read_id",
                    "sequence_length",
                    "decision",
                    "barcode_arrangement",
                    "mean_qscore",
                    "time_to_get_chunks",
                    "time_to_package_and_send_reads",
                    "time_to_basecall_reads",
                    "time_to_make_decisions",
                ]

                df = pd.DataFrame(read_until_data, columns=columns)
                df["time_to_send_decisions"] = round(time_to_send_decisions - time_to_make_decisions, 3)

                if self.mode != "barcode_balance":
                    # mean_qscore is present in non barcoding, but not needed for debug
                    df.drop(["mean_qscore", "barcode_arrangement"], axis=1, inplace=True)

                if not self.config["debug_metrics"]:
                    # Remove any columns starting with time_to as they are debug metrics
                    df.drop(list(df.columns[df.columns.str.startswith("time_to_")]), axis=1, inplace=True)

                df.to_csv(self.report_read_path, mode="a", header=not self.exists, index=False)
                self.exists = True

            if sam_data:
                with open(self.report_sam_path, "a") as f:
                    f.writelines(sam_data)

    def stop(self) -> None:
        """Signals to the thread (run method) to stop.

        This can take a while to actually stop because basecalling for that round
        will need to be finished before the stop signal will take effect
        """
        self.keep_going = False


class ReadUntilManager:
    """Class used to start/stop ReadUntil process as necessary"""

    def __init__(self, device: BaseDeviceInterface, custom_settings: dict[str, Any]):
        self.device = device
        self.custom_settings = custom_settings
        self.read_until_process = None
        self.logger = logging.getLogger(__name__)
        self.keep_going = True
        self.watch_acquisition_status = None
        self.watch_protocol_stream = None
        self.barcode_counts = defaultdict(int)  # Used to cache between read until starts
        self.exception = None
        self._lock = threading.RLock()

        self.read_until_client = None

        self.custom_settings = custom_settings
        if custom_settings["filter_type"] == "barcode_balance":
            custom_settings["reference_files"] = [""]

        if len(custom_settings["reference_files"]) != 1:
            raise RuntimeError("Only one alignment reference file is supported")

        # Get client handle
        self.guppy_client = PyBasecallClient(
            address=self._get_guppy_conn_str(),
            config=custom_settings["guppy_config"],
            align_ref=custom_settings["reference_files"][0],
            minimap_opt_string="-ax sr -k 15 -w 2 -n 1 -m 10",
            bed_file=custom_settings["bed_file"],
            server_file_load_timeout=180,
            priority=PyBasecallClient.high_priority,
            retries=3,
            throttle=0.01,
        )

    def _get_guppy_conn_str(self) -> str:
        manager = minknow_api.manager.Manager(port=self.device.provided_manager_port)
        guppy_info = manager.get_guppy_connection_info()
        if guppy_info.port is not None:
            return f"127.0.0.1:{guppy_info.port}"
        else:
            return f"ipc://{guppy_info.ipc_path}"

    def trigger_guppy_connect(self) -> None:
        if self.custom_settings["filter_type"] == "barcode_balance":
            # Set up barcoding options to guppy
            opts = self._get_barcoding_connection_options()
            self.guppy_client.set_params(opts)

            if "barcode_kits" not in opts:
                raise RuntimeError("Barcode balancing specified, but no barcode kits enabled")

        # This will block until the guppy server has loaded the model, index,
        #   and bed file to get around this, we could init the basecall client in another
        #   thread
        self.guppy_client.connect()
        self.logger.info("Guppy client connected")

    def _get_barcoding_connection_options(self) -> dict[str, Any]:
        # Gets the barcoding parameters to pass to guppy
        # Retrieves the options from MinKNOW
        basecall_info = self.device.connection.analysis_configuration.get_basecaller_configuration()
        basecall_config = json.loads(
            json_format.MessageToJson(
                basecall_info, preserving_proto_field_name=True, including_default_value_fields=True
            )
        )

        opts = basecall_config.get("barcoding_configuration", {})
        if opts:
            for (key, new_key) in GUPPY_PARAMS_RENAME.items():
                if key in opts:
                    opts[new_key] = opts.pop(key)

            if "trim_barcodes" in opts:
                # Guppy negated option
                opts["enable_trim_barcodes"] = opts.pop("trim_barcodes")

            if "min_score_barcode_front" in opts and "min_score_barcode_rear" not in opts:
                opts["min_score_barcode_rear"] = opts["min_score_barcode_front"]

            # Don't care about this option
            opts.pop("ignore_unspecified_barcodes", None)

        if "barcoding_kits" in self.custom_settings:
            opts["barcode_kits"] = self.custom_settings["barcoding_kits"]

        return opts

    def _exit(self, *args) -> None:
        # Can get called from multiple places as it's a registered handler
        with self._lock:
            self.logger.info("Received exit signal")
            self.keep_going = False

            if self.watch_protocol_stream:
                self.watch_protocol_stream.cancel()
            if self.watch_acquisition_status:
                self.watch_acquisition_status.stop()

            self._stop_read_until()

            self.guppy_client.disconnect()

    def _stop_read_until(self) -> None:
        # To be safe also wrap the one this calls in a lock (Rlock so can reenter for free)
        with self._lock:
            if self.read_until_process:
                self.logger.info("Stopping read until process")
                self.read_until_process.stop()

                if self.read_until_process:
                    self.logger.info("Joining read until process")
                    try:
                        self.read_until_process.join()
                    except RuntimeError:
                        pass  # We could get in the situation where we join with ourselves
                    if not self.exception:
                        self.exception = self.read_until_process.exception
                if self.read_until_client:
                    self.logger.info("Resetting read until client")
                    self.read_until_client.reset()
                    self.read_until_client = None

                self.logger.info("Read until process ended")
                self.read_until_process = None

    def _start_read_until(self) -> None:
        with self._lock:
            if self.read_until_process is None:

                # Make sure any reads waiting in basecaller are drained as we don't need them
                self.logger.info("Draining any straggling basecaller reads")
                while self.guppy_client.get_completed_reads():
                    pass

                self.read_until_client = ReadUntilClient(
                    mk_host="localhost",
                    mk_port=get_env_grpc_port(),
                    cache_type=AccumulatingCache,
                    one_chunk=False,
                    filter_strands=True,
                    calibrated_signal=False,
                    prefilter_classes=set(self.custom_settings.get("accepted_first_chunk_classifications", [])),
                )
                self.read_until_client.run(
                    first_channel=self.custom_settings.get("first_channel", 1),
                    last_channel=self.custom_settings.get("last_channel", self.device.channel_count),
                    max_unblock_read_length_samples=self.custom_settings.get("max_unblock_read_length_samples", 0),
                    accepted_first_chunk_classifications=self.custom_settings.get(
                        "accepted_first_chunk_classifications", []
                    ),
                )

                self.logger.info("Starting read until process")
                self.read_until_process = ReadUntil(
                    self.custom_settings,
                    self.device,
                    self.read_until_client,
                    self.guppy_client,
                    barcode_counts=self.barcode_counts,
                    exception_callback=self._exit,
                )
                self.read_until_process.start()

    def run(self) -> None:
        """This coordinates the start/stop of the ReadUntil class.

        It will ensure sequencing has started (Wait for a sequencing process to start and be acquiring data)
        And then follow that run.

        There are several things that can cause this class/method to exit:
        1. Hook into sigterm to allow read until processes to exit
        2. Terminates if acquisition status of the run changes from RUNNING
        3. Receive a CANCELLED/ABORTED

        There are 2 modes to running:
        * config["run_in_mux_scan"] = True -> Just run read until until one of the exit methods above get triggered
        * = False, Start ReadUntil only outside of mux scans (Watches the keystore to figure this out)
        """
        # Attach triggers to stop
        signal.signal(signal.SIGTERM, self._exit)

        if not wait_for_sequencing_to_start(self.device):
            return

        # Move this to after waiting for sequencing to start to ensure all config options have been applied
        self.trigger_guppy_connect()

        self.watch_acquisition_status = WatchAcquisitionStatus(self.device, self._exit)
        self.watch_acquisition_status.start()

        while self.keep_going:

            if not self.custom_settings["run_in_mux_scan"]:

                # Now watch the phase to start/stop the read until feature
                try:
                    self.watch_protocol_stream = self.device.connection.protocol.watch_current_protocol_run()
                    for msg in self.watch_protocol_stream:
                        self.logger.info("Received state %s" % (protocol_service.ProtocolPhase.Name(msg.phase),))

                        if msg.phase in {
                            protocol_service.PHASE_MUX_SCAN,
                            protocol_service.PHASE_PREPARING_FOR_MUX_SCAN,
                        }:
                            self.logger.info("Stopping Read Until as entering a mux scan")
                            self._stop_read_until()

                        elif msg.phase == protocol_service.PHASE_SEQUENCING:
                            self.logger.info("Starting Read Until as starting sequencing")
                            self._start_read_until()

                except RpcError as exception:
                    self.logger.info(f"Received exception {exception} on stream")

                    code = exception.code()

                    # Check if it wasn't a blip
                    if code == code.CANCELLED or code == code.ABORTED:
                        self._stop_read_until()
                        break

            else:
                self.logger.info("Starting Read Until")
                self._start_read_until()
                self.read_until_process.join()  # type: ignore
                self._stop_read_until()
                break

        # If inner read until process has a problem, re-raise
        if self.exception:
            raise self.exception


def main(config: dict[str, Any], device: Optional[BaseDeviceInterface] = None) -> None:
    if device is None:
        device = device_wrapper.create_grpc_client()

    if config["last_channel"] is None:
        config["last_channel"] = device.channel_count

    logger = logging.getLogger(__name__)
    logger.info(pprint.pformat(config))

    read_until_manager = ReadUntilManager(device, config)
    read_until_manager.run()


def eval_wrapper(string: str) -> Any:
    # Literal eval will gobble up an escaped sequence, so make sure to double escape
    return ast.literal_eval(string.replace("\\", "\\\\"))


def parse_barcode_list(string: str) -> list[str]:
    # Given: "[[1, 3], [5, 5], 7]" return [barcode01, barcode02, barcode03, barcode05, barcode07]"
    group_list = eval_wrapper(string)
    expanded_group_list = expand_grouped_list(group_list)

    return [f"barcode{x:02}" for x in expanded_group_list]


def parse_args(args: list) -> argparse.Namespace:
    parser = argparse.ArgumentParser()

    parser.add_argument("--filter_type", required=True, choices=["enrich", "deplete", "barcode_balance"])
    parser.add_argument(
        "--reference_files",
        help="List of alignments (Only one supported). Required if filter_type != barcode_balance",
        type=eval_wrapper,
    )
    parser.add_argument("--guppy_config", help="Guppy config to use to basecall reads", type=str, required=True)
    parser.add_argument("--batch_size", default=512, help="How many new reads to grab from read until client", type=int)
    parser.add_argument("--bed_file", help="Path to bed file", default="''", type=eval_wrapper)
    parser.add_argument("--report_sam", default=False, type=boolify, help="Whether to output sam data")
    parser.add_argument("--report_read", default=True, type=boolify, help="Whether to output read decisions")
    parser.add_argument("--run_in_mux_scan", default=False, type=boolify, help="Whether to run read until in mux scans")
    parser.add_argument(
        "--barcoding_kits", help="Override what barcode kits to pass to guppy for basecalling", type=eval_wrapper
    )
    parser.add_argument(
        "--enrich_unblock_min_sequence_length",
        type=int,
        default=200,
        help="Enrich will wait up to x bases to make the decision to unblock the read",
    )
    parser.add_argument(
        "--deplete_stop_receiving_min_sequence_length",
        type=int,
        default=4000,
        help="Deplete will wait up to x bases to make the decision to accept the read",
    )
    parser.add_argument(
        "--max_unblock_read_length_samples",
        type=int,
        default=0,
        help="Maximum read length MinKNOW will attempt to unblock (in samples)",
    )
    parser.add_argument(
        "--accepted_first_chunk_classifications",
        nargs="*",
        default=["strand", "adapter"],
        help="RU will only stream reads that start with one of these classifications. "
        + "All others will be _accepted_ + not streamed",
    )
    parser.add_argument(
        "--barcode_balance_after_x_s",
        default=1800,
        type=int,
        help="Don't consider balancing until x seconds have passed.",
    )
    parser.add_argument(
        "--barcode_balance_min_percent",
        default=1,
        type=int,
        help="Don't consider flicking until barcode has > x% in the median barcode counts." + "Default x=1 for 1%",
    )
    parser.add_argument(
        "--min_score_filter", type=float, default=60.0, help="Minimum score for barcodes to go through decision process"
    )
    parser.add_argument(
        "--barcode_balance_barcodes",
        type=parse_barcode_list,
        help="If present, only balance these barcodes (in the form: [[1, 3], 12, [13, 13]] means 1, 2, 3, 12, 13)",
    )
    parser.add_argument(
        "--first_channel",
        type=int,
        default=1,
        help="Start of the range of channels that read until will work with. Defaults to all channels",
    )
    parser.add_argument(
        "--last_channel",
        type=int,
        help="End of the range of channels that read until will work with. Defaults to all channels",
    )
    parser.add_argument("--unblock_duration", type=float, default=0.1, help="How long to unblock reads for")
    parser.add_argument(
        "--basecall_timeout", type=int, default=60, help="How long to wait for guppy to basecall a batch of reads"
    )
    parser.add_argument(
        "--flick_strand_once",
        default=False,
        type=boolify,
        help="If strand fails to flick, accept the strand instead of reflicking",
    )

    parser.add_argument(
        "--debug_metrics",
        type=boolify,
        default=False,
        help="Whether to include timing metrics if the report is enabled",
    )
    return parser.parse_args(args)


@minknow.main
def read_until_script(local: minknow.ScriptLocalState):
    args = parse_args(local.argv[1:])
    main(config=vars(args))