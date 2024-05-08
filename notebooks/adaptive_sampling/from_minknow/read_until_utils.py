from __future__ import annotations

import logging
import threading
import time
from collections.abc import Generator
from itertools import chain
from typing import Callable, Optional

import numpy as np
from bream4.device_interfaces.devices.base_device import BaseDeviceInterface
from grpc import RpcError
from minknow_api import protocol_service
from minknow_api.acquisition_pb2 import ACQUISITION_RUNNING
from pybasecall_client_lib.helper_functions import package_read
from pybasecall_client_lib.pyclient import PyBasecallClient

BREAM_KEYSTORE = "bream.protocol_communication"
logger = logging.getLogger(__name__)


# All of these indicate that a protocol is still happily going
# So we should keep waiting if it's any of these
# Brought out of function to make sure test can check for any excess states
_PROTOCOL_STATE_TO_WAIT = {
    protocol_service.PROTOCOL_RUNNING,
    protocol_service.PROTOCOL_WAITING_FOR_ACQUISITION,
    protocol_service.PROTOCOL_WAITING_FOR_TEMPERATURE,
}


def wait_for_sequencing_to_start(device: BaseDeviceInterface) -> bool:
    """This function blocks until the phase reported by minknow (via bream) is not unknown/initialising

    Returns true if moved out of unknown/initialising and false if protocol state has finished
    :param device: MinKNOW device wrapper
    """

    skip_phases = {
        protocol_service.PHASE_UNKNOWN,
        protocol_service.PHASE_INITIALISING,
    }

    for msg in device.connection.protocol.watch_current_protocol_run():
        logger.info(f"Phase: {msg.phase} State: {msg.state}")
        logger.debug("Received message: %s" % (msg,))

        if msg.phase not in skip_phases:
            return True
        if msg.state not in _PROTOCOL_STATE_TO_WAIT:
            # Something funny happening here. May be in a test case
            return False

    return False


class WatchAcquisitionStatus(threading.Thread):
    """Class used to stop read until if the acquisition terminates - By bream or other means"""

    def __init__(self, device: BaseDeviceInterface, callback: Callable[[], None]):
        self.device = device
        self.stream_handle = None
        self.logger = logging.getLogger(__name__)
        self.callback = callback

        super().__init__()

    def run(self) -> None:
        self.stream_handle = self.device.connection.acquisition.watch_current_acquisition_run()

        try:
            for msg in self.stream_handle:
                self.logger.info(f"State: {msg.state}")
                self.logger.debug("Received message: %s" % (msg,))

                if msg.state != ACQUISITION_RUNNING:
                    self.stream_handle = None
                    self.callback()
                    return

        except RpcError as exception:
            self.logger.info(f"Received RpcError {exception}. Stopping watcher")

    def stop(self):
        self.logger.info("Stopping watcher stream")
        if self.stream_handle:
            self.stream_handle.cancel()


def basecall(
    guppy_client: PyBasecallClient,
    reads: list,
    dtype: "np.dtype",
    daq_values: dict,
    sample_rate: float,
    basecall_timeout: Optional[int] = None,
) -> Generator[tuple[tuple[int, int], dict], None, None]:
    """Generator that sends and receives data from guppy

    :param caller: pybasecall_client_lib.pyclient.PyBasecallClient
    :param reads: List of reads from read_until
    :param dtype: Numpy dtype to cast buffer to
    :param daq_values: Dict of channel -> NamedTuple(offset/scaling)
    :param basecall_timeout: Wait up to x seconds for the basecaller to basecall reads.
                             Skip the rest of the reads. If 0 or none, no timeout.
    :param sample_rate: int

    :returns:
        - read_info (:py:class:`tuple`) - channel (int), read number (int)
        - read_data (:py:class:`dict`) - Data returned from Guppy
    :rtype: Iterator[ ((channel, read_number), BasecallInfoDict) ]
    """

    hold = {}
    packaged_reads = []
    reads_waiting = len(reads)
    attempts = 3

    # Package reads into a list to send to guppy
    for (channel, read) in reads:
        packaged = package_read(
            read_id=read.id,
            raw_data=np.frombuffer(read.raw_data, dtype),
            daq_offset=daq_values[channel].offset,
            daq_scaling=daq_values[channel].scaling,
            start_time=int(read.start_sample),
            sampling_rate=sample_rate,
        )
        packaged_reads.append(packaged)
        hold[read.id] = (channel, read.number)

    # Try 3 times to send to guppy
    for _ in range(attempts):
        if guppy_client.pass_reads(packaged_reads):
            break
    else:
        logger.info(f"Guppy rejected {len(packaged_reads)} reads {attempts} times. Giving up.")
        return

    start_time = time.monotonic()

    # Retrieve result for all that guppy acknowledged. If timeout reached, return
    while reads_waiting > 0:
        # Add a timeout in case of guppy problems
        if basecall_timeout and time.monotonic() - start_time > basecall_timeout:
            logger.info(
                f"Still waiting for {reads_waiting} reads from the basecaller. "
                + f"Timeout reached at {basecall_timeout} seconds, so skipping those reads."
            )
            return

        split_results = guppy_client.get_completed_reads()
        results = list(chain.from_iterable(split_results))
        basecall_time = time.monotonic()

        if not results:
            # Give guppy a moment to process more reads
            time.sleep(guppy_client.throttle)
            continue

        for read in results:
            read_id = read["metadata"]["read_id"]
            try:
                (channel, read_number) = hold[read_id]
                read["metadata"]["send_read_time"] = start_time
                read["metadata"]["receive_read_time"] = basecall_time
                reads_waiting -= 1
                yield ((channel, read_number), read)
            except KeyError:
                logger.info(f"Read not requested: {read_id}")
