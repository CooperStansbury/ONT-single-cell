
######## snakemake preamble start (automatically inserted, do not edit) ########
import sys; sys.path.extend(['/home/cstansbu/miniconda3/envs/ont_sc_rna/lib/python3.10/site-packages']); import pickle; snakemake = pickle.loads(b'\x80\x04\x95\x1b\x07\x00\x00\x00\x00\x00\x00\x8c\x10snakemake.script\x94\x8c\tSnakemake\x94\x93\x94)\x81\x94}\x94(\x8c\x05input\x94\x8c\x0csnakemake.io\x94\x8c\nInputFiles\x94\x93\x94)\x81\x94\x8ct/nfs/turbo/umms-indikar/shared/projects/HSC/data/10xBarcoded_SingleCell/P2_run5/blaze_outputs/matched_reads.fastq.gz\x94a}\x94(\x8c\x06_names\x94}\x94\x8c\x05fastq\x94K\x00N\x86\x94s\x8c\x12_allowed_overrides\x94]\x94(\x8c\x05index\x94\x8c\x04sort\x94eh\x12\x8c\tfunctools\x94\x8c\x07partial\x94\x93\x94h\x06\x8c\x19Namedlist._used_attribute\x94\x93\x94\x85\x94R\x94(h\x18)}\x94\x8c\x05_name\x94h\x12sNt\x94bh\x13h\x16h\x18\x85\x94R\x94(h\x18)}\x94h\x1ch\x13sNt\x94bh\x0eh\nub\x8c\x06output\x94h\x06\x8c\x0bOutputFiles\x94\x93\x94)\x81\x94(\x8cA/scratch/indikar_root/indikar1/cstansbu/HSC_p2/fastqc/report.html\x94\x8cG/scratch/indikar_root/indikar1/cstansbu/HSC_p2/fastqc/report_fastqc.zip\x94e}\x94(h\x0c}\x94(\x8c\x04html\x94K\x00N\x86\x94\x8c\x03zip\x94K\x01N\x86\x94uh\x10]\x94(h\x12h\x13eh\x12h\x16h\x18\x85\x94R\x94(h\x18)}\x94h\x1ch\x12sNt\x94bh\x13h\x16h\x18\x85\x94R\x94(h\x18)}\x94h\x1ch\x13sNt\x94bh*h&h,h\'ub\x8c\x06params\x94h\x06\x8c\x06Params\x94\x93\x94)\x81\x94\x8c\x07--quiet\x94a}\x94(h\x0c}\x94h\x10]\x94(h\x12h\x13eh\x12h\x16h\x18\x85\x94R\x94(h\x18)}\x94h\x1ch\x12sNt\x94bh\x13h\x16h\x18\x85\x94R\x94(h\x18)}\x94h\x1ch\x13sNt\x94bub\x8c\twildcards\x94h\x06\x8c\tWildcards\x94\x93\x94)\x81\x94}\x94(h\x0c}\x94h\x10]\x94(h\x12h\x13eh\x12h\x16h\x18\x85\x94R\x94(h\x18)}\x94h\x1ch\x12sNt\x94bh\x13h\x16h\x18\x85\x94R\x94(h\x18)}\x94h\x1ch\x13sNt\x94bub\x8c\x07threads\x94K\x08\x8c\tresources\x94h\x06\x8c\tResources\x94\x93\x94)\x81\x94(K\x08K\x01M\xfahM\xfah\x8c\x04/tmp\x94e}\x94(h\x0c}\x94(\x8c\x06_cores\x94K\x00N\x86\x94\x8c\x06_nodes\x94K\x01N\x86\x94\x8c\x06mem_mb\x94K\x02N\x86\x94\x8c\x07disk_mb\x94K\x03N\x86\x94\x8c\x06tmpdir\x94K\x04N\x86\x94uh\x10]\x94(h\x12h\x13eh\x12h\x16h\x18\x85\x94R\x94(h\x18)}\x94h\x1ch\x12sNt\x94bh\x13h\x16h\x18\x85\x94R\x94(h\x18)}\x94h\x1ch\x13sNt\x94bh^K\x08h`K\x01hbM\xfahhdM\xfahhfh[ub\x8c\x03log\x94h\x06\x8c\x03Log\x94\x93\x94)\x81\x94\x8c\x16logs/fastqc/report.log\x94a}\x94(h\x0c}\x94h\x10]\x94(h\x12h\x13eh\x12h\x16h\x18\x85\x94R\x94(h\x18)}\x94h\x1ch\x12sNt\x94bh\x13h\x16h\x18\x85\x94R\x94(h\x18)}\x94h\x1ch\x13sNt\x94bub\x8c\x06config\x94}\x94(\x8c\x0boutput_path\x94\x8c//scratch/indikar_root/indikar1/cstansbu/HSC_p2/\x94\x8c\x05fastq\x94\x8ct/nfs/turbo/umms-indikar/shared/projects/HSC/data/10xBarcoded_SingleCell/P2_run5/blaze_outputs/matched_reads.fastq.gz\x94\x8c\x08ref_path\x94\x8ci/nfs/turbo/umms-indikar/shared/projects/HSC/data/RefGenome/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz\x94\x8c\x08gtf_path\x94\x8cV/nfs/turbo/umms-indikar/shared/projects/HSC/data/RefGenome/Homo_sapiens.GRCh38.107.gtf\x94\x8c\rminimap2_args\x94\x8c"-ax splice -uf --secondary=no --MD\x94\x8c\x0cumi_distance\x94K\x01u\x8c\x04rule\x94\x8c\x06fastqc\x94\x8c\x0fbench_iteration\x94N\x8c\tscriptdir\x94\x8cFhttps://github.com/snakemake/snakemake-wrappers/raw/v1.29.0/bio/fastqc\x94ub.'); from snakemake.logging import logger; logger.printshellcmds = False; __real_file__ = __file__; __file__ = 'https://github.com/snakemake/snakemake-wrappers/raw/v1.29.0/bio/fastqc/wrapper.py';
######## snakemake preamble end #########
"""Snakemake wrapper for fastqc."""

__author__ = "Julian de Ruiter"
__copyright__ = "Copyright 2017, Julian de Ruiter"
__email__ = "julianderuiter@gmail.com"
__license__ = "MIT"


from os import path
import re
from tempfile import TemporaryDirectory

from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)


def basename_without_ext(file_path):
    """Returns basename of file path, without the file extension."""

    base = path.basename(file_path)
    # Remove file extension(s) (similar to the internal fastqc approach)
    base = re.sub("\\.gz$", "", base)
    base = re.sub("\\.bz2$", "", base)
    base = re.sub("\\.txt$", "", base)
    base = re.sub("\\.fastq$", "", base)
    base = re.sub("\\.fq$", "", base)
    base = re.sub("\\.sam$", "", base)
    base = re.sub("\\.bam$", "", base)

    return base


# Run fastqc, since there can be race conditions if multiple jobs
# use the same fastqc dir, we create a temp dir.
with TemporaryDirectory() as tempdir:
    shell(
        "fastqc {snakemake.params} -t {snakemake.threads} "
        "--outdir {tempdir:q} {snakemake.input[0]:q}"
        " {log}"
    )

    # Move outputs into proper position.
    output_base = basename_without_ext(snakemake.input[0])
    html_path = path.join(tempdir, output_base + "_fastqc.html")
    zip_path = path.join(tempdir, output_base + "_fastqc.zip")

    if snakemake.output.html != html_path:
        shell("mv {html_path:q} {snakemake.output.html:q}")

    if snakemake.output.zip != zip_path:
        shell("mv {zip_path:q} {snakemake.output.zip:q}")
