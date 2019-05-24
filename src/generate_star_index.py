"""
This script is a part of the ENCODE micro-rna-pipeline.
This script generates index for STAR aligner.
"""

__author__ = "Otto Jolanki"
__version__ = "0.1.0"
__license__ = "MIT"

import argparse
import logging
import os
import shlex
import shutil
import subprocess
import tarfile

# logging config is done in a primitive way to avoid passing config files around to wdl tasks

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
filehandler = logging.FileHandler("generate_star_index.log")
filehandler.setLevel(logging.DEBUG)
consolehandler = logging.StreamHandler()
consolehandler.setLevel(logging.INFO)
formatter = logging.Formatter("%(asctime)s | %(levelname)s | %(name)s: %(message)s")
filehandler.setFormatter(formatter)
consolehandler.setFormatter(formatter)
logger.addHandler(consolehandler)
logger.addHandler(filehandler)


def reset_tar_info(tarinfo):
    """Resets Tarinfo object.

    To be used as a filter argument when adding into a Tarfile archive.
    Resets attributes described below for the purpose of being able to get
    identical md5 sums for tar archives that have identical contents, but
    have been created at different times by possibly different users.

    Args:
        tarinfo: A Tarinfo object.

    Returns:
        A Tarinfo object with gid, uid and mtime set to 0, and gname and
        uname set to 'root'.
    """
    tarinfo.gid = tarinfo.uid = 0
    tarinfo.gname = tarinfo.uname = "root"
    tarinfo.mtime = 0
    return tarinfo


def make_tar_archive_from_dir(input_dir, output_filename, filter=None):
    """Make tar archive from a directory with option to apply a filter.

    Args:
        input_dir: Path to directory to be archived e.g. ~/inputs
        output_filename: filename where output archive is written to
        filter: function that takes TarInfo object as an input and returns
                a TarInfo object. If filter returns None instead the TarInfo
                object will be excluded from the archive.
    """
    with tarfile.open(output_filename, "w") as out_tar:
        out_tar.add(input_dir, filter=filter)
    return None


def compress_with_pigz(input_filename, output_filename, ncpus):
    """Run subprocess call to pigz.

    Args:
        input_filename
        output_filename
        ncpus
    """
    command = "pigz -c -n -p {ncpus} {input}".format(ncpus=ncpus, input=input_filename)
    with open(output_filename, "w") as f:
        subprocess.call(shlex.split(command), stdout=f)
    return None


def make_star_index(ncpus, annotation, refgenome, outfile):
    """Build index with STAR.

    Args:
        ncpus: Int number of threads
        annotation: Path to annotation .gtf file
        refgenome: Path to reference genome .fasta
        outfile: Output filename for index tar.gz archive
    """
    STAR_COMMAND_TEMPLATE = """STAR \
    --runThreadN {ncpus} \
    --runMode genomeGenerate \
    --genomeDir out \
    --sjdbGTFfile {annotation} \
    --sjdbOverhang 1 \
    --genomeFastaFiles {refgenome}
    """
    command = STAR_COMMAND_TEMPLATE.format(
        ncpus=ncpus, annotation=annotation, refgenome=refgenome
    )
    logger.info("Creating temporary directory out/")
    os.mkdir("out")
    logger.info("Running STAR command %s", command)
    subprocess.call(shlex.split(command))
    logger.info("Index building success.")
    output_tar = os.path.splitext(outfile)[0]
    logger.info("Building tar archive %s", output_tar)
    make_tar_archive_from_dir("out", output_tar, reset_tar_info)
    logger.info("Compressing tar %s into %s", output_tar, outfile)
    compress_with_pigz(output_tar, outfile, ncpus)
    logger.info("Removing the uncompressed tar archive and temporary directory.")
    os.remove(output_tar)
    shutil.rmtree("out")
    return None


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--ncpus", type=int, help="Number of threads")
    parser.add_argument(
        "--annotation_file", type=str, help="Filename for annotation .gtf"
    )
    parser.add_argument(
        "--genome_file", type=str, help="Filename for reference genome .fasta"
    )
    parser.add_argument(
        "--output_file",
        type=str,
        help="Filename for the output. Format is gzipped tar (tar.gz or .tgz are the recommended suffixes).",
    )
    args = parser.parse_args()

    make_star_index(
        args.ncpus, args.annotation_file, args.genome_file, args.output_file
    )
