"""
This script is a part of the ENCODE micro-rna-pipeline.
This script generates index for STAR aligner.
"""

__author__ = 'Otto Jolanki'
__version__ = '0.1.0'
__license__ = 'MIT'

import argparse
import gzip
import logging
import os
import shlex
import subprocess
import tarfile

# logging config is done in a primitive way to avoid passing config files around to wdl tasks

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
filehandler = logging.FileHandler('generate_star_index.log')
filehandler.setLevel(logging.DEBUG)
consolehandler = logging.StreamHandler()
consolehandler.setLevel(logging.INFO)
formatter = logging.Formatter(
    '%(asctime)s | %(levelname)s | %(name)s: %(message)s')
filehandler.setFormatter(formatter)
consolehandler.setFormatter(formatter)
logger.addHandler(consolehandler)
logger.addHandler(filehandler)

STAR_COMMAND_TEMPLATE = '''STAR \
    --runThreadN {ncpus} \
    --runMode genomeGenerate \
    --genomeDir out/ \
    --sjdbGTFfile {annotation} \
    --sjdbOverhang 1 \
    --genomeFastaFiles {refgenome}
'''


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
    tarinfo.gname = tarinfo.uname = 'root'
    tarinfo.mtime = 0
    return tarinfo


def make_tar_archive_from_dir(input_dir, output_filename, filter=None):
    with tarfile.open(output_filename, 'w') as out_tar:
        out_tar.add(input_dir, filter=filter)
    return None


def compress_without_header(input_filename, output_filename):
    with open(output_filename, 'wb') as out_gz:
        with open(input_filename, 'rb') as src_file:
            with gzip.GzipFile('', 'wb', fileobj=out_gz, mtime=0) as gz_dest:
                gz_dest.writelines((line for line in src_file))
    return None


def call_star(ncpus, annotation, refgenome, outfile):
    output_basename = os.path.splittext(outfile)[0]
    command = STAR_COMMAND_TEMPLATE.format(ncpus=ncpus, annotation=annotation,
                                           refgenome=refgenome)
    logger.info('Running STAR command %s', command)
    subprocess.call(shlex.split(command))
    output_tar = output_basename + '.tar'
    make_tar_archive_from_dir('out', output_tar, reset_tar_info)
    compress_without_header(output_tar, outfile)
    return None


def main(args):
    call_star(args.ncpus, args.annotation_file, args.genome_file,
              args.output_file)
    with open('staridx_test2.tar.gz', 'wb') as out_gz:
        with open('staridx_test2.tar', 'rb') as src_tar:
            with gzip.GzipFile('', 'wb', fileobj=out_gz, mtime=0) as gz_dest:
                gz_dest.writelines((line for line in src_tar))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--ncpus', type=int, help='Number of threads')
    parser.add_argument('--annotation_file', type=str,
                        help='Filename for annotation .gtf')
    parser.add_argument('--genome_file', type=str,
                        help='Filename for reference genome .fasta')
    parser.add_argument('--output_file', type=str, 
                        help='Basename for the output. Format is gzipped tar (tar.gz or .tgz are the recommended suffixes).')
    args = parser.parse_args()

    call_star(args.ncpus, args.annotation_file, args.genome_file,
              args.output_file)
