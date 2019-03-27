"""
Calculates and formats star qc metrics
"""

__author__ = 'Otto Jolanki'
__license__ = 'MIT'

import argparse
import json
import logging
import pandas as pd
from qc_utils import QCMetric, QCMetricRecord
from qc_utils.parsers import parse_starlog

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
filehandler = logging.FileHandler('star_qc.log')
filehandler.setLevel(logging.DEBUG)
consolehandler = logging.StreamHandler()
consolehandler.setLevel(logging.INFO)
formatter = logging.Formatter(
    '%(asctime)s | %(levelname)s | %(name)s: %(message)s')
filehandler.setFormatter(formatter)
consolehandler.setFormatter(formatter)
logger.addHandler(consolehandler)
logger.addHandler(filehandler)


def main(args):
    logger.info('Reading input tsv: %s' % args.quants)
    quants_tsv = pd.read_csv(args.quants, sep='\t', header=None, skiprows=4)
    # calculate number of mirnas expressed at cpm>2
    per_million = quants_tsv[1].sum() / 1000000
    quants_tsv['cpm'] = quants_tsv[1] / per_million
    cpm_gte2 = sum(quants_tsv['cpm'] >= 2)
    star_qc_record = QCMetricRecord()
    cpm_metric = QCMetric('expressed_mirnas',
                          {'expressed_mirnas': cpm_gte2})
    # get metrics from star log
    star_qc = QCMetric('star_qc_metric', args.star_log, parse_starlog)
    star_qc_record.add_all([cpm_metric, star_qc])
    # calculate number of reads (unique + multimapping)
    reads_mapped = int(star_qc.content['Uniquely mapped reads number']) + int(
        star_qc.content['Number of reads mapped to multiple loci'])
    reads_mapped_qc = QCMetric('aligned_reads', {'aligned_reads': reads_mapped})
    star_qc_record.add(reads_mapped_qc)
    logger.info('Writing output json %s' % args.output_filename)
    with open(args.output_filename, 'w') as fp:
        json.dump(star_qc_record.to_ordered_dict(), fp)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--quants',
        type=str,
        required=True,
        help='Tab separated file containing quantitations.')
    parser.add_argument(
        '--star_log',
        type=str,
        required=True,
        help='STAR Log.final.out ending log file.')
    parser.add_argument(
        '--output_filename',
        type=str,
        required=True,
        help='Name of the final qc json file')
    args = parser.parse_args()
    main(args)
