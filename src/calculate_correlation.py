"""
Calculates spearman correlation between the quants in replicates.
"""

__author__ = 'Otto Jolanki'
__license__ = 'MIT'

import argparse
import json
import pandas as pd
from qc_utils import QCMetric, QCMetricRecord


def main(args):
    quant1 = pd.read_csv(args.quants[0], sep='\t', header=None, skiprows=4)
    quant2 = pd.read_csv(args.quants[1], sep='\t', header=None, skiprows=4)
    spearman_correlation = quant1[1].corr(quant2[1], method='spearman')
    qc_record = QCMetricRecord()
    spearman_metric = QCMetric('spearman_correlation', {'spearman_correlation': spearman_correlation})
    qc_record.add(spearman_metric)
    with open(args.output_filename, 'w') as fp:
        json.dump(qc_record.to_ordered_dict(), fp)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--quants', nargs=2, help='Quantitation tsvs to correlate.')
    parser.add_argument('--output_filename', type=str)
    args = parser.parse_args()
    main(args)
