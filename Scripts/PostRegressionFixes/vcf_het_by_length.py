#!/usr/bin/env python
"""
Calculate length-based heterozygosity for STR VCF
"""

import argparse
import numpy as np
import scipy.stats
import sys
import vcf


def GetLengthHet(record):
    len_to_count = {}
    for sample in record:
        if sample["GT"] is None or sample["GT"] == ".": continue
        alleles = sample.gt_bases.split(sample.gt_phase_char())
        l1 = len(alleles[0])
        l2 = len(alleles[1])
        len_to_count[l1] = len_to_count.get(l1, 0) + 1
        len_to_count[l2] = len_to_count.get(l2, 0) + 1
    total = sum(len_to_count.values())
    x = 0
    for al in len_to_count.keys():
        x += (len_to_count[al]*1.0/total)**2
    return 1-x, len(len_to_count.keys())

def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--vcf", help="VCF file", type=str, required=True)
    args = parser.parse_args()

    reader = vcf.Reader(open(args.vcf, "rb"))
    for record in reader:
        num_alleles = len([item for item in record.alleles if item is not None])
        lenhet, num_len_alleles = GetLengthHet(record)
        sys.stdout.write("\t".join(map(str, [record.CHROM, record.POS, record.INFO["END"], \
                                             num_len_alleles, num_alleles,
                                             lenhet, record.heterozygosity]))+"\n")
        sys.stdout.flush()

if __name__ == "__main__":
    main()