#!/usr/bin/env python3

# Example
# ./RunCaviarGTEx.py --linreg_snp /storage/szfeupe/Runs/650GTEx_estr/Analysis_by_Tissue/Nerve-Tibial/SNP_Analysis/Lin_Reg_Out --linreg_str /storage/szfeupe/Runs/650GTEx_estr/Analysis_by_Tissue/Nerve-Tibial/Lin_Reg_Out --samples /storage/mgymrek/gtex-estrs/revision/caviar/samples/Nerve-Tibial.samples --strgt /storage/mgymrek/gtex-estrs/revision/caviar/genotypes/gtex_strgt.tab --snpgt /storage/mgymrek/gtex-estrs/revision/caviar/genotype/gtex_snpgt.tab --out test --genes ENSG00000160213.5 --tmpdir test/

# Notes:
# Make STR and SNP genotypes indexed, maybe all in one file to keep sample order same?

import argparse
import pandas as pd
import sys

def PROGRESS(msg, printit=True):
    if printit: # false for some messages when not in debug mode
        sys.stderr.write("%s\n"%msg.strip())

# Write LDFILE, ZFILE
def GenerateCAVIARFiles(gene, samples, strreg, snpreg, strgt, snpgt, tmpdir):
    pass # TODO

# Run CAVIAR using LDFILE and ZFILE in tmp/
# Write output to tmp/
def RunCAVIAR(gene, tmp):
    pass # TODO

# Write output. Include info on failed genes
def WriteOutput(gene, tmp, out):
    pass # TODO

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run CAVIAR on GTEx data")
    parser.add_argument("--linreg_snp", help="File with snp linear regression output", type=str, required=True)
    parser.add_argument("--linreg_str", help="File with str linear regression output", type=str, required=True)
    parser.add_argument("--samples", help="File with samples to process for this tissue", type=str, required=True)
    parser.add_argument("--strgt", help="File with noramlized STR genotypes", type=str, required=True)
    parser.add_argument("--snpgt", help="File with normalized SNP genotypes", type=str, required=True)
    parser.add_argument("--out", help="Write results to this file", type=str, required=True)
    parser.add_argument("--genes", help="Only process these genes", type=str, default="all")
    parser.add_argument("--precomputed", help="Use precomputed input files for CAVIAR", action="store_true")
    parser.add_argument("--tmpdir", help="Use this directory for temporary files", type=str, default="/tmp")
    args = parser.parse_args()

    # Load regression results
    PROGRESS("\nLoad snps regression")
    snpreg = pd.read_csv(args.linreg_snp, sep="\t", usecols=["gene","chrom","str.start","beta","beta.se"])
    PROGRESS("\nLoad strs regression")
    strreg = pd.read_csv(args.linreg_str, sep="\t", usecols=["gene","chrom","str.start","beta","beta.se"])

    # Get list of samples to process
    samples = [] # TODO get from args

    # Get list of genes to process
    if args.genes == "all": genes = set(strreg["gene"])
    else: genes = set(args.genes.split(","))

    # For each gene:
    # 1. Get intermediate files
    # 2. Run CAVIAR
    # 3. Generate output
    for gene in genes:
        if not args.precomputed: # Store in args.tmpdir/gene/ LDFILE, ZFILE
            GenerateCAVIARFiles(gene, samples, strreg, snpreg, args.strgt, args.snpgt, args.tmpdir)
        RunCAVIAR(gene, args.tmp)
        WriteOutput(gene, args.tmp, args.out)
