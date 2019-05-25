#!/usr/bin/env python3

# Example
# ./RunCaviarGTEx.py --linreg_snp /storage/szfeupe/Runs/650GTEx_estr/Analysis_by_Tissue/Nerve-Tibial/SNP_Analysis/Lin_Reg_Out --linreg_str /storage/szfeupe/Runs/650GTEx_estr/Analysis_by_Tissue/Nerve-Tibial/Lin_Reg_Out --samples /storage/mgymrek/gtex-estrs/revision/caviar/samples/Nerve-Tibial.samples --strgt /storage/mgymrek/gtex-estrs/revision/caviar/genotypes/gtex_strgt.tab --snpgt /storage/mgymrek/gtex-estrs/revision/caviar/genotype/gtex_snpgt.tab --out test --genes ENSG00000160213.5 --tmpdir test/

# Notes:
# Make STR and SNP genotypes indexed, maybe all in one file to keep sample order same?

import argparse
import os
import pandas as pd
import sys

def PROGRESS(msg, printit=True):
    if printit: # false for some messages when not in debug mode
        sys.stderr.write("%s\n"%msg.strip())

def LoadReg(linreg, genes, prefix=""):
    reg = pd.read_csv(linreg, sep="\t", usecols=["gene","chrom","str.start","beta","beta.se"])
    reg = reg[reg["gene"].apply(lambda x: x in genes)]
    reg["Z"] = reg.apply(lambda x: x["beta"]/x["beta.se"], 1)
    reg["ID"] = reg.apply(lambda x: prefix+x["chrom"]+":"+x["str.start"], 1)
    return reg

# Write LDFILE, ZFILE
def GenerateCAVIARFiles(gene, samples, strreg, snpreg, strgt, snpgt, tmpdir):
    if not os.path.exists(os.path.join(tmpdir, gene)): os.mkdir(os.path.join(tmpdir, gene))
    # 1. Get ZFILE
    zfile = (os.path.join(tmpdir, gene, "ZFILE"))
    strreg[strreg["gene"]==gene][["ID", "Z"]].to_csv(open(zfile, "w"), header=None, index=False)
    snpreg[snpreg["gene"]==gene][["ID", "Z"]].to_csv(open(zfile, "a"), header=None, index=False)
    # 2. Get LDFILE for only that set of variants
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

    # Get list of genes to process
    if args.genes == "all": genes = set(strreg["gene"])
    else: genes = set(args.genes.split(","))

    # Load regression results
    if not args.precomputed:
        PROGRESS("\nLoad snps regression")
        snpreg = LoadReg(args.linreg_snp, genes, prefix="SNP_")
        PROGRESS("\nLoad strs regression")
        strreg = LoadReg(args.linreg_str, genes, prefix="STR_")

        # Get list of samples to process
        samples = [] # TODO get from args

    # For each gene:
    # 1. Get intermediate files
    # 2. Run CAVIAR
    # 3. Generate output
    for gene in genes:
        PROGRESS("Processing gene %s"%gene)
        if not args.precomputed: # Store in args.tmpdir/gene/ LDFILE, ZFILE
             GenerateCAVIARFiles(gene, samples, strreg, snpreg, args.strgt, args.snpgt, args.tmpdir)
        RunCAVIAR(gene, args.tmp)
        WriteOutput(gene, args.tmp, args.out)
