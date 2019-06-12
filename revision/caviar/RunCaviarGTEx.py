#!/usr/bin/env python3

# Example
# ./RunCaviarGTEx.py --linreg_snp /storage/szfeupe/Runs/650GTEx_estr/Analysis_by_Tissue/Nerve-Tibial/SNP_Analysis/Lin_Reg_Out --linreg_str /storage/szfeupe/Runs/650GTEx_estr/Analysis_by_Tissue/Nerve-Tibial/Lin_Reg_Out --samples /storage/mgymrek/gtex-estrs/revision/caviar/samples/Nerve-Tibial.samples --strgt /storage/mgymrek/gtex-estrs/revision/caviar/genotypes/GTExNormalizedSTRGenotypes.table.gz --snpgt /storage/mgymrek/gtex-estrs/revision/caviar/genotypes/GTExNormalizedSNPGenotypes_chr21.table.gz --out test --genes ENSG00000160213.5 --tmpdir test/

# Notes:
# Make STR and SNP genotypes indexed, maybe all in one file to keep sample order same?

import argparse
import gzip
import numpy as np
import os
import pandas as pd
from subprocess import Popen, PIPE, DEVNULL
import sys
import tabix

def PROGRESS(msg, printit=True):
    if printit: # false for some messages when not in debug mode
        sys.stderr.write("%s\n"%msg.strip())

def LoadSamples(samplesfile):
    if samplesfile is None: return []
    return [item.strip() for item in open(samplesfile, "r").readlines()]

# Quickly load for small gene sets
def LoadReg(linreg, genes, prefix="", tempdir="/tmp"):
    # Write gene list to a file
    with open(os.path.join(tempdir, "genelist.txt"), "w") as f:
        for gene in genes: f.write(gene+"\n")
    # Grep for the gene or gene list from the linreg file
    newlinreg = os.path.join(tempdir, os.path.basename(linreg))
    if os.path.exists(newlinreg):
        PROGRESS("Warning: Intermediate linreg file %s exists. Overwriting\n"%os.path.join(tempdir, os.path.basename(linreg)))
    cmd1 = "head -n 1 %s > %s"%(linreg, newlinreg)
    cmd2 = "grep -f %s %s >> %s"%(os.path.join(tempdir, "genelist.txt"), \
                                linreg, newlinreg)
    os.system(cmd1+";"+cmd2)
    # Open the new reduced regression file
    reg = pd.read_csv(newlinreg, sep="\t", \
                      usecols=["gene","chrom","str.start","beta","beta.se","p.wald"])
    reg["Z"] = reg.apply(lambda x: x["beta"]/x["beta.se"], 1)
    reg["ID"] = reg.apply(lambda x: prefix+x["chrom"]+":"+str(x["str.start"]), 1)
    return reg

# Write LDFILE, ZFILE
def GenerateCAVIARFiles(gene, samples, strreg, snpreg, strgt, snpgt, \
                        use_topn_strs, use_topn_snps, \
                        str_gt_ind, snp_gt_ind, \
                        tmpdir):
    if not os.path.exists(os.path.join(tmpdir, gene)): os.mkdir(os.path.join(tmpdir, gene))
    strdata = strreg[strreg["gene"]==gene].sort_values("p.wald").head(use_topn_strs).sort_values("str.start")
    snpdata = snpreg[snpreg["gene"]==gene].sort_values("p.wald").head(use_topn_snps).sort_values("str.start")
    # 1. Get ZFILE
    zfile = os.path.join(tmpdir, gene, "ZFILE")
    strdata[["ID", "Z"]].to_csv(open(zfile, "w"), header=None, index=False, sep="\t")
    snpdata[["ID", "Z"]].to_csv(open(zfile, "a"), header=None, index=False, sep="\t")
    # 2. Get LDFILE for only that set of variants
    str_genotypes = LoadGenotypes(strgt, str_gt_ind, strdata)
    snp_genotypes = LoadGenotypes(snpgt, snp_gt_ind, snpdata)
    all_genotypes = pd.DataFrame(str_genotypes + snp_genotypes)
    ldmatrix = np.square(all_genotypes.transpose().corr())
    ldfile = os.path.join(tmpdir, gene, "LDFILE")
    ldmatrix.to_csv(ldfile, header=None, index=False, sep="\t")

def GetFloat(value):
    if value == "None": return np.nan
    else: return float(value)

def LoadGenotypes(gtfile, gtind, regdata):
    chrom = regdata["chrom"].values[0]
    start = min(regdata["str.start"])
    end = max(regdata["str.start"])
    positions = list(regdata["str.start"])
    loaded_positions = []
    tb = tabix.open(gtfile)
    records = tb.query(chrom, start-1, end+1)
    data = []
    for record in records:
        pos = int(record[1])
        if pos not in positions: continue
        loaded_positions.append(pos)
        data.append([GetFloat(record[i+2]) for i in gtind]) # first two cols are chrom, start
    assert(len(positions)==len(loaded_positions))
    assert([positions[i]==loaded_positions[i] for i in range(len(positions))])
    return data

def GetGenotypeIndices(strgtfile, snpgtfile, samples):
    str_samples = [item.decode('UTF-8') for item in (gzip.open(strgtfile, "r").readline().strip().split()[2:])]
    snp_samples = [item.decode('UTF-8') for item in gzip.open(snpgtfile, "r").readline().strip().split()[2:]]
    use_samples = (set(str_samples).intersection(set(snp_samples))).intersection(samples)
    str_ind = [str_samples.index(item) for item in use_samples]
    snp_ind = [snp_samples.index(item) for item in use_samples]
    return str_ind, snp_ind, use_samples

# Run CAVIAR using LDFILE and ZFILE in tmp/
# Write output to tmp/
def RunCAVIAR(gene, tmpdir, numcausal):
    zfile = os.path.join(tmpdir, gene, "ZFILE")
    ldfile = os.path.join(tmpdir, gene, "LDFILE")
    outfile = os.path.join(tmpdir, gene, "CAVIAR")
    cmd = "CAVIAR -o %s -l %s -z %s -c %s"%(outfile, ldfile, zfile, numcausal)
    p = Popen(cmd, shell=True, stdout=DEVNULL, stderr=DEVNULL)
    output = p.communicate()[0]
    if p.returncode != 0:
        PROGRESS("CAVIAR on %s failed"%gene)
        return False
    return True

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
    parser.add_argument("--use-topn-strs", help="Use top n STRs (by p-value)", type=int, default=1)
    parser.add_argument("--use-topn-snps", help="Use top n SNPs (by p-value)", type=int, default=1000000)
    parser.add_argument("--num-causal", help="Number of causal variants to consider", type=int, default=3)
    parser.add_argument("--tmpdir", help="Use this directory for temporary files", type=str, default="/tmp")
    args = parser.parse_args()

    # Get list of genes to process
    if args.genes == "all": genes = set(strreg["gene"])
    else: genes = set(args.genes.split(","))

    if not args.precomputed:
        # Load regression results
        PROGRESS("\nLoad snps regression")
        snpreg = LoadReg(args.linreg_snp, genes, prefix="SNP_", tempdir=args.tmpdir)
        PROGRESS("\nLoad strs regression")
        strreg = LoadReg(args.linreg_str, genes, prefix="STR_", tempdir=args.tmpdir)

        # Get list of samples to process
        samples = LoadSamples(args.samples)

        # Get sample indices for genotype data
        str_gt_ind, snp_gt_ind, samples = GetGenotypeIndices(args.strgt, args.snpgt, samples)

    # For each gene:
    # 1. Get intermediate files
    # 2. Run CAVIAR
    # 3. Generate output
    for gene in genes:
        PROGRESS("Processing gene %s"%gene)
        if not args.precomputed: # Store in args.tmpdir/gene/ LDFILE, ZFILE
             GenerateCAVIARFiles(gene, samples, strreg, snpreg, args.strgt, args.snpgt, \
                                 args.use_topn_strs, args.use_topn_snps, \
                                 str_gt_ind, snp_gt_ind, \
                                 args.tmpdir)
        if not RunCAVIAR(gene, args.tmpdir, args.num_causal): continue
        WriteOutput(gene, args.tmpdir, args.out)
