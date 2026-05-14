#!/bin/bash
pfile_path=$1
pheno_path=$2
outfile=$3
overlap_snps=$4

plink2 \
  --pfile $pfile_path \
  --extract $overlap_snps \
  --glm omit-ref allow-no-covars \
  --no-input-missing-phenotype \
  --pheno $pheno_path \
  --pheno-name tvec \
  --geno-counts \
  --out $outfile
