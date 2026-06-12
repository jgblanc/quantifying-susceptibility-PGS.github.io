## Format phenotype file for REGENIE

args=commandArgs(TRUE)

if(length(args)<2){stop("Rscript calc_R2_EO.R")}


library(data.table)
library(tidyverse)

ss_file = args[1]
snp_file = args[2]
pvar_file = args[3]
out_file = args[4]


# Read in SNPs
dfSNPs <- fread(snp_file)

# Read in pvar
dfPvar <- fread(pvar_file)

# Join files
dfSNPs <- inner_join(dfSNPs, dfPvar)

# Read in SS
dfSS <- fread(ss_file)
dfSS <- dfSS %>% select("CHROM", "GENPOS", "ID", "ALLELE0", "ALLELE1", "BETA", "LOG10P")
colnames(dfSS) <- c("#CHROM", "POS", "rsID", "REF", "ALT", "BETA", "LOG10_P")

# Join DFs
dfOut <- inner_join(dfSNPs, dfSS)
dfOut <- dfOut %>% select("#CHROM", "POS","ID", "REF", "ALT", "BETA", "LOG10_P")


# Write output
fwrite(dfOut, out_file,row.names = FALSE, sep = "\t", quote = FALSE,na = "NA")
