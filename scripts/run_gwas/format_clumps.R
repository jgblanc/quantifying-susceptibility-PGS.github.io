## Format phenotype file for REGENIE

args=commandArgs(TRUE)

if(length(args)<2){stop("Rscript calc_R2_EO.R")}


library(data.table)
library(tidyverse)

ss_file = args[1]
clump_file = args[2]
out_file = args[3]

# Read in SS
dfSNPs <- fread(ss_file)

# Read clumps
dfClump <- fread(clump_file)

# Join files
if (nrow(dfSNPs) > 0) {
   dfOut <- inner_join(dfSNPs, dfClump)
   dfOut <- dfOut %>% select("#CHROM", "ID", "ALLELE0", "ALLELE1", "BETA")
   colnames(dfOut) <- c("#CHROM", "ID", "REF", "ALT", "BETA")
   fwrite(dfOut, out_file,row.names = FALSE, sep = "\t", quote = FALSE,na = "NA")
}
