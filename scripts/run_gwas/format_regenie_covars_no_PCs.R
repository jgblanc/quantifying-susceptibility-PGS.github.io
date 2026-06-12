## Format phenotype file for REGENIE

args=commandArgs(TRUE)

if(length(args)<2){stop("Rscript calc_R2_EO.R")}


library(data.table)
library(tidyverse)

in_file = args[1]
out_file = args[2]


# Read in phenotypes
df <- fread(in_file)

# Add IID column
df <- df[, c(1, 1, 2:4)]

# Name columns
colnames(df) <- c("FID", "IID", "Age", "Sex", "Batch")
df <- df %>% mutate(Batch = replace(Batch, Batch < 0, 0)) %>%
  mutate(Batch = replace(Batch, Batch > 0, 1))


# Write output
fwrite(df, out_file,row.names = FALSE, sep = "\t", quote = FALSE,na = "NA")
