## Find the estimation error of FGr using a Jackknife approach

args=commandArgs(TRUE)

if(length(args)<2){stop("Rscript compute_error_jacknife.R <prefix to Tm chromosomes> <outfile>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
  library(tidyverse)
}))

infileVar = args[1]
ldFile = args[2]
outfile = args[3]

ld <- fread(ldFile)

df <- fread(infileVar)
df <- df %>% separate(ID, c("CHR", "POS"), remove = FALSE)
print(head(df))

# Assign SNPs to blocks
assign_SNP_to_block <- function(CHR, BP, block = ld) {

  # Filter blocks based on snp
  block_chr <- block %>% filter(chr == CHR)
  first_start <- as.numeric(block_chr[1, "start"])
  block_bp <- block_chr %>% filter( (start < BP & stop >= BP) | BP == first_start)

  # Assign
  block_num <- as.numeric(block_bp[,"block_number"])
  return(block_num)
}

# Add block info - takes a while
df <- df %>%
  mutate(block = apply(., MARGIN = 1, FUN = function(params)assign_SNP_to_block(as.numeric(params[3]), as.numeric(params[4])))) %>%
  drop_na() %>% select("#CHROM","ID",	"REF","ALT","BETA", "ID", "REF", "block")
print(paste0("Now df blocks has", nrow(df), " rows"))
print(head(df))


fwrite(df, outfile, row.names = F, col.names = T, quote = F, sep = "\t")







