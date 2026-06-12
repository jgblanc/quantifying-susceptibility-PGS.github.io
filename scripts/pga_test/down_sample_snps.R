# Script to run selection test with jacknife error

args=commandArgs(TRUE)

if(length(args)<2){stop("Provide <formatted betas> <output>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(tidyverse)
}))

outfile=args[1]
target_snps=args[2]

# Read in target SNPs
df <- fread(target_snps)
print(head(df))

# Loop through Chrs
nChr = length(args) - 2

# Do first chr
dfR <- fread(args[3])
print(head(dfR))
dfOut <- df %>%
  filter(`#CHROM` == 1) %>%
  left_join(dfR %>% select(ID, BETA_regenie = BETA), by = "ID") %>%
  mutate(BETA = coalesce(BETA_regenie, BETA)) %>%
  select(-BETA_regenie)
print(head(dfOut))


# Loop through CHRs
for (i in 1:(nChr-1)) {
  idx <- 3+i
  chr_idx <- i+1
  dfR <- fread(args[idx])

  tmp <- df %>%
    filter(`#CHROM` == chr_idx) %>%
    left_join(dfR %>% select(ID, BETA_regenie = BETA), by = "ID") %>%
    mutate(BETA = coalesce(BETA_regenie, BETA)) %>%
    select(-BETA_regenie)
  dfOut <- rbind(dfOut, tmp)

}


# Save output
fwrite(dfOut, outfile,col.names=T,row.names=F,quote=F,sep="\t")
