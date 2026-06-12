## Compute Add together all LD blocks

args=commandArgs(TRUE)

if(length(args)<3){stop("Rscript calc_FGr.R <prefix to plink files> <r prefix> <var prefix> <out prefix> <snp> <ids> <outfile>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(tidyverse)
}))

FGr_file = args[1]
out_Fhat = args[2]
r_prefix = args[3]
id_file = args[4]
snp_file =args[5]

# Read in Fmat
dfMat <- as.matrix(fread(FGr_file, drop = 1))
fhat_raw <- apply(dfMat, 1, sum)
print(length(fhat_raw))

# Read in r values

## First Chr
r_file_name <- paste0(r_prefix, 1, ".rvec")
df <- fread(r_file_name)
df$CHR <- 1

## All other Chrs
for (i in 2:22) {

  # Read in R file
  r_file_name <- paste0(r_prefix, i, ".rvec")
  tmp <- fread(r_file_name)
  tmp$CHR <- i

  # Combine
  df <- rbind(df, tmp)
}
print(paste0("There are ", nrow(df), " SNPs in all the R files"))

# Read in SNP file
dfSNP <- fread(snp_file) %>% select("ID", "block")
print(paste0("Number of PC SNPs ", nrow(dfSNP)))

# Combine SNP and R files
dfALL <- inner_join(df, dfSNP) %>% drop_na()
print(paste0("There are ", nrow(dfALL), " SNPs in all the R files combined with the pruned SNPs"))

## calculate r^\top t
r <- dfALL$r
rTr <- t(as.matrix(r)) %*% as.matrix(r)
print(rTr)

## calculate \hat{f}
fhat<- fhat_raw / c(rTr)

## format output
dfOut <- fread(id_file) %>% select("#FID", "IID")
dfOut$fhat <- fhat

# FGR output
fwrite(dfOut, out_Fhat, quote = F, row.names = F, sep = "\t")
