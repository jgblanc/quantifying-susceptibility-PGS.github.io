## Compute FGr using one PC SNPs

args=commandArgs(TRUE)

if(length(args)<7){stop("Rscript calc_FGr.R <prefix to plink files> <r prefix> <var prefix> <out prefix> <snp> <ids> <outfile>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(tidyverse)
}))

plink_prefix = args[1]
r_prefix = args[2]
out_prefix = args[3]
snp_file = args[4]
id_file = args[5]
out_file_FGr = args[6]
out_file_SNP = args[7]

# Read in IDs
dfIDs <- fread(id_file)

# Set data collectors
M <- nrow(dfIDs)

# Read in all values of r

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


# Combine SNP and R files
dfALL <- inner_join(df, dfSNP) %>% drop_na()
print(paste0("There are ", nrow(dfALL), " SNPs in all the R files combined with the pruned SNPs"))
L <- nrow(dfALL)
print(paste0("There are ", L, " snps total"))


# Set up data frame to collect SNP numbers
numBlocks <- length(unique(dfALL$block))
dfSNPs <- as.data.frame(matrix(NA, ncol = 3, nrow = numBlocks))
colnames(dfSNPs) <- c("Block", "nSNP", "CHR")
print(paste0("The total number of blocks is ", numBlocks))
dfFGr_mat <- matrix(NA, nrow = M, ncol = numBlocks)

# Subset SNP IDs
dfSNP_tmp <- dfALL %>% select("ID")
snp_name <- paste0(out_prefix, ".snp")
fwrite(dfSNP_tmp, snp_name, quote = F, row.names = F, sep = "\t")

# Get freq file
plink_cmd <- paste0("plink2 --pfile ", plink_prefix, " --keep ", id_file, " --extract ", snp_name ," --threads 8 ",
                      " --freq --out ", out_prefix)
system(plink_cmd)

# Set up plink command
freq_file <- paste0(out_prefix, ".afreq")
tmp_r_name <- paste0(out_prefix, ".rvec")
plink_cmd <- paste0("plink2 --pfile ", plink_prefix, " --keep ", id_file, " --extract ", snp_name ," --threads 8 --read-freq ", freq_file,
                      " --score ", tmp_r_name, " header-read variance-standardize cols=dosagesum,scoresums --out ", out_prefix)

# Individually score each of the 581 blocks
for (i in 1:numBlocks) {

  # Block num
  blockNum <- unique(dfALL$block)[i]
  print(blockNum)

  # Subset Rs and save
  dfR_tmp1 <- dfALL %>% filter(block == blockNum)
  chr <- dfR_tmp1$CHR[1]
  print(chr)
  dfR_tmp <- dfR_tmp1 %>% select("ID", "ALT", "r")
  fwrite(dfR_tmp, tmp_r_name, quote = F, row.names = F, sep = "\t")

  # Save number of SNPs
  nsnp_in_block <- nrow(dfR_tmp)
  dfSNPs[i,1] <- blockNum
  dfSNPs[i,2] <- nrow(dfR_tmp)
  dfSNPs[i,3] <- chr

  # Run plink
  system(plink_cmd)

  # Read in plink output
  df <- fread(paste0(out_prefix, ".sscore"))
  rawFGr <- as.matrix(df[,3])
  print(length(rawFGr))
  dfFGr_mat[,i] <- rawFGr

}

# Remove tmp files
rm_cmd <- paste0("rm ", out_prefix, ".*")
system(rm_cmd)

# Save Raw FGR
dfFGr <- as.data.frame(dfFGr_mat)
dfFGr <- cbind(IID = df[[1]], dfFGr)
fwrite(dfFGr, out_file_FGr, quote = F, row.names = F, sep = "\t")

# Save SNP file
fwrite(dfSNPs, out_file_SNP, quote = F, row.names = F, sep = "\t")



