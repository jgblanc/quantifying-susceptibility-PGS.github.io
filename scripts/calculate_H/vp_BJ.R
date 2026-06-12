## Compute H using one SNP per independent block

args=commandArgs(TRUE)

if(length(args)<6){stop("Rscript calc_H.R <prefix to plink files> <r prefix> <var prefix> <out prefix> <snp> <ids> <outfile>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(tidyverse)
}))

plink_prefix = args[1]
r_prefix = args[2]
out_prefix = args[3]
snp_file = args[4]
id_file = args[5]
out_file = args[6]
nsnp = as.numeric(args[7])


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
numBlocks <- length(unique(dfSNP$block))
nsnp_per_block <- floor(nsnp / numBlocks)
dftmp1 <- inner_join(dfSNP, df, by = "ID") %>% drop_na()
dftmp2 <- dftmp1 %>%
  group_by(block) %>%
  mutate(n_snp_block = n()) %>%
  sample_n(min(nsnp_per_block, n_snp_block)) %>%
  ungroup() %>% select(-n_snp_block)
makeup <- nsnp - nrow(dftmp2)
print(makeup)
dftmp3 <- dftmp1 %>% filter(!ID %in% dftmp2$ID) %>% sample_n(makeup)
dfALL <- rbind(dftmp2, dftmp3)
print(paste0("There are ", nrow(dfALL), " SNPs in all the R files combined with the pruned SNPs"))
L <- nrow(dfALL)
print(L)

# Standardize r values and divide by GWAS variance
dfALL$r <- dfALL$r / sd(dfALL$r)
print(paste0("The variance of r is ", var(dfALL$r)))
print(paste0("The mean of r is ", mean(dfALL$r)))

# Set up data frame to collect SNP numbers
numBlocks <- length(unique(dfALL$block))
dfSNPs <- as.data.frame(matrix(NA, ncol = 2, nrow = numBlocks))
colnames(dfSNPs) <- c("Block", "nSNP")
print(dim(dfSNPs))
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
  dfR_tmp <- dfALL %>% filter(block == blockNum) %>% select("ID", "ALT", "r")
  fwrite(dfR_tmp, tmp_r_name, quote = F, row.names = F, sep = "\t")

  # Save number of SNPs
  nsnp_in_block <- nrow(dfR_tmp)
  dfSNPs[i,1] <- blockNum
  dfSNPs[i,2] <- nrow(dfR_tmp)

  # Run plink
  system(plink_cmd)

  # Read in plink output
  df <- fread(paste0(out_prefix, ".sscore"))
  rawFGr <- as.matrix(df[,3])
  dfFGr_mat[,i] <- rawFGr

}

# Remove tmp files
rm_cmd <- paste0("rm ", out_prefix, ".*")
system(rm_cmd)

# Calculate FGr
FGr_raw <- apply(dfFGr_mat, 1, sum)
print(paste0("The raw var is ", var(FGr_raw)))

# Scale by 1/sqrt(L-1)
FGr <- FGr_raw * (1/(sqrt(L-1)))
print(paste0("The scaled var is ", var(FGr)))

# Calculate H
H <- (1/(M * (L-1))) * (t(FGr) %*% FGr)
print(paste0("H is ", H))
print(paste0("1/L is ", 1/L))

# Compute SE for H
allHs <- rep(NA, numBlocks)
allHi <- rep(NA, numBlocks)
for (i in 1:numBlocks) {

  print(paste0("This is rep number ",i)) 
  mi <- as.numeric(dfSNPs[i,2])
  FGri <- (FGr_raw - dfFGr_mat[,i]) * (1/sqrt(L-mi-1))
  Hi <- sum(FGri^2) * (1/M) * (1/(L-mi-1))
  allHs[i] <- ((L - mi)/mi) * (H - Hi)^2
  allHi[i] <- Hi 

}
varH <- mean(allHs)
meanH <- mean(allHi)

# P-value from sims
pvalSim <- pnorm(H ,mean = meanH, sd = sqrt(varH), lower.tail = FALSE)
pvalNorm <- pnorm(H, mean = 1/L, sd=sqrt(varH), lower.tail = FALSE)

# Construct output
dfOut <- data.frame(H = H, L = L, meanH = meanH, varH = varH, pvalNorm = pvalNorm, pvalSim = pvalSim)
fwrite(dfOut, out_file, quote = F, row.names = F, sep = "\t")



