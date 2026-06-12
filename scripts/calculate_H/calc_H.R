## Compute H using one SNP per independent block

args=commandArgs(TRUE)

if(length(args)<5){stop("Rscript calc_H.R <prefix to plink files> <r prefix> <var prefix> <out prefix> <snp> <ids> <outfile>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(tidyverse)
}))

r_prefix = args[1]
FGr_file = args[2]
snpnum_file = args[3]
tvec_file = args[4]
out_file = args[5]
snp_file = args[6]

####################################
########## Functions ###############
####################################


calc_fhat <- function(dfMat, r ) {

  r <- r - mean(r)
  fhat_raw <- apply(dfMat, 1, sum)
  rTr <- as.numeric(t(as.matrix(r)) %*% as.matrix(r))
  fhat <- fhat_raw / c(rTr)

  return(fhat)
}


calc_sigma2_f <- function(fhat, M) {

  fhat <- fhat - mean(fhat)
  numerator <- as.numeric(t(fhat) %*% fhat)
  out <- numerator / (M-1)

  return(out)
}


calc_sigma2_r <-function(r, L) {

  r <- r - mean(r)
  rTr <- as.numeric(t(as.matrix(r)) %*% as.matrix(r))
  out <- rTr / (L - 1)

  return(out)
}


####################################
########## Main  ###################
####################################



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
print(paste0("Number of PC SNPs ", nrow(dfSNP)))

# Combine SNP and R files
dfALL <- inner_join(df, dfSNP) %>% drop_na()
print(paste0("There are ", nrow(dfALL), " SNPs in all the R files combined with the pruned SNPs"))
L <- nrow(dfALL)
print(L)

# Read in Fmat
dfMat <- as.matrix(fread(FGr_file, drop=1))
M <- nrow(dfMat)
print(M)

### Calculate H Real
fhat <- calc_fhat(dfMat, dfALL$r)
sigma2F <- as.numeric(calc_sigma2_f(fhat, M))
print(paste0("Sigma2F is ", sigma2F))
sigma2r <- as.numeric(calc_sigma2_r(dfALL$r, L))
print(paste0("Sigma2r is ", sigma2r))
H <- as.numeric(sigma2F * sigma2r)
print(paste0("H is ", H))

# Read in block info
dfSNPs <- fread(snpnum_file)
numBlocks <- as.numeric(ncol(dfMat))

#### Calculate SE via block jackknife
allHs <- rep(NA, numBlocks)
for (i in 1:numBlocks) {

  print(i)

  # Block num
  blockNum <- as.numeric(dfSNPs[i,1])

  # Calc H
  dfR_i <- dfALL %>% filter(block == blockNum)
  mi <- nrow(dfR_i)
  dfR_not_i <- dfALL %>% filter(block != blockNum)
  fhat_i <- calc_fhat(dfMat[,-i],dfR_not_i$r)
  sigma2F_i <- as.numeric(calc_sigma2_f(fhat_i, M))
  sigma2r_i <- as.numeric(calc_sigma2_r(dfR_not_i$r,L-mi))
  Hi <- as.numeric(sigma2F_i * sigma2r_i)
  allHs[i] <- as.numeric(((L - mi)/mi) * (H - Hi)^2)

}

# Calculate SE
varH <- mean(allHs)

# P-value from jacknife
pvalNorm <- pnorm(H, mean = 1/L, sd=sqrt(varH), lower.tail = FALSE)

# Get N
dfTvec <- fread(tvec_file)
N <- nrow(dfTvec)

# Construct output
dfOut <- data.frame(H = H, L = L, varH = varH, pvalNorm = pvalNorm, sigma2r=sigma2r, sigma2F=sigma2F, N=N)
fwrite(dfOut, out_file, quote = F, row.names = F, sep = "\t")




