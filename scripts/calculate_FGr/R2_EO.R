## Compute FGr using one PC SNPs

args=commandArgs(TRUE)

if(length(args)<4){stop("Rscript calc_R2_EO.R")}


library(data.table)
library(tidyverse)
library(matrixStats)
library(future.apply)  # parallel bootstrap
plan(multisession)
plan(multicore, workers = 8)
options(future.globals.maxSize = 2 * 1024^3)

pca_file = args[1]
fgr_file = args[2]
snp_file = args[3]
out_file = args[4]
r_prefix = args[5]
pc_snp_file = args[6]

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


# Read inputs
PCs <- fread(pca_file)
PC_nums <- scale(as.matrix(PCs[, 3:ncol(PCs)]))
print(apply(PC_nums, 2, var))
dfSNP <- fread(snp_file)

# Get even chr
chrFGr <- seq(2, 22, 2)
dfSNP_filter <- dfSNP %>% filter(CHR %in% chrFGr)
blockIndex <- which(dfSNP$CHR %in% chrFGr)
numBlocks <- length(blockIndex)
print(paste0("There are ", numBlocks, "on even Chromosomes"))

# Read in R values

## First Chr
r_file_name <- paste0(r_prefix, 2, ".rvec")
df <- fread(r_file_name)
df$CHR <- 2

## All other Chrs
for (i in seq(4, 22, 2)) {

  # Read in R file
  r_file_name <- paste0(r_prefix, i, ".rvec")
  tmp <- fread(r_file_name)
  tmp$CHR <- i

  # Combine
  df <- rbind(df, tmp)
}
print(paste0("There are ", nrow(df), " SNPs in all the R files"))

# Read in SNP file
dfSNP <- fread(pc_snp_file)  %>% select("ID", "block")
print(paste0("Number of PC SNPs ", nrow(dfSNP)))

# Combine SNP and R files
dfALL <- inner_join(df, dfSNP) %>% drop_na() %>% filter(CHR %in% chrFGr)
print(paste0("There are ", nrow(dfALL), " SNPs in all the R files combined with the pruned SNPs"))
L <- nrow(dfALL)
print(L)

# Read in and compute FGr
dfMat <- as.matrix(fread(fgr_file, drop=1))[, blockIndex]
fhat <- calc_fhat(dfMat, dfALL$r)
M <- length(fhat)

### Calculate H Real
sigma2F <- as.numeric(calc_sigma2_f(fhat, M))
print(paste0("Sigma2F is ", sigma2F))
sigma2r <- as.numeric(calc_sigma2_r(dfALL$r, L))
print(paste0("Sigma2r is ", sigma2r))
H <- as.numeric(sigma2F * sigma2r)
print(paste0("H is ", H))

#### Calculate SE via block jackknife
allHs <- rep(NA, numBlocks)
for (i in 1:numBlocks) {


  # Block num
  blockNum <- as.numeric(dfSNP_filter[i,1])
  print(paste0("This is rep number ",i))
  print(paste0("This is block number ",blockNum))

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


### Calculate signal

# Leave one out f;s
jckFGr <- matrix(NA, nrow = M, ncol = numBlocks)
for (i in 1:numBlocks) {

  # Block num
  blockNum <- as.numeric(dfSNP_filter[i,1])
  print(paste0("This is rep number ",i))
  print(paste0("This is block number ",blockNum))
  dfR_not_i <- dfALL %>% filter(block != blockNum)
  dfR_i <- dfALL %>% filter(block == blockNum)
  mi <- nrow(dfR_i)
  fhat_i <- calc_fhat(dfMat[,-i],dfR_not_i$r)
  jckFGr[,i] <- ((L - mi)/mi) * (fhat - fhat_i)^2

}

# Compute numerator
tmp <- apply(jckFGr, 1, mean)
numerator <- mean(tmp)
print(paste0("The numerator is ",numerator))


# Compute Denominator
varFGr <- var(fhat)
print(paste0("The denominator is ",varFGr))

# Find Error
error <- numerator / varFGr

# Final signal
signal <- 1 - error

# Construct output
dfOut <- matrix(NA, nrow = ncol(PC_nums), ncol = 8)
colnames(dfOut) <- c("PC","H", "varH", "pvalNorm", "signal", "omega","R2", "Ratio")

# Loop through PCs
for (i in seq_len(ncol(PC_nums))) {

  # Get Single PC stats
  mod <- lm(fhat ~ PC_nums[,i])
  w <- mod$coefficients[2]

  # Fit all PCs
  mod  <- lm(fhat ~ PC_nums[,1:i])
  R2 <- summary(mod)$r.squared

  # Get Ratio
  Ratio <- R2 / signal

  dfOut[i,] <- c(i, H, varH, pvalNorm, signal, w, R2, Ratio)
  cat("Finished PC", i, "\n")
}


# Write output
fwrite(as.data.table(dfOut), out_file, sep = "\t", quote = FALSE)






