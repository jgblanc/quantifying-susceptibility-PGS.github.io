## Compute H using one SNP per independent block

args=commandArgs(TRUE)

if(length(args)<5){stop("Rscript calc_H.R <prefix to plink files> <r prefix> <var prefix> <out prefix> <snp> <ids> <outfile>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(tidyverse)
}))

resid_prefix = args[1]
r_prefix = args[2]
snp_file = args[3]
id_file = args[4]
out_file_FGr = args[5]
out_file_SNP = args[6]
out_file = out_file_FGr

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
print(paste0("Number of PC SNPs ", nrow(dfSNP)))

# Combine SNP and R files
dfALL <- inner_join(df, dfSNP) %>% drop_na()
print(paste0("There are ", nrow(dfALL), " SNPs in all the R files combined with the pruned SNPs"))
L <- nrow(dfALL)
print(L)

# Standardize r values variance
dfALL$r <- dfALL$r / sd(dfALL$r)
print(paste0("The variance of r is ", var(dfALL$r)))
print(paste0("The mean of r is ", mean(dfALL$r)))

# Set up data frame to collect SNP numbers
numBlocks <- length(unique(dfALL$block))
dfSNPs <- as.data.frame(matrix(NA, ncol = 2, nrow = numBlocks))
colnames(dfSNPs) <- c("Block", "nSNP")
print(paste0("The total number of blocks is ", numBlocks))
dfFGr_mat <- matrix(NA, nrow = M, ncol = numBlocks)



counter <- 1
for (i in 1:22) {

  # Get r values on Chr i
  dfR_chr <- dfALL %>% filter(CHR == i)
  chr_blocks <- unique(dfR_chr$block)
  num_blocks_chr <- length(chr_blocks)
  print(paste0("The number of blocks on chr ", i, " is ", num_blocks_chr))

  # Get SNP ids on Chr i
  traw_file <- paste0(resid_prefix, i, "_residual.traw ")
  cmd_file<- paste0( "cut -f1 ",traw_file)
  dfSNPchr <- fread(cmd = cmd_file,header = TRUE)

  for (b in chr_blocks) {

    print(paste0("We are on block", b))
    # Select r values on block b
    dfR_block <- dfR_chr %>% filter(block == b)
    r_ids <- dfR_block$ID

    # Get the right row values
    row_ids <- which(dfSNPchr$SNP %in% r_ids)
    print(paste0("There are ", length(row_ids), " SNPs on the block" ))

    # read only correct row IDs
    row_nums_str <- paste(row_ids + 1, collapse = ",")  # +1 for header row
    cmd_block <- paste0("awk 'NR==1 || NR==", gsub(",", " || NR==", row_nums_str), "' ", traw_file)
    tmpfile <- paste0(out_file, "_tmp.txt")
    system(paste(cmd_block, ">", tmpfile))
    k <- floor(length(row_ids)/3)
    print(k)
    dfBlock1 <- fread(tmpfile,nrow=k, header=TRUE)
    dfBlock2 <- fread(tmpfile, skip=k+1, nrow=k, header=FALSE)
    dfBlock3 <- fread(tmpfile, skip=(2*k)+1, header=FALSE)
    setnames(dfBlock2, names(dfBlock1))
    setnames(dfBlock3, names(dfBlock1))
    dfBlock <- rbind(dfBlock1, dfBlock2)
    dfBlock <- rbind(dfBlock, dfBlock3)
    print(nrow(dfBlock))
    system(paste("rm", tmpfile))


    # do matrix multiplication
    matBlock <- t(as.matrix(dfBlock[,4:ncol(dfBlock)]))
    matBlock <- scale(matBlock)
    print(paste0("The dim of matBlock is ", dim(matBlock)))

    matf <- matBlock %*% as.matrix(dfR_block$r)
    matf_raw <- apply(matf, 1, sum)
    dfFGr_mat[,counter] <- matf_raw

    # Save number of SNPs
    nsnp_in_block <- nrow(dfR_block)
    dfSNPs[counter,1] <- b
    dfSNPs[counter,2] <- nrow(dfR_block)
    counter <- counter + 1

  }
}

# Save Raw FGR
dfFGr <- as.data.frame(dfFGr_mat)
fwrite(dfFGr, out_file_FGr, quote = F, row.names = F, sep = "\t")

# Save SNP file
fwrite(dfSNPs, out_file_SNP, quote = F, row.names = F, sep = "\t")



