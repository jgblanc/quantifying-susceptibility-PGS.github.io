## Project Tvec per chromosome
## This script projects the test vector from the test to the gwas panel using genotypes for a single chromosome

args=commandArgs(TRUE)

if(length(args)<6){stop("Rscript project_Tvec_chr.R <test panel prefix> <gwas panel prefix> <test vec file> <outfile>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
}))

overlap_snps = args[1]
plink_prefix = args[2]


####################
## Functions #######
####################

# Compute t(X) %*% T
compute_b <- function(tvec_file, outfile) {

  # Compute t(X)T
  outpath <- strsplit(outfile, ".rvec")[[1]][1]
  outfile_XT <- paste0(outpath, "_xt_temp")
  cmd_XT <- paste("sh code/HGDP1KG/compute_XT.sh", plink_prefix, tvec_file, outfile_XT, overlap_snps, sep = " ")
  system(cmd_XT)

  # Adjust Betas to account for variance in x

  # Read in betas and genotype counts
  beta_plink <- fread(paste0(outpath, "_xt_temp.tvec.glm.linear"))
  count_plink <- fread(paste0(outpath, "_xt_temp.gcount"))

  # Calculate length of mean centered genotypes from counts
  nOBS <- (count_plink$HOM_REF_CT + count_plink$HET_REF_ALT_CTS + count_plink$TWO_ALT_GENO_CTS)
  counts <- (count_plink$HOM_REF_CT * 0) + (count_plink$HET_REF_ALT_CTS * 1) + (count_plink$TWO_ALT_GENO_CTS * 2)
  mean_gc <- counts / nOBS
  length_mc_genos <- (count_plink$HOM_REF_CT * (-1 * mean_gc)^2) + (count_plink$HET_REF_ALT_CTS * (1 - mean_gc)^2) +  (count_plink$TWO_ALT_GENO_CTS * (2 - mean_gc)^2)

  # Fix betas
  betas_plink_norm <- beta_plink$BETA * sqrt(length_mc_genos/nOBS)

  #  Re-write .linear file with correct betas
  beta_plink$BETA <- betas_plink_norm
  beta_reformat <- beta_plink %>% dplyr::select("#CHROM","ID", "REF", "ALT",  "BETA")
  beta_reformat[is.na(beta_reformat)] <- 0

  # Remove tmp files
  cmd_rm <- paste0("rm ", outfile_XT, "*")
  system(cmd_rm)

  return(beta_reformat)
}


#####################
##     Main       ###
#####################

# Compute EAS-NFE
r = compute_b(args[3], args[13])
colnames(r)[5] <- "r"
fwrite(r, args[13], row.names = F, col.names = T, quote = F, sep = "\t")

# Compute EAS-SAS
r = compute_b(args[4], args[14])
colnames(r)[5] <- "r"
fwrite(r, args[14], row.names = F, col.names = T, quote = F, sep = "\t")

# Compute EAS-AFR
r = compute_b(args[5], args[15])
colnames(r)[5] <- "r"
fwrite(r, args[15], row.names = F, col.names = T, quote = F, sep = "\t")

# Compute EAS-AMR
r = compute_b(args[6], args[16])
colnames(r)[5] <- "r"
fwrite(r, args[16], row.names = F, col.names = T, quote = F, sep = "\t")

# Compute NFE-SAS
r = compute_b(args[7], args[17])
colnames(r)[5] <- "r"
fwrite(r, args[17], row.names = F, col.names = T, quote = F, sep = "\t")

# Compute NFE-AFR
r = compute_b(args[8], args[18])
colnames(r)[5] <- "r"
fwrite(r, args[18], row.names = F, col.names = T, quote = F, sep = "\t")

# Compute NFE-AMR
r = compute_b(args[9], args[19])
colnames(r)[5] <- "r"
fwrite(r, args[19], row.names = F, col.names = T, quote = F, sep = "\t")

# Compute SAS-AFR
r = compute_b(args[10], args[20])
colnames(r)[5] <- "r"
fwrite(r, args[20], row.names = F, col.names = T, quote = F, sep = "\t")

# Compute SAS-AMR
r = compute_b(args[11], args[21])
colnames(r)[5] <- "r"
fwrite(r, args[21], row.names = F, col.names = T, quote = F, sep = "\t")

# Compute AFR-AMR
r = compute_b(args[12], args[22])
colnames(r)[5] <- "r"
fwrite(r, args[22], row.names = F, col.names = T, quote = F, sep = "\t")






































