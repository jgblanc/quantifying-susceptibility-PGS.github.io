## Compute FGr using one PC SNPs

args=commandArgs(TRUE)

if(length(args)<3){stop("Rscript calc_FGr.R <prefix to plink files> <r prefix> <var prefix> <out prefix> <snp> <ids> <outfile>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(tidyverse)
}))

plink_prefix = args[1]
out_prefix = args[2]
id_file = args[3]

# Read in IDs
dfIDs <- fread(id_file)

# Set data collectors
M <- nrow(dfIDs)


# Do chromsome 1 separately

# Set up parameters
chr_str <- paste0("1")
plink_out <- paste0(plink_prefix, "_1")

# Run plink
plink_cmd <- paste0("plink2 --pfile ", plink_prefix, " --keep ", id_file, " --chr ", chr_str ," --threads 8 --pca 40 approx --memory 38000 ", "--out ", plink_out)
system(plink_cmd)


# Loop through Chromosome
for (i in 2:22) {

  print(paste0("Running plinke for chromosomes 1 through ", i))

  # Set up parameters
  chr_str <- paste0("1-", i)
  plink_out <- paste0(plink_prefix, "_", "i")

  # Run plink
  plink_cmd <- paste0("plink2 --pfile ", plink_prefix, " --keep ", id_file, " --chr ", chr_str ," --threads 8 --pca 40 approx --memory 38000 ", "--out ", plink_out)
  system(plink_cmd)

}




