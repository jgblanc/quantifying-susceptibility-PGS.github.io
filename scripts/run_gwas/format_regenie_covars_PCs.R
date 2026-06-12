## Format phenotype file for REGENIE

args=commandArgs(TRUE)

if(length(args)<2){stop("Rscript calc_R2_EO.R")}


library(data.table)
library(tidyverse)

info_file = args[1]
common_file = args[2]
rare_file = args[3]
out_file = args[4]


# Read in phenotypes
df <- fread(info_file)

# Add IID column
df <- df[, c(1, 1, 2:4)]

# Name columns
colnames(df) <- c("FID", "IID", "Age", "Sex", "Batch")
df <- df %>% mutate(Batch = replace(Batch, Batch < 0, 0)) %>%
  mutate(Batch = replace(Batch, Batch > 0, 1))

# Read in PCs
dfCommon <- fread(common_file)
colnames(dfCommon)[1] <- "FID"
dfRare <- fread(rare_file)
colnames(dfRare)[1] <- "FID"
names <- colnames(dfRare)[3:42]
new_names <- paste0("R", names)
colnames(dfRare)[3:42] <- new_names
dfPCs <- cbind(dfCommon, dfRare[,3:42])

# Join Dataframes
dfOut <- inner_join(df, dfPCs)


# Write output
fwrite(dfOut, out_file,row.names = FALSE, sep = "\t", quote = FALSE,na = "NA")
