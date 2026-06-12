## Format phenotype file for REGENIE

args=commandArgs(TRUE)

if(length(args)<2){stop("Rscript calc_R2_EO.R")}


library(data.table)
library(tidyverse)

in_file = args[1]
out_file = args[2]


# Read in phenotypes
df <- fread(in_file)

# Add IID column
df <- df[, c(1, 1, 2:18)]

# Name columns
colnames(df) <- c("FID", "IID", "Standing_Height", "Alkaline_Phosphate", "Aspartate_Aminotransferase", "Basophill_Percentage",
                  "Body_Weight", "Cholesterol", "Eosinophill_Percentage", "Glucose", "HbA1c", "HDL", "LDL", "MCHC",
                  "MCV", "Platelet_Count", "SBP", "Total_Protein", "Triglycerides")


# Write output
write.table(df, file = out_file,
            sep = "\t", na = "NA", quote = FALSE, row.names = FALSE, fileEncoding = "UTF-8")

