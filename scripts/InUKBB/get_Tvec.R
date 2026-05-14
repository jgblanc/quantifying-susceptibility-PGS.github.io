# This script computes the allele frequency difference between two contintental "populations"

args=commandArgs(TRUE)

if(length(args)<11){stop("Rscript compute_test_vec.R")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
}))


## Load in IDs
dfIDs <- fread(args[11])
dfIDs <- dfIDs %>% select("#FID", IID, POP)
colnames(dfIDs) <- c("FID", "IID", "pop")

# Go through all pairwise comparisons

# ENG-NI
df <- dfIDs %>% filter(pop %in% c("England", "NorthernIreland")) %>% mutate(tvec = case_when(pop == "England" ~ 1, pop == "NorthernIreland" ~ -1))
df$tvec<- scale(df$tvec)
print(paste0("The mean of Tvec is ", mean(df$tvec)))
print(paste0("The variance of Tvec is ", var(df$tvec)))
fwrite(df, args[1] ,row.names=F,quote=F,sep="\t", col.names = T)


# ENG-ROI
df <- dfIDs %>% filter(pop %in% c("England", "RepublicOfIrelend"))%>% mutate(tvec = case_when(pop == "England" ~ 1, pop == "RepublicOfIrelend" ~ -1))
df$tvec<- scale(df$tvec)
print(paste0("The mean of Tvec is ", mean(df$tvec)))
print(paste0("The variance of Tvec is ", var(df$tvec)))
fwrite(df, args[2] ,row.names=F,quote=F,sep="\t", col.names = T)

# ENG-SCT
df <- dfIDs %>% filter(pop %in% c("England", "Scotland"))%>% mutate(tvec = case_when(pop == "England" ~ 1, pop == "Scotland" ~ -1))
df$tvec<- scale(df$tvec)
print(paste0("The mean of Tvec is ", mean(df$tvec)))
print(paste0("The variance of Tvec is ", var(df$tvec)))
fwrite(df, args[3] ,row.names=F,quote=F,sep="\t", col.names = T)

# EAS-AMR
df <- dfIDs %>% filter(pop %in% c("England", "Wales"))%>% mutate(tvec = case_when(pop == "England" ~ 1, pop == "Wales" ~ -1))
df$tvec<- scale(df$tvec)
print(paste0("The mean of Tvec is ", mean(df$tvec)))
print(paste0("The variance of Tvec is ", var(df$tvec)))
fwrite(df, args[4] ,row.names=F,quote=F,sep="\t", col.names = T)


# NFE-SAS
df <- dfIDs %>% filter(pop %in% c("NorthernIreland", "RepublicOfIrelend")) %>% mutate(tvec = case_when(pop == "NorthernIreland" ~ 1, pop == "RepublicOfIrelend" ~ -1))
df$tvec<- scale(df$tvec)
print(paste0("The mean of Tvec is ", mean(df$tvec)))
print(paste0("The variance of Tvec is ", var(df$tvec)))
fwrite(df, args[5] ,row.names=F,quote=F,sep="\t", col.names = T)

# NFE-AFR
df <- dfIDs %>% filter(pop %in% c("NorthernIreland", "Scotland")) %>% mutate(tvec = case_when(pop == "NorthernIreland" ~ 1, pop == "Scotland" ~ -1))
df$tvec<- scale(df$tvec)
print(paste0("The mean of Tvec is ", mean(df$tvec)))
print(paste0("The variance of Tvec is ", var(df$tvec)))
fwrite(df, args[6] ,row.names=F,quote=F,sep="\t", col.names = T)

# NFE-AMR
df <- dfIDs %>% filter(pop %in% c("NorthernIreland", "Wales")) %>% mutate(tvec = case_when(pop == "NorthernIreland" ~ 1, pop == "Wales" ~ -1))
df$tvec<- scale(df$tvec)
print(paste0("The mean of Tvec is ", mean(df$tvec)))
print(paste0("The variance of Tvec is ", var(df$tvec)))
fwrite(df, args[7] ,row.names=F,quote=F,sep="\t", col.names = T)

# SAS-AFR
df <- dfIDs %>% filter(pop %in% c("RepublicOfIrelend", "Scotland")) %>% mutate(tvec = case_when(pop == "RepublicOfIrelend" ~ 1, pop == "Scotland" ~ -1))
df$tvec<- scale(df$tvec)
print(paste0("The mean of Tvec is ", mean(df$tvec)))
print(paste0("The variance of Tvec is ", var(df$tvec)))
fwrite(df, args[8] ,row.names=F,quote=F,sep="\t", col.names = T)

# SAS-AMR
df <- dfIDs %>% filter(pop %in% c("RepublicOfIrelend", "Wales")) %>% mutate(tvec = case_when(pop == "RepublicOfIrelend" ~ 1, pop == "Wales" ~ -1))
df$tvec<- scale(df$tvec)
print(paste0("The mean of Tvec is ", mean(df$tvec)))
print(paste0("The variance of Tvec is ", var(df$tvec)))
fwrite(df, args[9] ,row.names=F,quote=F,sep="\t", col.names = T)

# AFR-AMR
df <- dfIDs %>%  filter(pop %in% c("Scotland", "Wales")) %>% mutate(tvec = case_when(pop == "Scotland" ~ 1, pop == "Wales" ~ -1))
df$tvec<- scale(df$tvec)
print(paste0("The mean of Tvec is ", mean(df$tvec)))
print(paste0("The variance of Tvec is ", var(df$tvec)))
fwrite(df, args[10] ,row.names=F,quote=F,sep="\t", col.names = T)



















