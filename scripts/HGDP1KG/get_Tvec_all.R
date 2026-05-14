# This script computes the allele frequency difference between two contintental "populations"

args=commandArgs(TRUE)

if(length(args)<11){stop("Rscript compute_test_vec.R")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
}))


## Load in IDs
dfIDs <- fread(args[11])
dfIDs <- dfIDs %>% select(FID, IID, pop)


# Go through all pairwise comparisons

# EAS-NFE
df <- dfIDs %>% filter(pop %in% c("eas", "nfe")) %>% mutate(tvec = case_when(pop == "eas" ~ 1, pop == "nfe" ~ -1))
df$tvec<- scale(df$tvec)
print(paste0("The mean of Tvec is ", mean(df$tvec)))
print(paste0("The variance of Tvec is ", var(df$tvec)))
fwrite(df, args[1] ,row.names=F,quote=F,sep="\t", col.names = T)


# EAS-SAS
df <- dfIDs %>% filter(pop %in% c("eas", "sas"))%>% mutate(tvec = case_when(pop == "eas" ~ 1, pop == "sas" ~ -1))
df$tvec<- scale(df$tvec)
print(paste0("The mean of Tvec is ", mean(df$tvec)))
print(paste0("The variance of Tvec is ", var(df$tvec)))
fwrite(df, args[2] ,row.names=F,quote=F,sep="\t", col.names = T)

# EAS-AFR
df <- dfIDs %>% filter(pop %in% c("eas", "afr"))%>% mutate(tvec = case_when(pop == "eas" ~ 1, pop == "afr" ~ -1))
df$tvec<- scale(df$tvec)
print(paste0("The mean of Tvec is ", mean(df$tvec)))
print(paste0("The variance of Tvec is ", var(df$tvec)))
fwrite(df, args[3] ,row.names=F,quote=F,sep="\t", col.names = T)

# EAS-AMR
df <- dfIDs %>% filter(pop %in% c("eas", "amr"))%>% mutate(tvec = case_when(pop == "eas" ~ 1, pop == "amr" ~ -1))
df$tvec<- scale(df$tvec)
print(paste0("The mean of Tvec is ", mean(df$tvec)))
print(paste0("The variance of Tvec is ", var(df$tvec)))
fwrite(df, args[4] ,row.names=F,quote=F,sep="\t", col.names = T)


# NFE-SAS
df <- dfIDs %>% filter(pop %in% c("nfe", "sas")) %>% mutate(tvec = case_when(pop == "nfe" ~ 1, pop == "sas" ~ -1))
df$tvec<- scale(df$tvec)
print(paste0("The mean of Tvec is ", mean(df$tvec)))
print(paste0("The variance of Tvec is ", var(df$tvec)))
fwrite(df, args[5] ,row.names=F,quote=F,sep="\t", col.names = T)

# NFE-AFR
df <- dfIDs %>% filter(pop %in% c("nfe", "afr")) %>% mutate(tvec = case_when(pop == "nfe" ~ 1, pop == "afr" ~ -1))
df$tvec<- scale(df$tvec)
print(paste0("The mean of Tvec is ", mean(df$tvec)))
print(paste0("The variance of Tvec is ", var(df$tvec)))
fwrite(df, args[6] ,row.names=F,quote=F,sep="\t", col.names = T)

# NFE-AMR
df <- dfIDs %>% filter(pop %in% c("nfe", "amr")) %>% mutate(tvec = case_when(pop == "nfe" ~ 1, pop == "amr" ~ -1))
df$tvec<- scale(df$tvec)
print(paste0("The mean of Tvec is ", mean(df$tvec)))
print(paste0("The variance of Tvec is ", var(df$tvec)))
fwrite(df, args[7] ,row.names=F,quote=F,sep="\t", col.names = T)

# SAS-AFR
df <- dfIDs %>% filter(pop %in% c("sas", "afr")) %>% mutate(tvec = case_when(pop == "sas" ~ 1, pop == "afr" ~ -1))
df$tvec<- scale(df$tvec)
print(paste0("The mean of Tvec is ", mean(df$tvec)))
print(paste0("The variance of Tvec is ", var(df$tvec)))
fwrite(df, args[8] ,row.names=F,quote=F,sep="\t", col.names = T)

# SAS-AMR
df <- dfIDs %>% filter(pop %in% c("sas", "amr")) %>% mutate(tvec = case_when(pop == "sas" ~ 1, pop == "amr" ~ -1))
df$tvec<- scale(df$tvec)
print(paste0("The mean of Tvec is ", mean(df$tvec)))
print(paste0("The variance of Tvec is ", var(df$tvec)))
fwrite(df, args[9] ,row.names=F,quote=F,sep="\t", col.names = T)

# AFR-AMR
df <- dfIDs %>%  filter(pop %in% c("afr", "amr")) %>% mutate(tvec = case_when(pop == "afr" ~ 1, pop == "amr" ~ -1))
df$tvec<- scale(df$tvec)
print(paste0("The mean of Tvec is ", mean(df$tvec)))
print(paste0("The variance of Tvec is ", var(df$tvec)))
fwrite(df, args[10] ,row.names=F,quote=F,sep="\t", col.names = T)



















