## Overlapping SNPs
## This script takes two sets of .afreq files and subsets each to variants over 1% and then gets the intersection of the two lists

args=commandArgs(TRUE)

if(length(args)<3){stop("Rscript overlapping_snps.R <ukbb.freq> <test panel.freq> <outfile>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
}))

ukbb_file = args[1]
tp_file = args[2]
outfile = args[3]

# Load freqency files
ukbb <- fread(ukbb_file)
tp <- fread(tp_file)
tp <- tp %>% group_by(ID) %>% filter(n() == 1) %>% ungroup()
print(nrow(tp))
print(nrow(ukbb))

# Subset to SNPs with missingness rate < 5%
#tp <- subset(tp, tp$OBS_CT > (0.95 * max(tp$OBS_CT)))   
#ukbb <- subset(ukbb, ukbb$OBS_CT > (0.95 *  max(ukbb$OBS_CT)))


# Get overlapping SNPs with same alt/ref
matched <- inner_join(ukbb, tp, by = c("#CHROM", "ID", "ALT", "REF"))
matched <- matched %>% select(`#CHROM`,ID,REF,ALT )
print(nrow(matched))

# Get overlapping SNPs with opposite alt/ref
tmp1 <- inner_join(ukbb, tp, by = c("#CHROM", "ID"))
flipped <- subset(tmp1, tmp1$REF.x == tmp1$ALT.y & tmp1$ALT.x == tmp1$REF.y)
flipped <- flipped %>% select("#CHROM", ID, REF.x, REF.y)
colnames(flipped) <- c("#CHROM", "ID", "REF", "ALT")
print(nrow(flipped))

# Get complete list of SNPs and list the reference allele in UKBB
out <- rbind(matched, flipped) %>% select(ID, REF) 
fwrite(out,outfile, row.names = F, col.names = T, quote = F, sep = "\t")

