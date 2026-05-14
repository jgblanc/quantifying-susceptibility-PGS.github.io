## Get Lat/Long Test Vector
## This script uses HGDP metadata and a psam to format test vectors for latitide and longitude

args=commandArgs(TRUE)

if(length(args)<3){stop("Rscript get_LatLong_Tvec.R <psam> <outfile Lat> <outfile Long>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
}))

md_file = args[1]
outfile_lat = args[2]
outfile_long = args[3]

# Read in metadata
df <- fread(md_file)

# Set up outfiles
dfLat <- df %>% select("FID","IID", "lat")
colnames(dfLat)[3] <- "tvec"
dfLat$tvec <- scale(dfLat$tvec)
print(paste0("The mean of Tvec is ", mean(dfLat$tvec)))
print(paste0("The variance of Tvec is ", var(dfLat$tvec)))

dfLong <- df %>% select("FID","IID", "long")
colnames(dfLong)[3] <- "tvec"
dfLong$tvec <- scale(dfLong$tvec)
print(paste0("The mean of Tvec is ", mean(dfLong$tvec)))
print(paste0("The variance of Tvec is ", var(dfLong$tvec)))

# Save output
fwrite(dfLat,outfile_lat, row.names = F, col.names = T, quote = F, sep = "\t")
fwrite(dfLong,outfile_long, row.names = F, col.names = T, quote = F, sep = "\t")
