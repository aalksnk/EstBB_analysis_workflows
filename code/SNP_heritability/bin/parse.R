#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
library(data.table)
library(arrow)
library(stringr)

input_file <- args[1]
snp_ref <- args[2]
hapmap_list <- args[3]
sample_count_file <- args[4]

pheno <- str_replace(args[1], "\\.parquet.snappy", "")
pheno <- str_replace(pheno, "results_concat_", "")

output_file <- paste0(pheno, "_processed.txt")

# Read the input parquet file
sumstats <- read_parquet(input_file)
message("Input read!")
# Convert to a data.table and calculate P from LOG10P
sumstats <- data.table(SNP = sumstats$ID, BETA = sumstats$BETA, SE = sumstats$SE, P = sumstats$LOG10P)
setkey(sumstats, "SNP")
message("Indexed!")
# Filter based on HapMap variants
hapmap <- fread(args[3])
sumstats <- sumstats[SNP %in% hapmap$SNP]
sumstats$P <- 10^(-sumstats$P)
message("Sumstats filtered!")
# Add alleles
message("Reading in SNP reference...")
snpref <- read_parquet(args[2])
setkey(snpref, "ID")
message("Indexed!")
snpref <- snpref[ID %in% hapmap$SNP]
message("Ref filtered!")
sumstats <- merge(sumstats, snpref, by.x = "SNP", by.y = "ID")
message("Merged!")
# Extract sample prevalence
prev <- fread(args[4])
prev <- prev[code %in% pheno]
N <- prev$cases + prev$controls

sumstats <- data.table(
    SNP = sumstats$SNP,
    A1 = sumstats$ALLELE1,
    A2 = sumstats$ALLELE0,
    N = N,
    BETA = sumstats$BETA,
    SE = sumstats$SE,
    P = sumstats$P,
    EAF = sumstats$A1FREQ,
    INFO = sumstats$INFO
)
message("File prepared!")
prev <- prev$cases / (prev$controls + prev$cases)
write(prev, stdout())
# Write the processed data to a text file
fwrite(sumstats, output_file, sep = "\t")
message("Output written!")
