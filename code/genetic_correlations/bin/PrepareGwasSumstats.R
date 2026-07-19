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

if (input_file == "NA" || !file.exists(input_file)) {
    stop("Error: input_file does not exist or is set to 'NA'")
}

if (snp_ref == "NA" || !file.exists(snp_ref)) {
    stop("Error: snp_ref does not exist or is set to 'NA'")
}

if (hapmap_list == "NA" || !file.exists(hapmap_list)) {
    stop("Error: hapmap_list does not exist or is set to 'NA'")
}

if (sample_count_file == "NA" || !file.exists(hapmap_list)) {
    stop("Error: case_control_file does not exist or is set to 'NA'")
}

# Read the input parquet file
sumstats <- read_parquet(input_file)
message("Input read!")
# Convert to a data.table and calculate P from LOG10P
sumstats <- data.table(SNP = sumstats$ID, BETA = sumstats$BETA, SE = sumstats$SE, P = sumstats$LOG10P)
setkey(sumstats, "SNP")
message("Indexed!")
# Filter based on HapMap variants
hapmap <- fread(args[3])
#sumstats <- sumstats[SNP %in% hapmap$SNP]
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


gwas <- data.table(
    SNP = sumstats$SNP,
    A1 = sumstats$REF,
    A2 = sumstats$ALT,
    N = N,
    BETA = sumstats$BETA *(-1),
    SE = sumstats$SE,
    P = sumstats$P,
    EAF = sumstats$ALT_AF,
    INFO = sumstats$INFO
)

message("File prepared!")
prev <- prev$cases / (prev$controls + prev$cases)
write(prev, stdout())
# Write the processed data to a text file
fwrite(gwas, output_file, sep = "\t")
message("Output written!")

