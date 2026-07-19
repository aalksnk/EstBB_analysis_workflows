#!/usr/bin/env Rscript
library(arrow)
library(dplyr)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 8) {
  stop("Usage: Rscript process_regions_by_pheno.R <input_folder> <snp_info.gz> <output_folder> <regions_file> <chr> <start> <end> <pheno>")
}


input_folder  <- args[1]
snp_info_file <- args[2]
output_folder <- args[3]
regions_file <- args[4]
chr   <- as.numeric(args[5])
start <- as.numeric(args[6])
end   <- as.numeric(args[7])
pheno <- args[8]

message("Loading regions file: ", regions_file)
regions <- data.table::fread(regions_file)

required_cols <- c("phenotype", "chr", "start", "end")
missing_cols <- setdiff(required_cols, names(regions))
if (length(missing_cols) > 0) {
  stop("Regions file is missing required columns: ", paste(missing_cols, collapse = ", "))
}

# Coerce essential columns
regions[, phenotype := as.character(phenotype)]
regions[, chr       := as.integer(chr)]
regions[, start     := as.numeric(start)]
regions[, end       := as.numeric(end)]

# Drop obviously bad rows
regions <- regions[!is.na(phenotype) & !is.na(chr) & !is.na(start) & !is.na(end)]
if (nrow(regions) == 0) {
  stop("No valid rows found in regions file after cleaning.")
}

# Load reference once
message("Loading SNP reference: ", snp_info_file)
snp_info <- data.table::fread(snp_info_file)

# Ensure output dir exists
dir.create(output_folder, showWarnings = FALSE, recursive = TRUE)

# Helper to filter chromosome robustly (handles numeric/character and X/23)
chrom_match_filter <- function(df, chr_num) {
  # CHROM may be integer, numeric, or character
  chrom_chr <- as.character(df$CHROM)

  if (!is.na(chr_num) && chr_num == 23) {
    # accept "23" or "X"
    chrom_chr %in% c("23", "X")
  } else {
    chrom_chr == as.character(chr_num)
  }
}

# Process phenotype by phenotype to avoid re-reading parquet repeatedly
phenos <- unique(regions$phenotype)

for (pheno in phenos) {
  message("\n=== Phenotype: ", pheno, " ===")

  parquet_file <- file.path(input_folder, paste0("results_concat_", pheno, ".parquet.snappy"))
  if (!file.exists(parquet_file)) {
    warning("Parquet file not found for phenotype ", pheno, ": ", parquet_file, " | Skipping.")
    next
  }

  message("Reading parquet: ", parquet_file)
  dat <- arrow::read_parquet(parquet_file, as_data_frame = TRUE)

  # Merge with reference
  merged <- dplyr::left_join(dat, snp_info, by = "ID")

  # Get all regions for this phenotype
  regs <- regions[phenotype == pheno]

  # Loop regions
  for (i in seq_len(nrow(regs))) {
    chr_i   <- regs$chr[i]
    start_i <- regs$start[i]
    end_i   <- regs$end[i]

    message("Processing region chr", chr_i, ":", start_i, "-", end_i)

    # Robust chromosome filter
    is_chr_match <- chrom_match_filter(merged, chr_i)

    cropped <- merged %>%
      dplyr::filter(
        is_chr_match,
        GENPOS >= start_i,
        GENPOS <= end_i
      )

    if (nrow(cropped) == 0) {
      message("No variants found in region chr", chr_i, ":", start_i, "-", end_i, " for ", pheno)
      next
    }

    # --- If this is chr 23, change CHROM to "X" BEFORE creating ID_long ---
    if (!is.na(chr_i) && chr_i == 23) {
      cropped <- cropped %>% mutate(CHROM = "X")
    }

    # Transformations (kept from your original script)
    if ("BETA" %in% names(cropped)) {
      cropped <- cropped %>% mutate(BETA = -BETA)
    }

    # Safe swap of ALLELE0 and ALLELE1
    if (all(c("ALLELE0", "ALLELE1") %in% names(cropped))) {
      tmpA0 <- cropped$ALLELE0
      cropped$ALLELE0 <- cropped$ALLELE1
      cropped$ALLELE1 <- tmpA0
    }

    if ("A1FREQ" %in% names(cropped)) {
      cropped <- cropped %>% mutate(A1FREQ = 1 - A1FREQ)
    }

    if ("LOG10P" %in% names(cropped)) {
      cropped <- cropped %>% mutate(P = 10^(-LOG10P))
    }

    # --- Create ID_long and drop duplicates ---
    if (all(c("CHROM", "GENPOS", "ALLELE0", "ALLELE1") %in% names(cropped))) {
      cropped <- cropped %>%
        mutate(
          CHROM   = as.character(CHROM),
          ID_long = paste(CHROM, GENPOS, ALLELE0, ALLELE1, sep = "_")
        ) %>%
        distinct(ID_long, .keep_all = TRUE)

      message("Created ID_long and removed duplicates: remaining rows = ", nrow(cropped))
    } else {
      warning("Some columns missing for ID_long creation in region chr", chr_i, ":", start_i, "-", end_i)
    }

    # Optional MAF filter if A1FREQ exists
    if ("A1FREQ" %in% names(cropped)) {
      cropped <- cropped %>%
        mutate(A1FREQ = as.numeric(A1FREQ)) %>%
        filter(!is.na(A1FREQ) & A1FREQ >= 0.001)
    }

    if (nrow(cropped) == 0) {
      message("All variants filtered out after A1FREQ filter in chr", chr_i, ":", start_i, "-", end_i, " for ", pheno)
      next
    }

    # Output only parquet
    out_file <- file.path(
      output_folder,
      paste0("processed_", chr_i, "_", start_i, "_", end_i, "_", pheno, ".parquet.snappy")
    )

    arrow::write_parquet(cropped, out_file, compression = "snappy")
    message("Saved: ", out_file, " (", nrow(cropped), " rows)")
  }
}

message("\nDone.")
