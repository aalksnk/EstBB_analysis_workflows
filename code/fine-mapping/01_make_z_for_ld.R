#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Usage: Rscript make_z_from_reference.R <snp_info.gz> <output_folder> <chr> <ld_start> <ld_end>")
}

snp_info_file <- args[1]
output_folder <- args[2]
chr   <- as.integer(args[3])
start <- as.numeric(args[4])  # ld_start
end   <- as.numeric(args[5])  # ld_end

message("Processing reference LD region chr", chr, ":", start, "-", end)


dir.create(output_folder, showWarnings = FALSE, recursive = TRUE)

out_file <- file.path(
  output_folder,
  paste0("reference_LD_", chr, "_", start, "_", end, ".z")
)

if (file.exists(out_file)) {
  message("File already exists for chr", chr, ":", start, "-", end,
          " -> ", out_file, " ; skipping.")
  quit(save = "no", status = 0)
}

message("Processing reference LD region chr", chr, ":", start, "-", end)


# Load reference
snp_info <- data.table::fread(snp_info_file)

# Basic required columns to define a region
needed <- c("ID", "CHROM", "GENPOS")
missing <- setdiff(needed, names(snp_info))
if (length(missing) > 0) {
  stop("snp_info missing required columns: ", paste(missing, collapse = ", "))
}

# Robust chromosome match (handles X/23)
chrom_chr <- as.character(snp_info$CHROM)
is_chr_match <- if (!is.na(chr) && chr == 23) {
  chrom_chr %in% c("23", "X")
} else {
  chrom_chr == as.character(chr)
}

cropped <- snp_info %>%
  dplyr::filter(
    is_chr_match,
    GENPOS >= start,
    GENPOS <= end
  )

if (nrow(cropped) == 0) {
  message("No variants found in reference for region ", chr, ":", start, "-", end)
  quit(save = "no", status = 0)
}

# --- If this is chr 23, change CHROM to "X" BEFORE creating ID_long ---
if (!is.na(chr) && chr == 23) {
  cropped <- cropped %>% mutate(CHROM = "X")
}

# ---- Keep your previous manipulation logic, but only if columns exist ----
# (These likely won't exist in snp_info; this is safe-guarded.)
if ("BETA" %in% names(cropped)) {
  cropped <- cropped %>% mutate(BETA = -BETA)
}
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

# --- Create ID_long and drop duplicates if possible ---
if (all(c("CHROM", "GENPOS", "ALLELE0", "ALLELE1") %in% names(cropped))) {
  cropped <- cropped %>%
    mutate(
      CHROM   = as.character(CHROM),
      ID_long = paste(CHROM, GENPOS, ALLELE0, ALLELE1, sep = "_")
    ) %>%
    distinct(ID_long, .keep_all = TRUE)

  message("Created ID_long and removed duplicates: remaining rows = ", nrow(cropped))
}

# Optional A1FREQ filter if present
if ("A1FREQ" %in% names(cropped)) {
  cropped <- cropped %>%
    mutate(A1FREQ = as.numeric(A1FREQ)) %>%
    filter(!is.na(A1FREQ) & A1FREQ >= 0.001)
}

if (nrow(cropped) == 0) {
  message("All variants filtered out after A1FREQ filter in chr", chr, ":", start, "-", end)
  quit(save = "no", status = 0)
}

# Ensure beta/se columns exist (reference likely doesn't have them)
if (!("BETA" %in% names(cropped))) cropped$BETA <- NA_real_
if (!("SE"   %in% names(cropped))) cropped$SE   <- NA_real_
if (!("A1FREQ" %in% names(cropped))) cropped$A1FREQ <- NA_real_

# Build z in the format you want
# Note: This is now a "reference-style z" unless your snp_info already has BETA/SE.
z <- cropped %>%
  transmute(
    rsid       = ID,
    chromosome = as.integer(CHROM),
    position   = as.integer(GENPOS),
    allele1    = if ("ALLELE0" %in% names(cropped)) as.character(ALLELE0) else NA_character_,
    allele2    = if ("ALLELE1" %in% names(cropped)) as.character(ALLELE1) else NA_character_,
    maf        = as.numeric(A1FREQ),
    beta       = as.numeric(BETA),
    se         = as.numeric(SE)
  )

# Optional maf filter only when maf is available
if (any(!is.na(z$maf))) {
  z <- z %>% filter(is.na(maf) | maf >= 0.001)
}


# Make chr X => 23 in z
z$chromosome[is.na(z$chromosome)] <- 23

write.table(
  z,
  file = out_file,
  sep = " ",
  quote = FALSE,
  row.names = FALSE
)

message("Saved: ", out_file, " (", nrow(z), " rows)")
