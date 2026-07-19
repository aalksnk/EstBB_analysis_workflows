#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
options(error = function() { traceback(2); q(status = 1) })

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(readr)
  library(arrow)
  library(susieR)
  library(dotgen)
  library(Rfast)
  library(stringr)
})

# ---------- helpers ----------
clean_str <- function(x) {
  x <- gsub("\r", "", x, fixed = TRUE)     # CR
  x <- gsub("\u00A0", "", x, fixed = TRUE) # NBSP
  trimws(x)
}


normalize_chr_to_int <- function(chr) {
  chr <- clean_str(as.character(chr))
  chr <- sub("^chr", "", chr, ignore.case = TRUE)
  if (chr %in% c("X", "x", "23")) return(23L)
  if (chr %in% c("Y", "y", "24")) return(24L)
  if (chr %in% c("MT", "Mt", "mt", "M", "m", "25")) return(25L)
  out <- suppressWarnings(as.integer(chr))
  if (is.na(out)) stop("Could not parse chromosome from '", chr, "'.")
  out
}

chr_label_set <- function(chr_int) {
  if (chr_int == 23L) return(c("23", "X", "x", "chr23", "chrX", "chrx"))
  if (chr_int == 24L) return(c("24", "Y", "y", "chr24", "chrY", "chry"))
  if (chr_int == 25L) return(c("25", "MT", "Mt", "mt", "M", "m", "chr25", "chrMT", "chrM"))
  c(as.character(chr_int), paste0("chr", chr_int))
}

make_region_stubs <- function(chr_int, start, end) {
  unique(paste0(chr_label_set(chr_int), ":", format(start, scientific = FALSE), "-", format(end, scientific = FALSE)))
}

normalize_variant_chr_prefix <- function(x) {
  x <- gsub("\r|\u00A0", "", x)
  x <- trimws(x)
  x <- sub("^chr", "", x, ignore.case = TRUE)
  x <- sub("^(X|x)_",  "23_", x)
  x <- sub("^(Y|y)_",  "24_", x)
  x <- sub("^(MT|Mt|mt|M|m)_", "25_", x)
  x
}

# Force an object to be a true 2D numeric matrix. This avoids errors such as
# "attempt to set 'colnames' on an object with less than two dimensions" when
# R accidentally simplifies a one-row/one-column object to a vector.

set_colnames_2d <- function(x, value, object_name = "object") {
  if (is.null(dim(x)) || length(dim(x)) != 2) {
    stop(object_name, " is not two-dimensional before assigning column names. class=",
         paste(class(x), collapse = ","), "; length=", length(x))
  }
  if (length(value) != ncol(x)) {
    stop(object_name, " column-name length mismatch: ncol=", ncol(x),
         "; names=", length(value))
  }
  colnames(x) <- value
  x
}

set_rownames_2d <- function(x, value, object_name = "object") {
  if (is.null(dim(x)) || length(dim(x)) != 2) {
    stop(object_name, " is not two-dimensional before assigning row names. class=",
         paste(class(x), collapse = ","), "; length=", length(x))
  }
  if (length(value) != nrow(x)) {
    stop(object_name, " row-name length mismatch: nrow=", nrow(x),
         "; names=", length(value))
  }
  rownames(x) <- value
  x
}

force_numeric_matrix <- function(x, object_name = "matrix") {
  if (is.data.frame(x) || data.table::is.data.table(x)) {
    x <- as.matrix(x)
  }

  if (is.null(dim(x))) {
    len <- length(x)
    root <- sqrt(len)
    if (len > 0 && root == floor(root)) {
      x <- matrix(as.numeric(x), nrow = root, ncol = root, byrow = TRUE)
    } else {
      stop(object_name, " has less than two dimensions and cannot be converted safely to a square matrix. Length = ", len)
    }
  }

  x <- as.matrix(x)
  storage.mode(x) <- "numeric"
  x
}

compute_yty <- function(beta, se, p, R, n, k) {
  beta_s <- beta * sqrt(2 * p * (1 - p))
  se_s   <- se   * sqrt(2 * p * (1 - p))
  XjtXj  <- (n - 1) * diag(R)
  median(beta_s**2 * XjtXj + se_s**2 * XjtXj * (n - k))
}

summarize.susie.cs <- function(object, orig_vars, R, ..., low_purity_threshold = 0.5) {
  if (is.null(object$sets)) stop("Cannot summarize SuSiE object because credible set information is not available")
  variables <- data.frame(cbind(1:length(object$pip), object$pip, -1, NA, NA, NA))
  variables <- set_colnames_2d(variables, c("variable","variable_prob","cs","cs_specific_prob","low_purity","lead_r2"), "summary variables")
  if (object$null_index > 0) variables <- variables[-object$null_index, ]
  added_vars <- c()
  if (!is.null(object$sets$cs)) {
    cs <- data.frame(matrix(NA, length(object$sets$cs), 5))
    cs <- set_colnames_2d(cs, c("cs","cs_log10bf","cs_avg_r2","cs_min_r2","variable"), "summary credible sets")
    for (i in seq_along(object$sets$cs)) {
      if (any(object$sets$cs[[i]] %in% added_vars)) next else added_vars <- c(added_vars, object$sets$cs[[i]])
      in_cs_idx <- which(variables$variable %in% object$sets$cs[[i]])
      variables$cs[in_cs_idx] <- object$sets$cs_index[[i]]
      variables$cs_specific_prob[in_cs_idx] <- object$alpha[object$sets$cs_index[[i]], object$sets$cs[[i]]]
      variables$low_purity[in_cs_idx] <- object$sets$purity$min.abs.corr[i] < low_purity_threshold
      lead_pip_idx <- in_cs_idx[which.max(variables$variable_prob[in_cs_idx])]
      variables$lead_r2 <- R[lead_pip_idx, ]^2
      cs$cs[i]         <- object$sets$cs_index[[i]]
      cs$cs_log10bf[i] <- log10(exp(object$lbf[cs$cs[i]]))
      cs$cs_avg_r2[i]  <- object$sets$purity$mean.abs.corr[i]^2
      cs$cs_min_r2[i]  <- object$sets$purity$min.abs.corr[i]^2
      cs$variable[i]   <- paste(object$sets$cs[[i]], collapse = ",")
    }
    variables <- variables[order(variables$variable_prob, decreasing = TRUE), , drop = FALSE]
  } else cs <- NULL
  list(vars = variables, cs = na.omit(cs))
}

susie_ss_wrapper <- function(df, R, n, L,
                             estimate_residual_variance = TRUE,
                             var_y = 1, prior_weights = NULL,
                             min_abs_corr = 0.5, low_purity_threshold = 0.5) {
  fitted_bhat <- susie_rss(
    bhat = df$beta, shat = df$se, R = R, n = n, var_y = var_y, L = L,
    prior_weights = prior_weights, scaled_prior_variance = 0.1,
    estimate_residual_variance = estimate_residual_variance, estimate_prior_variance = TRUE,
    standardize = TRUE, check_input = FALSE, min_abs_corr = min_abs_corr
  )
  cs_summary <- summarize.susie.cs(fitted_bhat, df, R, low_purity_threshold = low_purity_threshold)
  variables  <- cs_summary$vars %>% as.data.frame() %>%
    dplyr::rename(prob = variable_prob) %>%
    arrange(variable) %>%
    mutate(mean = susie_get_posterior_mean(fitted_bhat),
           sd   = susie_get_posterior_sd(fitted_bhat))
  cs <- cs_summary$cs

  sets_95 <- fitted_bhat$sets
  fitted_bhat$sets <- susieR::susie_get_cs(fitted_bhat, coverage = 0.99, Xcorr = R, min_abs_corr = min_abs_corr)
  cs_summary_99 <- summarize.susie.cs(fitted_bhat, df, R, low_purity_threshold = low_purity_threshold)
  fitted_bhat$sets_99 <- fitted_bhat$sets
  fitted_bhat$sets <- sets_95

  variables_99 <- cs_summary_99$vars %>% as.data.frame() %>%
    dplyr::rename(prob = variable_prob) %>%
    arrange(variable) %>%
    mutate(mean = susie_get_posterior_mean(fitted_bhat),
           sd   = susie_get_posterior_sd(fitted_bhat))
  variables_99 <- set_colnames_2d(variables_99, paste0(colnames(variables_99), "_99"), "variables_99")

  list(
    susie_obj = fitted_bhat,
    variables = variables,
    variables_99 = variables_99,
    cs = cs,
    cs_99 = cs_summary_99$cs
  )
}

report_fit <- function(fit, R, prefix) {
  warns <- attr(fit, "warnings")
  if (length(warns)) {
    message(prefix, ": SuSiE warnings (", length(warns), "):")
    message(paste0("  - ", unique(warns)), collapse = "\n")
  }

  message(prefix, ": alpha dim: ", paste(dim(fit$alpha), collapse="x"))
  message(prefix, ": max PIP: ", signif(max(fit$pip, na.rm=TRUE), 4),
          " | #PIP>0.01: ", sum(fit$pip > 0.01, na.rm=TRUE))

  # Per-effect concentration
  mx <- apply(fit$alpha, 1, max, na.rm = TRUE)
  message(prefix, ": per-effect max(alpha): ",
          paste(signif(mx, 3), collapse = ", "))

  # Null weight if present
  if (!is.null(fit$null_index) && fit$null_index > 0) {
    null_mass <- fit$alpha[, fit$null_index]
    message(prefix, ": per-effect null alpha: ",
            paste(signif(null_mass, 3), collapse = ", "))
  }

  # CS counts
  n_cs <- if (is.null(fit$sets) || is.null(fit$sets$cs)) 0 else length(fit$sets$cs)
  message(prefix, ": CS count: ", n_cs)

  if (!is.null(fit$sets) && !is.null(fit$sets$purity)) {
    message(prefix, ": CS purity min.abs.corr: ",
            paste(signif(fit$sets$purity$min.abs.corr, 3), collapse=", "))
  }
}



ensure_cols <- function(df, cols) {
  for (cc in cols) {
    if (!cc %in% names(df)) df[[cc]] <- NA
  }
  df
}

# ---------- args ----------
# New interface (11 args):
# 1 chr
# 2 ld_start
# 3 ld_end
# 4 gwas_start
# 5 gwas_end
# 6 LD_folder
# 7 results_folder
# 8 phenos_raw
# 9 n_samples
# 10 gwas_folder
# 11 vars_folder
if (length(args) < 11) {
  stop("Usage: Rscript run_susie_ld_gwas_map.R <chr> <ld_start> <ld_end> <gwas_start> <gwas_end> <LD_folder> <results_folder> <phenos> <n_samples> <gwas_folder> <vars_folder>")
}

locus_chr_raw <- clean_str(args[1])
ld_start      <- as.numeric(args[2])
ld_end        <- as.numeric(args[3])
gwas_start    <- as.numeric(args[4])
gwas_end      <- as.numeric(args[5])

LD_region_folder     <- clean_str(args[6])
coloc_results_folder <- clean_str(args[7])
phenos_raw           <- clean_str(args[8])         # "PH1;PH2" or "PH1"
n_samples_raw        <- clean_str(args[9])
gwas_folder          <- clean_str(args[10])
vars_folder          <- clean_str(args[11])

locus_chr_int <- normalize_chr_to_int(locus_chr_raw)
locus_chr_labels <- chr_label_set(locus_chr_int)
message("==== CHR DEBUG: INPUT ARGS ====")
message("args[1] raw chr: '", args[1], "'")
message("locus_chr_raw after clean_str: '", locus_chr_raw, "'")
message("locus_chr_int: ", locus_chr_int)
message("Accepted chromosome labels for file search: ", paste(locus_chr_labels, collapse = ", "))
message("ld_start: ", ld_start, " ld_end: ", ld_end)
message("gwas_start: ", gwas_start, " gwas_end: ", gwas_end)

# ---------- tunables ----------
max_causal_SNPs <- 10
n_covariates    <- 13
min_cs_corr     <- 0.5
low_purity_threshold <- 0.5
GRCh <- 38
gwas_suffix <- ".parquet.snappy"
LD_suffix   <- ".ld.gz"
VARS_suffix <- ".vars"

# ---------- parse ----------
n_samples <- suppressWarnings(readr::parse_number(n_samples_raw))
if (is.na(n_samples)) stop("Could not parse n_samples from '", n_samples_raw, "'.")

phenos <- unlist(strsplit(phenos_raw, ";", fixed = TRUE))
phenos <- clean_str(phenos)
phenos <- phenos[nchar(phenos) > 0]

message("LD region  chr", locus_chr_int, ":", ld_start, "-", ld_end)
message("GWAS region chr", locus_chr_int, ":", gwas_start, "-", gwas_end)
message("LD folder: ", LD_region_folder)
message("Vars folder: ", vars_folder)
message("Results: ", coloc_results_folder, " | GWAS folder: ", gwas_folder)
message("Phenotypes: ", paste(phenos, collapse = ", "), " | N = ", n_samples)

ld_qc_folder <- file.path(coloc_results_folder, "ld_qc")
dir.create(ld_qc_folder, showWarnings = FALSE, recursive = TRUE)

# ---------- LD files ----------
region_stubs_ld <- make_region_stubs(locus_chr_int, ld_start, ld_end)
region_stub_ld <- region_stubs_ld[1]

message("==== CHR DEBUG: LD REGION ====")
message("Candidate LD region stubs: ", paste(region_stubs_ld, collapse = ", "))

ld_candidates <- file.path(LD_region_folder, paste0(region_stubs_ld, LD_suffix))
ld_candidates <- ld_candidates[file.exists(ld_candidates)]

vars_candidates <- unlist(lapply(region_stubs_ld, function(stub) {
  list.files(
    vars_folder,
    pattern = paste0("^", gsub("([.|()\\^{}+$*?]|\\[|\\])", "\\\\\\1", stub), ".*\\.vars$"),
    full.names = TRUE
  )
}), use.names = FALSE)
vars_candidates <- unique(vars_candidates)

if (length(ld_candidates) == 0 || length(vars_candidates) == 0) {
  message("Could not find LD inputs for any candidate chromosome label.")
  message("Expected LD candidates:")
  print(ld_candidates)
  message("First 30 files in LD folder:")
  print(utils::head(list.files(LD_region_folder, full.names = FALSE), 30))
  message("First 30 files in vars folder:")
  print(utils::head(list.files(vars_folder, full.names = FALSE), 30))
  message("Skipping LD region.")
  quit(save = "no", status = 0)
}

LD_file <- ld_candidates[1]
vars_file <- vars_candidates[1]

message("Using LD_file:   ", LD_file)
message("Using vars_file: ", vars_file)
message("LD exists?   ", file.exists(LD_file))
message("vars exists? ", file.exists(vars_file))

if (!file.exists(LD_file) || !file.exists(vars_file)) {
  message("Missing LD inputs for LD region. Skipping.")
  quit(save = "no", status = 0)
}

# Read vars + LD matrix
# Important: do NOT pass duplicated vars directly as row.names to read.table().
# Some .vars files can contain duplicate IDs, or duplicates can be created after
# normalizing X_... to 23_.... We read the LD matrix without names first, then
# remove duplicated variants from both vars and the matching LD rows/columns.
vars_raw <- readr::read_tsv(
  vars_file,
  col_names = FALSE,
  show_col_types = FALSE
)[[1]]

vars <- normalize_variant_chr_prefix(vars_raw)

message("==== CHR DEBUG: VARS ====")
message("Number of vars: ", length(vars))
message("First 10 raw vars:")
print(utils::head(vars_raw, 10))
message("First 10 cleaned vars:")
print(utils::head(vars, 10))
message("Cleaned variant chromosome/prefix table:")
print(table(sub("_.*$", "", vars), useNA = "ifany"))

LD_matrix_raw <- read.table(
  LD_file,
  fill = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

LD_matrix <- force_numeric_matrix(LD_matrix_raw, object_name = "LD_matrix")

if (is.null(dim(LD_matrix)) || length(dim(LD_matrix)) != 2) {
  stop("LD matrix was not read as a 2D matrix. Check LD file format: ", LD_file)
}

message("Raw LD matrix dim: ", paste(dim(LD_matrix), collapse = " x "))
message("Raw vars length: ", length(vars))

if (nrow(LD_matrix) != length(vars) || ncol(LD_matrix) != length(vars)) {
  stop(
    "LD matrix dimensions do not match vars length before duplicate removal. LD dim = ",
    paste(dim(LD_matrix), collapse = " x "),
    "; vars length = ", length(vars),
    ". Check that the .ld.gz and .vars files come from the same LD region. ",
    "If LD dim is N x 1, your LD file may not be a whitespace-delimited full matrix."
  )
}

if (anyDuplicated(vars)) {
  dup_ids <- unique(vars[duplicated(vars)])
  message("Duplicated variant IDs in vars after cleaning: ", length(dup_ids))
  message("First duplicated IDs: ", paste(utils::head(dup_ids, 20), collapse = ", "))

  dup_report <- data.frame(
    row_index = seq_along(vars),
    raw_id = vars_raw,
    cleaned_id = vars,
    stringsAsFactors = FALSE
  )
  dup_report <- dup_report[dup_report$cleaned_id %in% dup_ids, , drop = FALSE]

  dup_report_file <- file.path(
    ld_qc_folder,
    paste0("duplicated_vars_", locus_chr_int, "_", ld_start, "_", ld_end, ".tsv")
  )
  data.table::fwrite(dup_report, dup_report_file, sep = "\t")
  message("Duplicate vars report written: ", dup_report_file)

  # Keep the first occurrence of each cleaned variant ID and drop later duplicates
  # from BOTH the vars vector and the corresponding LD matrix rows/columns.
  keep_idx <- !duplicated(vars)
  vars <- vars[keep_idx]
  LD_matrix <- LD_matrix[keep_idx, keep_idx, drop = FALSE]
  LD_matrix <- force_numeric_matrix(LD_matrix, object_name = "LD_matrix_after_dedup")

  message("After removing duplicate variants: LD matrix dim = ",
          paste(dim(LD_matrix), collapse = " x "),
          "; vars length = ", length(vars))
}

if (is.null(dim(LD_matrix)) || length(dim(LD_matrix)) != 2) {
  stop("LD matrix lost its 2D dimensions before assigning row/column names.")
}
if (nrow(LD_matrix) != length(vars) || ncol(LD_matrix) != length(vars)) {
  stop("Cannot assign LD row/column names: LD dim = ",
       paste(dim(LD_matrix), collapse = " x "),
       "; vars length = ", length(vars))
}

LD_matrix <- set_rownames_2d(LD_matrix, vars, "LD_matrix")
LD_matrix <- set_colnames_2d(LD_matrix, vars, "LD_matrix")
message("LD matrix dim: ", paste(dim(LD_matrix), collapse = " x "))

# ---------- per phenotype ----------
for (pheno in phenos) {

  pheno_clean <- gsub("\u00A0", "", pheno, fixed = TRUE)
  pheno_clean <- gsub("\\s+", "", pheno_clean)

  # GWAS file based on SMALL GWAS region
  gwas_file <- file.path(
    gwas_folder,
    paste0(
      "processed_", locus_chr_int, "_",
      format(gwas_start, scientific = FALSE), "_",
      format(gwas_end,   scientific = FALSE), "_",
      pheno_clean, gwas_suffix
    )
  )

  message("\n[", pheno_clean, "] GWAS file: ", gwas_file)

  # Output prefix includes both LD and GWAS coords
  output_prefix <- file.path(
    coloc_results_folder,
    paste0(
      pheno_clean,
      "_LD_", locus_chr_int, "_", format(ld_start, scientific = FALSE), "_", format(ld_end, scientific = FALSE),
      "_GWAS_", format(gwas_start, scientific = FALSE), "_", format(gwas_end, scientific = FALSE)
    )
  )

  if (!file.exists(gwas_file)) {
    cand <- list.files(gwas_folder, pattern = "\\.parquet\\.snappy$", full.names = TRUE)
    tags <- paste0(
      "processed_", chr_label_set(locus_chr_int), "_",
      format(gwas_start, scientific = FALSE), "_",
      format(gwas_end, scientific = FALSE), "_",
      pheno_clean
    )
    cand <- cand[Reduce(`|`, lapply(tags, function(tag) grepl(tag, basename(cand), fixed = TRUE)))]
    if (length(cand) >= 1) {
      gwas_file <- cand[1]
      message("Using matched GWAS file: ", gwas_file)
    } else {
      message("GWAS not found; writing NULL outputs.")
      file.create(paste0(output_prefix, "_NULL_coloc5.tsv"))
      file.create(paste0(output_prefix, "_NULL_coloc3.tsv"))
      file.create(paste0(output_prefix, "_NULL_clpp.tsv"))
      next
    }
  }

  gwas <- arrow::read_parquet(gwas_file, as_data_frame = TRUE)

  message("==== CHR DEBUG: GWAS RAW [", pheno_clean, "] ====")
  message("GWAS columns:")
  print(names(gwas))
  if ("CHROM" %in% names(gwas)) {
    message("Raw CHROM class: ", paste(class(gwas$CHROM), collapse = ","))
    message("Raw CHROM values table:")
    print(table(gwas$CHROM, useNA = "ifany"))
  }

  # --- Robust CHROM normalization ---
  if (!"CHROM" %in% names(gwas)) stop("GWAS file missing CHROM column.")
gwas$CHROM <- vapply(gwas$CHROM, normalize_chr_to_int, integer(1))

  message("==== CHR DEBUG: GWAS NORMALIZED [", pheno_clean, "] ====")
  message("Normalized CHROM class: ", paste(class(gwas$CHROM), collapse = ","))
  message("Normalized CHROM values table:")
  print(table(gwas$CHROM, useNA = "ifany"))
  message("Looking for locus_chr_int = ", locus_chr_int)
  message("Rows matching chromosome before position filter: ",
          sum(gwas$CHROM == locus_chr_int, na.rm = TRUE))

  # position column
  pos_col <- if ("GENPOS" %in% names(gwas)) "GENPOS" else if ("BP" %in% names(gwas)) "BP" else if ("POS" %in% names(gwas)) "POS" else NA
  if (is.na(pos_col)) stop("No position column (GENPOS/BP/POS) in GWAS.")

  rng <- range(gwas[[pos_col]], na.rm = TRUE)
message("Position col=", pos_col, " range: ", paste(rng, collapse=" - "))
  message("Position range by chromosome:")
  print(
    gwas %>%
      dplyr::group_by(CHROM) %>%
      dplyr::summarise(
        n = dplyr::n(),
        min_pos = suppressWarnings(min(.data[[pos_col]], na.rm = TRUE)),
        max_pos = suppressWarnings(max(.data[[pos_col]], na.rm = TRUE)),
        .groups = "drop"
      )
  )

  if (!"A1FREQ" %in% names(gwas)) stop("GWAS missing A1FREQ.")
  if (!"ID" %in% names(gwas)) stop("GWAS missing ID column for LD alignment.")
  if (!all(c("ALLELE0","ALLELE1") %in% names(gwas))) stop("GWAS missing ALLELE0/ALLELE1.")
  if (!all(c("BETA","SE") %in% names(gwas))) stop("GWAS missing BETA/SE.")

  # Filter by SMALL GWAS region + MAF
  gwas <- gwas %>%
    mutate(A1FREQ = as.numeric(A1FREQ)) %>%
    filter(
      CHROM == locus_chr_int,
      dplyr::between(.data[[pos_col]], gwas_start, gwas_end),
      dplyr::between(A1FREQ, 0.001, 0.999)
    ) %>%
    arrange(CHROM, .data[[pos_col]])


  message("Rows after GWAS region+MAF: ", nrow(gwas))

  if (nrow(gwas) == 0) {
    file.create(paste0(output_prefix, "_NULL_coloc5.tsv"))
    file.create(paste0(output_prefix, "_NULL_coloc3.tsv"))
    file.create(paste0(output_prefix, "_NULL_clpp.tsv"))
    next
  }

  # ---------- Align GWAS to LD vars/order ----------
  common_ids <- intersect(gwas$ID, vars)

  if (length(common_ids) < 2) {
    message("Too few overlapping variants between GWAS and LD for ", pheno_clean, ". Writing NULL outputs.")
    file.create(paste0(output_prefix, "_NULL_coloc5.tsv"))
    file.create(paste0(output_prefix, "_NULL_coloc3.tsv"))
    file.create(paste0(output_prefix, "_NULL_clpp.tsv"))
    next
  }

  gwas <- gwas %>% filter(ID %in% common_ids)
  ord  <- match(gwas$ID, vars)
  gwas <- gwas[order(ord), ]

  ids_ordered <- gwas$ID

  # LD submatrix (by aligned IDs)
  R <- as.matrix(LD_matrix[ids_ordered, ids_ordered, drop = FALSE])

  if (nrow(R) == 0 || ncol(R) == 0) {
    file.create(paste0(output_prefix, "_NULL_coloc5.tsv"))
    file.create(paste0(output_prefix, "_NULL_coloc3.tsv"))
    file.create(paste0(output_prefix, "_NULL_clpp.tsv"))
    next
  }

  R[upper.tri(R)] <- t(R)[upper.tri(t(R))]
  R[is.na(R)] <- 0

  # ---------- LD QC ----------
  offdiag_ones <- sum(R[row(R) != col(R)] == 1)

  eigvals <- tryCatch(
    eigen(R, symmetric = TRUE, only.values = TRUE)$values,
    error = function(e) NA_real_
  )
  min_eig <- if (all(is.na(eigvals))) NA_real_ else min(eigvals, na.rm = TRUE)

  qc_df <- data.frame(
    chr  = locus_chr_int,
    start = ld_start,
    end   = ld_end,
    n_offdiag_ones = offdiag_ones,
    min_eigenvalue = min_eig
  )

  qc_file <- file.path(
    ld_qc_folder,
    paste0("ld_qc_", locus_chr_int, "_", ld_start, "_", ld_end, "_", pheno_clean, ".tsv")
  )
  readr::write_tsv(qc_df, qc_file)
  message("LD QC written: ", qc_file)

  # ---------- analysis df ----------
  gwas <- gwas %>%
    mutate(
      A1FREQ = as.numeric(.data[["A1FREQ"]]),
      maf = pmin(A1FREQ, 1 - A1FREQ)
    ) %>%
    rename(
      beta = BETA, se = SE,
      rsid = ID,
      position = !!pos_col,
      chromosome = CHROM,
      allele1 = "ALLELE1",
      allele2 = "ALLELE0"
    ) %>%
    mutate(CHR = chromosome)


  n <- n_samples
  L <- min(max_causal_SNPs, nrow(gwas))
  if (L < 1) {
    file.create(paste0(output_prefix, "_NULL_coloc5.tsv"))
    file.create(paste0(output_prefix, "_NULL_coloc3.tsv"))
    file.create(paste0(output_prefix, "_NULL_clpp.tsv"))
    next
  }

  message("Step 1. Variants in analysis (after LD overlap): ", nrow(gwas))
  message("Step 2. Compute var_y")
  yty <- compute_yty(beta = gwas$beta, se = gwas$se, p = gwas$maf, R = R, n = n, k = n_covariates)
  var_y <- yty / (n - 1)

  message("Step 3. Run SuSiE RSS")
  res <- tryCatch(
    susie_ss_wrapper(
      df = gwas, R = R, n = n, L = L,
      estimate_residual_variance = TRUE, var_y = var_y,
      prior_weights = NULL, min_abs_corr = min_cs_corr,
      low_purity_threshold = low_purity_threshold
    ),
    error = function(cond) { message("SuSiE error: ", conditionMessage(cond)); NULL }
  )

  if (is.null(res)) {
    file.create(paste0(output_prefix, "_NULL_coloc5.tsv"))
    file.create(paste0(output_prefix, "_NULL_coloc3.tsv"))
    file.create(paste0(output_prefix, "_NULL_clpp.tsv"))
    next
  }

  report_fit(res$susie_obj, R, prefix = pheno_clean)

  # ---------- PIP/CS 95 & 99 ----------
  vars95 <- res$variables
  vars95$rsid <- gwas$rsid[vars95$variable]
  vars95 <- vars95[, c("rsid","variable","prob","cs","cs_specific_prob","low_purity","mean","sd"), drop = FALSE]
  data.table::fwrite(vars95, paste0(output_prefix, ".pip_95.tsv"), sep = "\t")

  if (!is.null(res$cs) && nrow(res$cs) > 0) {
    data.table::fwrite(res$cs, paste0(output_prefix, ".cs_95.tsv"), sep = "\t")
  }

  vars99 <- res$variables_99
  vars99$rsid <- gwas$rsid[vars99$variable_99]
  vars99 <- vars99[, c("rsid","variable_99","prob_99","cs_99","cs_specific_prob_99","low_purity_99","mean_99","sd_99"), drop = FALSE]
  data.table::fwrite(vars99, paste0(output_prefix, ".pip_99.tsv"), sep = "\t")

  if (!is.null(res$cs_99) && nrow(res$cs_99) > 0) {
    data.table::fwrite(res$cs_99, paste0(output_prefix, ".cs_99.tsv"), sep = "\t")
  }

  # ---------- coloc-style tables ----------
  susie_obj <- res$susie_obj
  fit <- res$susie_obj
  fit$sets <- susieR::susie_get_cs(fit, coverage = 0.95, Xcorr = R, min_abs_corr = min_cs_corr)
  message(pheno_clean, ": CS@0.95 count: ", ifelse(is.null(fit$sets$cs), 0, length(fit$sets$cs)))

  fit$sets99 <- susieR::susie_get_cs(fit, coverage = 0.99, Xcorr = R, min_abs_corr = min_cs_corr)
  message(pheno_clean, ": CS@0.99 count: ", ifelse(is.null(fit$sets99$cs), 0, length(fit$sets99$cs)))

  lbf_variables <- susie_obj$lbf_variable
  if (is.null(lbf_variables)) {
    stop("SuSiE object has no lbf_variable matrix; cannot create coloc5 output.")
  }
  lbf_variables <- as.matrix(lbf_variables)
  if (is.null(dim(lbf_variables)) || length(dim(lbf_variables)) != 2) {
    stop("susie_obj$lbf_variable is not two-dimensional; cannot create coloc5 output.")
  }

  transposed_lbf <- as.data.frame(t(lbf_variables))
  if (ncol(transposed_lbf) < 1) {
    stop("No lbf_variable columns found after transposing SuSiE lbf_variable.")
  }
  transposed_lbf <- set_colnames_2d(transposed_lbf, paste0("lbf_variable", seq_len(ncol(transposed_lbf))), "transposed_lbf")

  transposed_lbf$rsid         <- gwas$rsid
  transposed_lbf$old_position <- as.numeric(gwas$position)
  transposed_lbf$A1           <- toupper(trimws(gwas$allele1))
  transposed_lbf$A2           <- toupper(trimws(gwas$allele2))

  old_rsid <- gwas$rsid
  old_position <- gwas$position

  gwas2 <- gwas %>%
    rename(
      SNP = rsid, BP = position,
      A1 = allele1, A2 = allele2
    )

  gwas2$rsid <- old_rsid
  gwas2$old_position <- old_position

  if (!"CHR" %in% names(gwas2)) {
    if ("chromosome" %in% names(gwas2)) {
      gwas2$CHR <- gwas2$chromosome
    } else {
      stop("Neither CHR nor chromosome column exists in gwas2.")
    }
  }

  message("gwas2 chromosome values before final label conversion: ",
          paste(utils::head(unique(gwas2$CHR), 20), collapse = ", "))

  # Keep chromosome numeric internally. Convert 23/24/25 only for final output labels.
  gwas2$CHR <- as.character(gwas2$CHR)
  gwas2$CHR[gwas2$CHR == "23"] <- "X"
  gwas2$CHR[gwas2$CHR == "24"] <- "Y"
  gwas2$CHR[gwas2$CHR == "25"] <- "MT"

  message("gwas2 chromosome values after final label conversion: ",
          paste(utils::head(unique(gwas2$CHR), 20), collapse = ", "))

  if (GRCh == 38) {
    df_hg38 <- gwas2
  } else {
    # only call if package available in your container
    df_hg38 <- MungeSumstats::liftover(gwas2, convert_ref_genome = "hg38", ref_genome = "hg19")
  }

  message("N before liftover: ", nrow(gwas2), " | after: ", nrow(df_hg38))

  df_hg38$molecular_trait_id <- pheno_clean
  df_hg38$region <- paste0("chr", unique(df_hg38$CHR), ":", min(df_hg38$BP), "-", max(df_hg38$BP))
  df_hg38$SNP <- paste0("chr", df_hg38$CHR, "_", df_hg38$BP, "_", df_hg38$A2, "_", df_hg38$A1)

  df_hg38$old_position <- as.numeric(df_hg38$old_position)
  df_hg38$A1 <- toupper(trimws(df_hg38$A1))
  df_hg38$A2 <- toupper(trimws(df_hg38$A2))

  final_lbf <- dplyr::inner_join(
    df_hg38,
    transposed_lbf,
    by = c("rsid","old_position","A1","A2")
  )
  final_lbf <- final_lbf[order(final_lbf$BP), , drop = FALSE]

  lbf_cols <- grep("^lbf_variable\\d+$", names(final_lbf), value = TRUE)
  if (nrow(final_lbf) == 0L || length(lbf_cols) == 0L) {
    file.create(paste0(output_prefix, "_NULL_coloc5.tsv"))
    message("coloc5: no rows after join (or no LBF cols); wrote NULL.")
  } else {
    data_coloc5 <- final_lbf[, c("molecular_trait_id","region","SNP","CHR","BP", lbf_cols), drop = FALSE]
    names(data_coloc5)[1:5] <- c("molecular_trait_id","region","variant","chromosome","position")
    write.table(data_coloc5, paste0(output_prefix, "_coloc5.tsv"),
                row.names = FALSE, sep = "\t", quote = FALSE)
    message("Coloc5 written (", nrow(data_coloc5), " rows, ", length(lbf_cols), " LBF cols).")
  }

  # coloc3 (be robust to missing INFO/LOG10P)
  df_hg38 <- ensure_cols(df_hg38, c("maf","beta","se","LOG10P","INFO"))

  data_coloc3 <- df_hg38[, c("molecular_trait_id","region","SNP","A1","A2","CHR","BP","maf","beta","se","LOG10P","INFO"), drop = FALSE]
  names(data_coloc3) <- c("molecular_trait_id","region","variant","ref","alt",
                           "chromosome","position","maf","beta","se","log10p","info")
  write.table(data_coloc3, paste0(output_prefix, "_coloc3.tsv"),
              row.names = FALSE, sep = "\t", quote = FALSE)
  message("Coloc3 written.")

  # ---------- CLPP ----------
  variables <- cbind(rsid = gwas$rsid, res$variables)
  variables$old_position <- gwas$position
  variables$A1 <- gwas$allele1
  variables$A2 <- gwas$allele2
  variables$LOG10P <- if ("LOG10P" %in% names(gwas)) gwas$LOG10P else NA_real_

  alpha_mat <- susie_obj$alpha
  if (is.null(alpha_mat)) {
    stop("SuSiE object has no alpha matrix; cannot create CLPP output.")
  }
  alpha_mat <- as.matrix(alpha_mat)
  if (is.null(dim(alpha_mat)) || length(dim(alpha_mat)) != 2) {
    stop("susie_obj$alpha is not two-dimensional; cannot create CLPP output.")
  }

  pip <- as.data.frame(t(alpha_mat))
  if (ncol(pip) < 1) {
    stop("No alpha columns found after transposing SuSiE alpha.")
  }
  pip <- set_colnames_2d(pip, paste0("alpha", seq_len(ncol(pip))), "pip alpha table")
  pip$rsid <- rownames(pip)
  pip$pip <- susie_obj$pip
  pip$old_position <- gwas$position
  pip$A1 <- gwas$allele1
  pip$A2 <- gwas$allele2

  variables <- merge(variables, pip, by = c("rsid","old_position","A1","A2"))
  variables_clpp <- variables %>% filter(cs > 0)

  cs <- res$cs
  if (is.null(cs)) {
    write.table(variables_clpp, paste0(output_prefix, "_NULL_clpp.tsv"),
                row.names = FALSE, sep = "\t")
  } else {
    variables_clpp <- merge(variables_clpp, cs, by = "cs")
    variables_clpp_hg38 <- dplyr::inner_join(
      df_hg38,
      variables_clpp,
      by = c("rsid","old_position","A1","A2")
    )

    message("CS before join: ", nrow(variables_clpp), " | after join: ", nrow(variables_clpp_hg38))

    if (nrow(variables_clpp_hg38) == 0) {
      file.create(paste0(output_prefix, "_NULL_clpp.tsv"))
    } else {
      if ("LOG10P" %in% names(variables_clpp_hg38)) {
        pvals <- 10^-(variables_clpp_hg38$LOG10P)
      } else if ("P" %in% names(variables_clpp_hg38)) {
        pvals <- variables_clpp_hg38$P
      } else {
        pvals <- NA_real_
      }

      if (all(!is.finite(pvals))) {
        variables_clpp_hg38$z <- NA_real_
      } else {
        variables_clpp_hg38$z <- zsc(pvals, variables_clpp_hg38$beta)
      }

      variables_clpp_hg38$cs_index <- paste0("L", variables_clpp_hg38$cs)
      variables_clpp_hg38$cs_id <- paste0(
        variables_clpp_hg38$molecular_trait_id, "_",
        variables_clpp_hg38$region, "_",
        variables_clpp_hg38$cs_index
      )

      alpha_cols <- grep("^alpha\\d+$", names(variables_clpp_hg38), value = TRUE)
      data_clpp <- variables_clpp_hg38[, c(
        "molecular_trait_id","region","SNP","CHR","BP","A1","A2",
        "cs_id","cs_index", alpha_cols, "pip","z"
      ), drop = FALSE]

      names(data_clpp) <- c(
        "molecular_trait_id","region","variant","chromosome","position",
        "ref","alt","cs_id","cs_index", alpha_cols, "pip","z"
      )

      write.table(data_clpp, paste0(output_prefix, "_clpp.tsv"),
                  row.names = FALSE, sep = "\t", quote = FALSE)
      message("CLPP written.")
    }
  }

  rm(R)
  gc()
}

message(date())
message("Done")
