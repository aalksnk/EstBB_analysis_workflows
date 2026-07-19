# process_results.R
args <- commandArgs(trailingOnly = TRUE)
library(data.table)
library(arrow)
library(stringr) 

final_results <- fread(args[1], stringsAsFactors = FALSE)
sample_count_file <- args[2]
final_results$Phenotype <- str_replace(final_results$Phenotype, "_processed.txt", "")

# Read and calculate prevalence
prev <- fread(sample_count_file)
prev$prevalence <- round(prev$cases / (prev$controls + prev$cases), 6)

# Define the ZtoP function
ZtoP <- function(Z, largeZ = FALSE, log10P = TRUE) {
  if (!is.numeric(Z)) {
    message("Some of the Z-scores are not numbers! Please check why!")
    Z <- as.numeric(Z)
  }
  if (largeZ) {
    P <- log(2) + pnorm(abs(Z), lower.tail = FALSE, log.p = TRUE)
    if (log10P) {
      P <- -(P * log10(exp(1)))
    }
  } else {
    P <- 2 * pnorm(abs(Z), lower.tail = FALSE)
    if (min(P) == 0) {
      P[P == 0] <- .Machine$double.xmin
      message("Some Z-score indicates very significant effect and P-value is truncated on 2.22e-308. Consider using largeZ = TRUE and logarithmed P-values.")
    }
  }
  return(P)
}

# Calculate P-values using h2 and se_h2 columns
final_results$P_value <- ZtoP(final_results$h2 / final_results$se_h2)
final_results <- as.data.frame(final_results)
prev <- as.data.frame(prev)

# Merge cases, controls, and prevalence data
final_data <- merge(final_results, prev, by.x = "Phenotype", by.y = "code", all.x = TRUE)

# Write to a file
fwrite(final_data, args[3], sep = "\t", row.names = FALSE, quote = FALSE)
