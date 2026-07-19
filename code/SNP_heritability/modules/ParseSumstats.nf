#!/usr/bin/env nextflow

process ParseSumstats {
    // conda 'environment.yml'
    container 'quay.io/urmovosa/ld_overlap_image:v0.2'
    tag "${file.baseName}"
    label 'R'  // Assuming 'R' label is defined in your nextflow.config for R-specific resources

    input:
        tuple path(file), path(snpref), path(sample_count_file), path(w_ld_chr)

    output:
        tuple env(phenoname), env(prevalence), path("*_processed.txt")

    script:
    """
    prevalence=\$(Rscript --vanilla ${baseDir}/bin/parse.R \
    ${file} \
    ${snpref} \
    ${w_ld_chr}/w_hm3.snplist \
    ${sample_count_file})

    # Parse phenotype name
    phenoname=\$(ls *_processed* | sed 's/^results_concat_\\(.*\\)\\(\\.parquet\\)*\\.snappy\$/\\1/')
    """
}
