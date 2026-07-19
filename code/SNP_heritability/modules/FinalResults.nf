#!/bin/bash nextflow 

process FinaliseResults {
    container 'quay.io/urmovosa/ld_overlap_image:v0.2'
    publishDir "${params.outputDir}", mode: 'copy', overwrite: true

    input:
        path results
        path case_control_file

    output:
        path("results_h2_liability.tsv")
        path("results_h2_observed.tsv")

    script:
    """
    Rscript --vanilla ${baseDir}/bin/process_results.R final_results_table_liability.tsv ${case_control_file} results_h2_liability.tsv
    Rscript --vanilla ${baseDir}/bin/process_results.R final_results_table_observed.tsv ${case_control_file} results_h2_observed.tsv
    """
}
