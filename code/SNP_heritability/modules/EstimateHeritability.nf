#!/usr/bin/env nextflow

process EstimateHeritability {
    container 'quay.io/cawarmerdam/ldsc:v0.3'
    tag "${name}"
    
    // Define inputs
    input:
        tuple val(name), val(prevalence), file(sumstats), file(w_ld_chr)

    // Define outputs
    output:
        tuple path("${name}_heritability_liability.log"), path("${name}_heritability_observed.log")

    // Script to execute
    script:
    """
    /ldsc/ldsc.py \
    --h2 ${sumstats} \
    --ref-ld-chr ${w_ld_chr}/ \
    --w-ld-chr ${w_ld_chr}/ \
    --chisq-max 10000 \
    --samp-prev ${prevalence} \
    --pop-prev ${prevalence} \
    --out ${name}_heritability_liability

    /ldsc/ldsc.py \
    --h2 ${sumstats} \
    --ref-ld-chr ${w_ld_chr}/ \
    --w-ld-chr ${w_ld_chr}/ \
    --chisq-max 10000 \
    --out ${name}_heritability_observed
    """
}
