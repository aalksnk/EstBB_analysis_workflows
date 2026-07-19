#!/bin/bash nextflow

process ProcessGwas {
    tag "${name}"
    container 'quay.io/cawarmerdam/ldsc:v0.3'

    input:
        tuple val(name), val(prevalence), path(sumstats), path(w_ld_chr_ch), path(ldsc)  // Receives phenotype name and sumstats path from previous process

    output:
        tuple val(name), val(prevalence), path("${name.replaceAll("\\s","_")}.processed.sumstats.gz"), path(w_ld_chr_ch)

    script:
        """
        ./${ldsc}/munge_sumstats.py \
        --sumstats ${sumstats} \
        --out ${name.replaceAll("\\s","_")}.processed \
        --snp SNP \
        --merge-alleles ${w_ld_chr_ch}/w_hm3.snplist \
        --chunksize 500000
        """
}
