#!/usr/bin/env nextflow

/*
 * enables modules
 */
nextflow.enable.dsl = 2


def helpmessage() {

log.info"""

EqtlGwasLdsc v${workflow.manifest.version}"
===========================================================
Pipeline for running LDSC genetic correlation analyses analyses (https://www.nature.com/articles/s41467-020-20885-8) between GWAS summary statistics. Based on the pipeline https://github.com/eQTLGen/eQTLGenEqtlGwasLdsc/tree/main. 

Usage:

nextflow run main.nf 
--eqtl_files \
--allele_info \
--ldsc_ref \
--OutputDir

Mandatory arguments:
--eqtl_files                eQTLGen parquet dataset.
--gwas_files                GWAS parquet dataset. Needs to be harmonised to follow same format as eQTL dataset.
--allele_info               Parquet file with alleles and SNP positions for eQTL dataset.
--snplist                   File with HapMap3 variants matching to eQTL dataset.
--ldsc_folder               Folder with LDSC scripts (https://github.com/bulik/ldsc).
--ldsc_ref                  Folder with LDSC reference files, as shared by LDSC authors.
--remove_hla                Whether to remove HLA region from the input eQTL file (yes/no). Default is no.

Optional arguments:
--OutputDir                 Output directory. Defaults to "results".

""".stripIndent()

}

if (params.help){
    helpmessage()
    exit 0
}

// Default parameters
params.OutputDir = 'results'
params.remove_hla = 'no'


//Show parameter values
log.info """=======================================================
GWAS genetic correlation pipeline v${workflow.manifest.version}"
======================================================="""
def summary = [:]
summary['Pipeline Version']                         = workflow.manifest.version
summary['Current user']                             = "$USER"
summary['Current home']                             = "$HOME"
summary['Current path']                             = "$PWD"
summary['Working dir']                              = workflow.workDir
summary['Script dir']                               = workflow.projectDir
summary['Config Profile']                           = workflow.profile
summary['Container Engine']                         = workflow.containerEngine
if(workflow.containerEngine) summary['Container']   = workflow.container
summary['Output directory']                         = params.OutputDir
summary['GWAS folder']                              = params.inputDir
//summary['GWAS manifest file']                       = params.gwas_manifest
summary['LDSC folder']                              = params.LdscDir
summary['LDSC reference']                           = params.wLdChr





// import modules
include { PREPAREGWAS; MUNGEGWAS; LDSC; COLLECTGWAS; PrepareGwas; MungeGwas; Ldsc; CollectGwas } from './modules/EqtlGwasLdsc.nf'

log.info summary.collect { k,v -> "${k.padRight(21)}: $v" }.join("\n")
log.info "======================================================="

// Define argument channels

Channel
    .fromPath("${params.inputDir}/*.parquet.snappy", checkIfExists: true)
    .ifEmpty { exit 1, "Sumstats directory is empty!" }
    .set { input_files_ch }
Channel
    .fromPath("${params.SnpRefFile}/*", checkIfExists: true)
    .ifEmpty { exit 1, "SNP reference file is missing!" }
    .set { snp_ref_ch }
Channel
    .fromPath(params.wLdChr, checkIfExists: true, type: 'dir')
    .ifEmpty { exit 1, "Weighted LD directory is empty or files are missing!" }
    .set { w_ld_chr_ch }
Channel
    .fromPath(params.LdscDir, checkIfExists: true, type: 'dir')
    .ifEmpty { exit 1, "LDSC directory is empty or files are missing!" }
    .set { ldsc_ch }
Channel
    .fromPath(params.CaseControlFile, checkIfExists: true)
    .ifEmpty { exit 1, "Case-control file is missing!" }
    .set { case_control_ch }

//gwas_input_ch = input_files_ch.combine(snp_ref_ch).combine(ldsc_ch).combine(w_ld_chr_ch)
/*workflow {
        PREPAREGWAS(input_files_ch.combine(snp_ref_ch).combine(case_control_ch).combine(w_ld_chr_ch)) 
        MUNGEGWAS(PREPAREGWAS.out.combine(ldsc_ch).combine(w_ld_chr_ch))
        
        collectgwas_input_ch = MUNGEGWAS.out.map { it[1] }.collect()
        COLLECTGWAS(collectgwas_input_ch)

        LDSC(MUNGEGWAS.out.combine(COLLECTGWAS.out).combine(ldsc_ch).combine(w_ld_chr_ch))

        LDSC.out.map { it[0] }.collectFile(name: 'EqtlGwasLdscResults.txt', keepHeader: true, sort: true, storeDir: "${params.OutputDir}")
        LDSC.out.map { it[1] }.collectFile(name: 'EqtlHeritabilityLdscResults.txt', keepHeader: true, sort: true, storeDir: "${params.OutputDir}")

        }


workflow.onComplete {
    println ( workflow.success ? "Pipeline finished!" : "Something crashed...debug!" )
}*/


workflow {
    prepare_ch = PREPAREGWAS(input_files_ch.combine(snp_ref_ch).combine(case_control_ch).combine(w_ld_chr_ch))

    munged_ch = MUNGEGWAS(prepare_ch.combine(ldsc_ch).combine(w_ld_chr_ch))

    // Collect all GWAS sumstats into one folder
    collectgwas_input_ch = munged_ch.map { it[1] }.collect()

// Run CollectGwas to create the folder
    collected_gwas_ch = COLLECTGWAS(collectgwas_input_ch)

// Pass collected_gwas_ch (the folder) directly to LDSC
    LDSC(munged_ch.combine(collected_gwas_ch).combine(ldsc_ch).combine(w_ld_chr_ch))


    LDSC.out.map { it[0] }.collectFile(name: 'EqtlGwasLdscResults.txt', keepHeader: true, sort: true, storeDir: "${params.OutputDir}")
    LDSC.out.map { it[1] }.collectFile(name: 'EqtlHeritabilityLdscResults.txt', keepHeader: true, sort: true, storeDir: "${params.OutputDir}")
}



