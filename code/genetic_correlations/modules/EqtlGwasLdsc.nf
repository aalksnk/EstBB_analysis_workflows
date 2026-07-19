#!/bin/bash nextflow


process PrepareGwas {

    container 'quay.io/urmovosa/ld_overlap_image:v0.2'
   // tag "${file.baseName}"
    label 'R' 
    
   // scratch true

   // publishDir "${params.OutputDir}", mode: 'copy', overwrite: true, pattern: "*.log"
    //container 'quay.io/urmovosa/pqtlvseqtl:v0.1'

    tag "${file.baseName}"

    input:
        tuple path(file), path(snp_list), path(sample_count_file), path(w_ld_chr)

    output:
        tuple file("*_processed.txt"), env(phenoname)

    script:
        """
        ##############################################################
        # Convert to LDSC .txt format and filter to HapMap3 variants #
        ##############################################################
        echo "GWAS file: ${file}"
        echo "SNP list: ${snp_list}"
        echo "LDSC reference: ${w_ld_chr}/w_hm3.snplist"

        Rscript --vanilla ${baseDir}/bin/PrepareGwasSumstats.R \
        ${file} \
        ${snp_list} \
        ${w_ld_chr}/w_hm3.snplist \
        ${sample_count_file}

        phenoname=\$(ls *_processed* | sed 's/^results_concat_\\(.*\\)\\(\\.parquet\\)*\\.snappy\$/\\1/')
        """
}

process MungeGwas {
    //scratch true
    container 'manninglab/ldsc'

    tag "$phenotype"

    input:
        tuple path(file), val(phenotype), path(ldsc), path(w_ld_chr_ch)

    output:
        tuple val(phenotype), path("*.sumstats.gz")

    script:
        """
        #!/bin/bash
        # Assign Nextflow input to a Bash variable
        phenotype="${phenotype}"

        # Strip extensions from phenotype name (remove .processed.txt)
        clean_phenotype=\$(basename "\$phenotype" | sed 's/.processed.txt//')

        echo "Using clean phenotype name: \$clean_phenotype"

        /gpfs/space/GI/GV/Projects/SAMPO/gencorrs/input/munge_sumstats.py \
        --sumstats ${file} \
        --out \${clean_phenotype} \
        --snp SNP \
        --merge-alleles ${w_ld_chr_ch}/w_hm3.snplist \
        """
}



process CollectGwas {
    scratch true

    input:
        path(gwas)

    output:
        path("gwas_folder")

    script:
        """
        mkdir -p gwas_folder
        mv ${gwas} gwas_folder/

        echo "=== GWAS files collected ==="
        ls -lh gwas_folder
        """
}





process Ldsc {
    //scratch true
    container 'manninglab/ldsc'
    
    tag "$phenotype"

    input:
        tuple val(phenotype), path(gwas_sumstats), path(all_gwas_sumstats), path(ldsc), path(ldsc_ref)

    output:
        tuple path("*vsAll.txt"), path("*_heritability.txt")

    script:
        """
        ###################################
        # Get the string of ALL GWAS sumstats #
        ###################################
       
        gwas_sumstats_string=\$(find -L ${all_gwas_sumstats} -name '*sumstats.gz' | tr '\\n' ',' | sed 's/,\$//')

        echo "Phenotype: ${phenotype}"
        echo "GWAS Sumstats File: ${gwas_sumstats}"
        echo "All GWAS Sumstats: \${gwas_sumstats_string}"

        if [ -z "\${gwas_sumstats_string}" ]; then
            echo "Error: No sumstats files found in ${all_gwas_sumstats}"
            exit 1
        fi


        ################################
        # Run LDSC genetic correlation #
        ################################

        /gpfs/space/GI/GV/Projects/SAMPO/gencorrs/input/ldsc.py \
            --rg ${gwas_sumstats},\${gwas_sumstats_string} \
            --ref-ld-chr ${ldsc_ref}/ \
            --w-ld-chr ${ldsc_ref}/ \
            --out ${phenotype}vsAll
        
        # rm ${phenotype}.sumstats.gz

        ##########################################
        # Extract results from the LDSC log file #
        ##########################################

        parse_ldsc_heritability.sh ${phenotype}vsAll.log > ${phenotype}_heritability.txt
        parse_ldsc_rg.sh ${phenotype}vsAll.log > ${phenotype}vsAll.txt
        """
}



workflow PREPAREGWAS {
    take:
        data

    main:
        loci_ch = PrepareGwas(data)
        
    emit:
        loci_ch

}

workflow MUNGEGWAS {
    take:
        input_ch

    main:
        munged = MungeGwas(input_ch)
    
    emit:
        munged
}


workflow COLLECTGWAS {
    take:
        data

    main:
        loci_ch = CollectGwas(data)
        
    emit:
        loci_ch

}


workflow LDSC {
    take:
        input_ch

    main:
        output_ch = Ldsc(input_ch)
    
    emit:
        output_ch
}



