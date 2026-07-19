#!/bin/bash nextflow 

process ConsolidateResults {
    container 'quay.io/cawarmerdam/ldsc:v0.3'

    input:
        path(logs)
        

    output:
        path("final_results_table_*.tsv")

    script:
    """
    # Initialize CSV files with headers
    echo 'Phenotype,h2,se_h2,intercept,se_intercept,ratio,mean_chi2,lambda_gc,nr_snps,nr_snps_removed_chi2_filter' > final_results_table_liability.tsv
    echo 'Phenotype,h2,se_h2,intercept,se_intercept,ratio,mean_chi2,lambda_gc,nr_snps,nr_snps_removed_chi2_filter' > final_results_table_observed.tsv

    # Process Liability Logs
    for log in \$(ls | grep "liability.log")
    do
        phenotype=\$(basename \${log} '_heritability_liability.log')
        h2=\$(grep 'Total Liability scale h2:' \${log} | cut -d ' ' -f5)
        se_h2=\$(grep 'Total Liability scale h2:' \${log} | cut -d ' ' -f6 | sed 's/[()]//g')
        intercept=\$(grep 'Intercept:' \${log} | cut -d ' ' -f2)
        se_intercept=\$(grep 'Intercept:' \${log} | cut -d ' ' -f3 | sed 's/[()]//g')
        ratio=\$(grep 'Ratio' \${log} | cut -d ' ' -f2 | sed 's/[<]/<0/g')
        mean_chi2=\$(grep 'Mean Chi' \${log} | cut -d ' ' -f3)
        lambda_gc=\$(grep 'Lambda GC:' \${log} | cut -d ' ' -f3)
        nr_snps=\$(grep 'After merging with regression SNP LD' \${log} | cut -d ' ' -f7)
        nr_snps_removed_chi2_filter=\$(grep ' SNPs with chi' \${log} | cut -d ' ' -f2)
        echo "\${phenotype},\${h2},\${se_h2},\${intercept},\${se_intercept},\${ratio},\${mean_chi2},\${lambda_gc},\${nr_snps},\${nr_snps_removed_chi2_filter}" >> final_results_table_liability.tsv
    done

    # Process Observed Logs
    for log in \$(ls | grep "observed.log")
    do
        phenotype=\$(basename \${log} '_heritability_observed.log')
        h2=\$(grep 'Total Observed scale h2:' \${log} | cut -d ' ' -f5)
        se_h2=\$(grep 'Total Observed scale h2:' \${log} | cut -d ' ' -f6 | sed 's/[()]//g')
        intercept=\$(grep 'Intercept:' \${log} | cut -d ' ' -f2)
        se_intercept=\$(grep 'Intercept:' \${log} | cut -d ' ' -f3 | sed 's/[()]//g')
        ratio=\$(grep 'Ratio' \${log} | cut -d ' ' -f2 | sed 's/[<]/<0/g')
        mean_chi2=\$(grep 'Mean Chi' \${log} | cut -d ' ' -f3)
        lambda_gc=\$(grep 'Lambda GC:' \${log} | cut -d ' ' -f3)
        nr_snps=\$(grep 'After merging with regression SNP LD' \${log} | cut -d ' ' -f7)
        nr_snps_removed_chi2_filter=\$(grep ' SNPs with chi' \${log} | cut -d ' ' -f2)
        echo "\${phenotype},\${h2},\${se_h2},\${intercept},\${se_intercept},\${ratio},\${mean_chi2},\${lambda_gc},\${nr_snps},\${nr_snps_removed_chi2_filter}" >> final_results_table_observed.tsv
    done
    """
}
