

set -f

nextflow_path=/Users/urmovosa/Documents/temp_tools/

inputDir=/Users/urmovosa/Documents/projects/2023/sampo_pipelines/ldoverlapnf/tests/data/sumstats/
outputDir=/Users/urmovosa/Documents/projects/2023/sampo_pipelines/HeritabilitySampoNF/tests/output/

NXF_VER=21.10.6 ${nextflow_path}/nextflow run main.nf  \
--inputDir ${inputDir} \
--wLdChr /Users/urmovosa/Documents/projects/2023/sampo_pipelines/HeritabilitySampoNF/tests/eur_w_ld_chr/ \
--SnpRefFile /Users/urmovosa/Documents/projects/2023/sampo_pipelines/HeritabilitySampoNF/tests/snp_reference.parquet \
--CaseControlFile /Users/urmovosa/Documents/projects/2023/sampo_pipelines/HeritabilitySampoNF/tests/phenos_cases_controls \
--LdscDir /Users/urmovosa/Documents/projects/2023/sampo_pipelines/HeritabilitySampoNF/tests/ldsc \
-profile docker \
-resume 