#!/bin/bash

#SBATCH --time=06:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=10G
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

# Here load needed system tools (Java 1.8 is required, one of singularity or anaconda - python 2.7 are needed,
# depending on the method for dependancy management)

module load jdk/16.0.1
module load openjdk/11.0.2
module load squashfs
module load singularity
module load nextflow/23.09.3-edge

set -f

inputDir="/path"  
w_ld_chr="/path"
sample_counts="/path"
snp_list="/path"
output_folder="/path" 

NXF_VER=23.09.3-edge nextflow run main.nf  \
--inputDir ${inputDir} \
--SnpRefFile ${snp_list} \
--p_thresh 5e-8 \
--OutputDir ${output_folder} \
--CaseControlFile ${sample_counts} \
--wLdChr ${w_ld_chr} \
--LdscDir "/path" \
-profile slurm,singularity \
-resume \
