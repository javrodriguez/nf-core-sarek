#!/bin/bash
#SBATCH -J sarek
#SBATCH --mem=20gb
#SBATCH --time=48:00:00
#SBATCH --output=logs-sarek/%J.log
#SBATCH --error=logs-sarek/%J.out

module load anaconda3/cpu
conda activate env_nf

module add singularity/3.9.8
module add nextflow

nextflow run nf-core/sarek --input ./samplesheet.csv --outdir ./results --genome GATK.GRCh37 -c ./nyu_ultraviolet.config
