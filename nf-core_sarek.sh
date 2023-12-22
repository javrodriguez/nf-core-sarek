#!/bin/bash
#SBATCH -J sarek
#SBATCH --mem=5gb
#SBATCH --time=72:00:00
#SBATCH --output=logs-sarek/%J.log
#SBATCH --error=logs-sarek/%J.out

module add singularity/3.9.8
module add nextflow

nextflow run nf-core/sarek --input ./samplesheet.csv --outdir ./results --genome GATK.GRCh37 -c ./nyu_ultraviolet.config
