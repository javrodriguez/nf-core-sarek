module load anaconda3/cpu
conda activate env_nf

module add singularity/3.9.8
module add nextflow

nextflow run nf-core/sarek --input ./samplesheet.csv --outdir ./results --genome GATK.GRCh37 -c ./nyu_ultraviolet.config
