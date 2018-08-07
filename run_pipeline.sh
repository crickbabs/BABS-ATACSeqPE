#!/bin/bash

#### (1) Create design file and place in same directory as this script

#### (2) Download pipeline
git clone https://github.com/crickbabs/BABS-ATACSeqPE
cd BABS-ATACSeqPE

#### (3) Set correct MODULEPATH (Uncomment if running at The Francis Crick Institute)
# source ./conf/babs_modules_path.sh
# ml purge

#### (4) Load Nextflow (Uncomment if running at The Francis Crick Institute or change as required to run on your system)
# ml nextflow/0.30.2

#### (5) Customise 'conf/genomes.config' to specify custom paths to genomes

#### (6) Customise other aspects of the pipeline config if required

###########################################

#### (7) Run pipeline with software specified in 'conf/babs_modules.config' via SLURM (Uncomment if running at The Francis Crick Institute)
# nextflow run main.nf --design ../design.csv --genome hg19 -profile babs_modules --outdir ../results/

################# OR ######################

#### (7) Run pipeline with Conda using software specified in 'environment.yaml' via SLURM. Need to load Conda first.
# ml Anaconda2/5.1.0			## Uncomment if running at The Francis Crick Institute or change as required to run on your system
# nextflow run main.nf --design ../design.csv --genome hg19 -profile conda --outdir ../results/

