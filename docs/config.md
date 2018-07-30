## Software configuration

Depending on where and how you would like to run the pipeline, nextflow needs to be configured to use the [software](https://github.com/crickbabs/BABS-ATACSeqPE/blob/master/docs/config.md#software-requirements) and associated versions.

### Environment modules

If you are running this pipeline at The Francis Crick Institute all of the required software can be loaded with the environment module system on CAMP. A [module config file](https://github.com/crickbabs/BABS-ATACSeqPE/blob/master/conf/babs_modules.config) specifying the software modules to load has been been created for this purpose. If you want to change the versions of software just edit this file before running the pipeline, however, this may require some testing to make sure the nextflow processes requiring the software can still be run with the same command-line parameters.  

By default, the pipeline will be executed using the `slurm` job submission system. You will need an account to use the HPC cluster on CAMP if you want to run this at The Francis Crick Institute. If you prefer to run the pipeline locally in serial just replace `executor = 'slurm'` to `executor = 'local'` in the module config file. However, the pipeline may take a very long time to complete and this is only recommended for testing purposes.

When running the pipeline just specify `-profile babs_modules` in order to use this configuration.

### Conda

A [Conda environment file](https://github.com/crickbabs/BABS-ATACSeqPE/blob/master/environment.yaml) is provided with the pipeline for those that wish to run the pipeline on another system without having to painstakingly install all the software and associated dependencies. This will require an internet connection on the command-line, and installation of [Anaconda or Miniconda](https://conda.io/docs/user-guide/install/index.html). Nextflow will rather amazingly create the Conda environment by downloading and installing all the required software before execution of the pipeline. This could take up to 45 minutes.

A [Conda config file](https://github.com/crickbabs/BABS-ATACSeqPE/blob/master/conf/conda.config) can be found in the `conf/` directory. By default, the pipeline will be executed using the `slurm` job submission system. You will need an account to use the HPC cluster on CAMP if you want to run this at The Francis Crick Institute. If you prefer to run the pipeline locally in serial just replace `executor = 'slurm'` to `executor = 'local'` in the conda config file. However, the pipeline may take a very long time to complete and this is only recommended for testing purposes.  

Before the submission of each nextflow process the `conda` command will need to be available on the command-line to activate the Conda environment. If you are not running the pipeline at The Francis Crick Institute you will need to edit `beforeScript = 'module purge && ml Anaconda2/5.1.0'` in the config file to load/use Conda.

When running the pipeline just specify `-profile conda` in order to use this configuration.

### Local install

If you really need to run the pipeline in serial you can specify `-profile standard` when running the pipeline. This is the most basic profile and will require all the software to be available at the command-line.

## Config files

Nextflow offers the flexibility to store parameters required by the pipeline in config files. This offers the flexibility to create custom implementations of the pipeline for running on various compute platforms, and to define pipeline/process-specific parameters in a single file.

| Path                                                                                                                   | Description                                                                                                                                                                                                                                                                                                                          |
| -----------------------------------------------------------------------------------------------------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| [`nextflow.config`](https://github.com/crickbabs/BABS-ATACSeqPE/tree/master/nextflow.config)                           | Config file used for the definition of nextflow profiles and other standard parameters. If you decide to use the `genomes.config` to pre-define the genome parameters for the pipeline, you will only need to change the `genome_base` parameter in this file to reflect the top-level directory where your genome files are stored. |
|                                                                                                                        |                                                                                                                                                                                                                                                                                                                                      |
| [`conf/genomes.config`](https://github.com/crickbabs/BABS-ATACSeqPE/tree/master/conf/genomes.config)                   | Config file used to store genome parameters for multiple assemblies. You will need to customise this file to use the specified genome parameters without needing to specify them individually at the command-line. This also increases the portability of the pipeline for use on other systems where the files are exactly the same but may differ in the path to the top-level directory. |
| [`conf/base.config `](https://github.com/crickbabs/BABS-ATACSeqPE/tree/master/conf/base.config)                        | Config file defining resources used by each process.                                                                                                                                                                                                                                                                                 |
| [`conf/conda.config`](https://github.com/crickbabs/BABS-ATACSeqPE/tree/master/conf/conda.config)                       | Config file specifying parameters required to build Conda environment. You will need to customise the `beforeScript` parameter in this file to load the relevant Conda distribution available on your system.                                                                                                                                                                                                                                                                |
| [`conf/babs_modules.config`](https://github.com/crickbabs/BABS-ATACSeqPE/tree/master/conf/babs_modules.config)         | Config file for environment modules that can be used at The Francis Crick Institute. You can customise this file to use specific versions of software for each process in the pipeline or you can use it as a template to build your own for use outside The Francis Crick Institute.                                                                                                                                                                                                                                                 |
| [`conf/fastq_screen.conf.txt`](https://github.com/crickbabs/BABS-ATACSeqPE/tree/master/conf/fastq_screen.conf.txt)     | Config file specifying paths to bowtie2 indexed genomes that will be used for contaminant screening. If you want the pipeline to run `fastq_screen` this file will need to be customised to reflect the paths to the pre-created bowtie2 indices with which to run FastQ Screen.                                                                                                                                                                                                                                 |
| [`conf/multiqc_config.yaml`](https://github.com/crickbabs/BABS-ATACSeqPE/tree/master/conf/multiqc_config.yaml)         | Config file for generating customised MultiQC report for BABS-ATACSeqPE pipeline. You wont need to customise this file unless you want to change which aspects of the QC are reported by MultiQC.                                                                                                                                                                                                                                                    |
| [`conf/bamtools_filter_pe.json`](https://github.com/crickbabs/BABS-ATACSeqPE/tree/master/conf/bamtools_filter_pe.json) | Config file used by BamTools for filtering. You wont need to customise this file unless you want to amend specific aspects of the read filtering criteria.                                                                                                                                                                                                                                                                                         |

## Software requirements

The software below is required to run the pipeline:

|                                                                                  |                                                                       |                                                                  |
|----------------------------------------------------------------------------------|-----------------------------------------------------------------------|------------------------------------------------------------------|
| [nextflow](https://www.nextflow.io/)                                             | [cutadapt](http://cutadapt.readthedocs.io/en/stable/guide.html)       | [Kent_tools](http://hgdownload.soe.ucsc.edu/admin/exe/)          |
| [R](https://www.r-project.org/)                                                  | [BWA](https://sourceforge.net/projects/bio-bwa/files/)                | [MACS2](https://github.com/taoliu/MACS)                          |
| [Python](https://www.python.org/downloads/)                                      | [picard](https://broadinstitute.github.io/picard/)                    | [HOMER](http://homer.ucsd.edu/homer/download.html)               |
| [Java](https://java.com/en/download/)                                            | [SAMtools](https://sourceforge.net/projects/samtools/files/samtools/) | [featureCounts](http://bioinf.wehi.edu.au/featureCounts/)        |
| [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)             | [BamTools](https://github.com/pezmaster31/bamtools)                   | [Pysam](http://pysam.readthedocs.io/en/latest/installation.html) |
| [FastQ Screen](https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/) | [BEDTools](https://github.com/arq5x/bedtools2/)                       | [MultiQC](http://multiqc.info/)                                  |

### R libraries

The following R libraries may need to be installed if unavailable. You can test this be loading the `R` module (if required), typing `R` at the command prompt and attempting to load the packages below e.g. `> library(optparse)` and so on. The pipeline assumes the correct R library path is set in order find the installed packages. If not, you can set this in the `.Rprofile` file in the user home directory or add a line which extends the `R` [libPaths](https://stat.ethz.ch/R-manual/R-devel/library/base/html/libPaths.html) in the executable R scripts in the `bin/` directory.

|                                                                         |                                                                         |                                                                                 |
|-------------------------------------------------------------------------|-------------------------------------------------------------------------|---------------------------------------------------------------------------------|
| [optparse](https://cran.r-project.org/web/packages/optparse/index.html) | [scales](https://cran.r-project.org/web/packages/scales/index.html)     | [RColorBrewer](https://cran.r-project.org/web/packages/RColorBrewer/index.html) |
| [ggplot2](https://ggplot2.tidyverse.org/)                               | [pheatmap](https://cran.r-project.org/web/packages/pheatmap/index.html) | [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)       |
| [reshape2](https://cran.r-project.org/web/packages/reshape2/index.html) | [lattice](https://cran.r-project.org/web/packages/lattice/index.html)   | [vsn](https://bioconductor.org/packages/release/bioc/html/vsn.html)             |

### Linux utilities

Standard linux tools including `cut`, `awk`, `sort`, `mv`, `touch`, `echo`, `mkdir`, `paste`, `cp`, `ln`, `grep` are also used throughout the pipeline.

<!---
Add information on how to customise each of these files to get pipeline running see NGI-RNASeq
-->
