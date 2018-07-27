
## Load/Install Nextflow

To run nextflow it needs to be available in your path. On an environment module system such as the one at The Francis Crick Institute this can be done by running the following command:

```bash
module load nextflow/0.30.2
```

**You need Nextflow version >= 0.30.2 to run this pipeline.**

See [nextflow.io](https://www.nextflow.io/) for further information on how to install Nextflow.

## Download the pipeline from github

To obtain the pipeline run the following command in a directory where you want to perform the analysis:

```bash
git clone https://github.com/crickbabs/BABS-ATACSeqPE
```

The nextflow pipeline and associated config and executable files will appear in the `BABS-ATACSeqPE` directory. It can now be configured to run on a Linux system of your choice.

```bash
cd BABS-ATACSeqPE
```

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

The following R libraries may need to be installed if unavailable. You can test this be loading the `R` module (if required), typing `R` at the command prompt and attempting to load the packages below e.g. `> library(DESeq2)` and so on. The pipeline assumes the correct R library path is set in order find the installed packages. If not, you can set this in the `.Rprofile` file in the user home directory or add a line which extends the `R` [libPaths](https://stat.ethz.ch/R-manual/R-devel/library/base/html/libPaths.html) in the executable R scripts in the `bin/` directory.

|                                                                         |                                                                         |                                                                                 |
|-------------------------------------------------------------------------|-------------------------------------------------------------------------|---------------------------------------------------------------------------------|
| [optparse](https://cran.r-project.org/web/packages/optparse/index.html) | [scales](https://cran.r-project.org/web/packages/scales/index.html)     | [RColorBrewer](https://cran.r-project.org/web/packages/RColorBrewer/index.html) |
| [ggplot2](https://ggplot2.tidyverse.org/)                               | [pheatmap](https://cran.r-project.org/web/packages/pheatmap/index.html) | [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)       |
| [reshape2](https://cran.r-project.org/web/packages/reshape2/index.html) | [lattice](https://cran.r-project.org/web/packages/lattice/index.html)   | [vsn](https://bioconductor.org/packages/release/bioc/html/vsn.html)             |

### Linux utilities

Standard linux tools including `cut`, `awk`, `sort`, `mv`, `touch`, `echo`, `mkdir`, `paste`, `cp`, `ln`, `grep` are also used throughout the pipeline.

## Software configuration

Depending on where and how you would like to run the pipeline, nextflow needs to be configured to use the software and associated versions.

### Environment modules

If you are running this pipeline at The Francis Crick Institute all of the required software can be loaded with the environment module system on CAMP. A [module config file](https://github.com/crickbabs/BABS-ATACSeqPE/blob/master/conf/babs_modules.config) specifying the software modules to load has been been created for this purpose. If you want to change the versions of software just edit this file before running the pipeline, however, this may require some testing to make sure the nextflow processes requiring the software can still be run with the same command-line parameters.  

By default, the pipeline will be run using the `slurm` job submission system. You will need an account to use the HPC cluster on CAMP if you want to run this at The Francis Crick Institute. If you prefer to run the pipeline locally in batch just replace `executor = 'slurm'` to `executor = 'local'` in the module config file (the pipeline may take a very long time and is only recommended for testing purposes).

When running the pipeline you will have to specify `-profile babs_modules` in order to use this configuration.

### Conda

A [Conda environment file](https://github.com/crickbabs/BABS-ATACSeqPE/blob/master/environment.yaml) is provided with the pipeline for those that wish to run the pipeline on another system without having to painstakingly install all the software and associated dependencies. This will require an internet connection on the command-line and installation of [Anaconda or Miniconda](https://conda.io/docs/user-guide/install/index.html). Nextflow will create the Conda environment by downloading and installing all the required software before execution of the pipeline.

A [conda config file](https://github.com/crickbabs/BABS-ATACSeqPE/blob/master/conf/conda.config) can be found in the `conf` directory. By default, the pipeline will be run using the `slurm` job submission system. You will need an account to use the HPC cluster on CAMP if you want to run this at The Francis Crick Institute. If you prefer to run the pipeline locally in batch just replace `executor = 'slurm'` to `executor = 'local'` in the conda config file (the pipeline may take a very long time and is only recommended for testing purposes). Before the submission of each nextflow process the relevant conda package will need to be available on the command-line to activate the Conda environment. If you are not running the pipeline at The Francis Crick Institute you will need to edit `beforeScript = 'module purge && ml Anaconda2/5.1.0'` in the config file to load/use Conda.

When running the pipeline you will have to specify `-profile conda` in order to use this configuration.

### Local install

If you really need to run the pipeline in serial you can specify `-profile standard` when running the pipeline. This is the most basic profile and will require all the software to be available at the command-line.

<!---

-->
