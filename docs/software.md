
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

If you are running this pipeline at The Francis Crick Institute all of the software below can be loaded with the environment module system on CAMP. The pipeline was tested with the software versions listed below but should work with other releases. This may require some testing to make sure the parameters are still applicable.  

|                                                                                  |        |                                                                       |                |
|----------------------------------------------------------------------------------|--------|-----------------------------------------------------------------------|----------------|
| [nextflow](https://www.nextflow.io/)                                             | 0.30.2 | [SAMtools](https://sourceforge.net/projects/samtools/files/samtools/) | 1.3.1          |
| [R](https://www.r-project.org/)                                                  | 3.3.1  | [BamTools](https://github.com/pezmaster31/bamtools)                   | 2.4.0          |
| [Python](https://www.python.org/downloads/)                                      | 2.7.12 | [BEDTools](https://github.com/arq5x/bedtools2/)                       | 2.26.0         |
| [Java](https://java.com/en/download/)                                            | 1.8.0 | [Kent_tools](http://hgdownload.soe.ucsc.edu/admin/exe/)               | 20161115       |
| [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)             | 0.11.5 | [MACS2](https://github.com/taoliu/MACS)                               | 2.1.1.20160309 |
| [FastQ Screen](https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/) | 0.9.3  | [HOMER](http://homer.ucsd.edu/homer/download.html)                    | 4.8            |
| [cutadapt](http://cutadapt.readthedocs.io/en/stable/installation.html)           | 1.9.1  | [featureCounts](http://bioinf.wehi.edu.au/featureCounts/)             | 1.5.1          |
| [BWA](https://sourceforge.net/projects/bio-bwa/files/)                           | 0.6.2  | [MultiQC](http://multiqc.info/)                                       | 1.3            |
| [picard](https://broadinstitute.github.io/picard/)                               | 2.1.1  | [Pysam](http://pysam.readthedocs.io/en/latest/installation.html)      | 0.9.0          |

### R libraries

The following R libraries may need to be installed if unavailable. You can test this be loading the `R` module (if required), typing `R` at the command prompt and attempting to load the packages below e.g. `> library(DESeq2)` and so on. The pipeline assumes the correct R library path is set in order find the installed packages. If not, you can set this in the `.Rprofile` file in the user home directory or add a line which extends the `R` [libPaths](https://stat.ethz.ch/R-manual/R-devel/library/base/html/libPaths.html) in the executable R scripts in the `bin/` directory.

[DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)    
[argparse](https://cran.r-project.org/web/packages/argparse/index.html)  
[vsn](https://bioconductor.org/packages/release/bioc/html/vsn.html)  
[ggplot2](https://ggplot2.tidyverse.org/)  
[reshape2](https://cran.r-project.org/web/packages/reshape2/index.html)  
[scales](https://cran.r-project.org/web/packages/scales/index.html)  
[RColorBrewer](https://cran.r-project.org/web/packages/RColorBrewer/index.html)  
[pheatmap](https://cran.r-project.org/web/packages/pheatmap/index.html)  
[lattice](https://cran.r-project.org/web/packages/lattice/index.html)  

### Linux utilities

Standard linux tools including `cut`, `awk`, `sort`, `mv`, `touch`, `echo`, `mkdir`, `paste`, `cp`, `ln`, `grep` are also used throughout the pipeline.

## Software configuration

Once installed, nextflow needs to be configured to use the software and associated versions.

### Environment modules

See [BABS module definitions](https://github.com/crickbabs/BABS-ATACSeqPE/blob/master/config/babs.config).

### Conda

See [Conda environment file](https://github.com/crickbabs/BABS-ATACSeqPE/blob/master/environment.yml).

### Local install

Local will rely on software being in path.


<!---
Add information on how to download file without internet connection see NGI-RNASeq
Add section on how to edit config files

## conda config --add envs_dirs /camp/stp/babs/working/patelh/code/conda/envs/
## conda config --add pkgs_dirs /camp/stp/babs/working/patelh/code/conda/pkgs/

## Download the environment.yml file
# curl https://raw.githubusercontent.com/crickbabs/BABS-ATACSeqPE/master/environment.yml -o environment.yml
#
## Load the Anaconda module
# ml purge && ml Anaconda2/5.1.0
#
## Create a new conda environment using it
# conda env create -f environment.yml
#
## Activate the new conda environment
# source activate BABS-ATACSeqPE
#
## Deactivate the conda environment
# source deactivate

##
## ml purge && ml Anaconda2/5.1.0
## conda env create -f BABS-ATACSeqPE_environment.yml
## conda remove --name BABS-ATACSeqPE --all && conda info --envs

-->
