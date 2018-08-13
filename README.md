
# ![BABS-ATACSeqPE](https://raw.githubusercontent.com/crickbabs/BABS-ATACSeqPE/master/docs/images/BABS-ATACSeqPE_logo.png)

## Introduction

A [Nextflow](https://www.nextflow.io/) pipeline for processing paired-end Illumina ATACSeq sequencing data.

The pipeline was written by [The Bioinformatics & Biostatistics Group](https://www.crick.ac.uk/research/science-technology-platforms/bioinformatics-and-biostatistics/) at [The Francis Crick Institute](https://www.crick.ac.uk/), London.

## Pipeline summary

1. Raw read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/), [`Fastq Screen`](https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/))
2. Adapter trimming ([`cutadapt`](http://cutadapt.readthedocs.io/en/stable/installation.html))
3. Alignment ([`BWA`](https://sourceforge.net/projects/bio-bwa/files/))
4. Mark duplicates ([`picard`](https://broadinstitute.github.io/picard/))
5. Filtering to remove:
    * reads mapping to mitochondrial DNA ([`SAMtools`](https://sourceforge.net/projects/samtools/files/samtools/))
    * reads mapping to blacklisted regions ([`SAMtools`](https://sourceforge.net/projects/samtools/files/samtools/), [`BEDTools`](https://github.com/arq5x/bedtools2/))
    * reads that are marked as duplicates ([`SAMtools`](https://sourceforge.net/projects/samtools/files/samtools/))
    * reads that arent marked as primary alignments ([`SAMtools`](https://sourceforge.net/projects/samtools/files/samtools/))
    * reads that are unmapped ([`SAMtools`](https://sourceforge.net/projects/samtools/files/samtools/))
    * reads that map to multiple locations ([`SAMtools`](https://sourceforge.net/projects/samtools/files/samtools/))
    * reads containing > 3 mismatches in either read of the pair ([`BAMTools`](https://github.com/pezmaster31/bamtools))
    * reads that have an insert size > 2kb ([`BAMTools`](https://github.com/pezmaster31/bamtools))
    * reads that are soft-clipped ([`BAMTools`](https://github.com/pezmaster31/bamtools))
    * reads that map to different chromosomes ([`Pysam`](http://pysam.readthedocs.io/en/latest/installation.html))
    * reads that arent in FR orientation ([`Pysam`](http://pysam.readthedocs.io/en/latest/installation.html))
    * reads where only one read of the pair fails the above criteria ([`Pysam`](http://pysam.readthedocs.io/en/latest/installation.html))
6. Merge alignments at replicate and sample level ([`picard`](https://broadinstitute.github.io/picard/))
    * Re-mark duplicates ([`picard`](https://broadinstitute.github.io/picard/))
    * Remove duplicate reads ([`SAMtools`](https://sourceforge.net/projects/samtools/files/samtools/))
    * Create normalised bigWig files scaled to 1 million mapped read pairs ([`BEDTools`](https://github.com/arq5x/bedtools2/), [`wigToBigWig`](http://hgdownload.soe.ucsc.edu/admin/exe/))
    * Call broad peaks ([`MACS2`](https://github.com/taoliu/MACS))
    * Annotate peaks relative to gene features ([`HOMER`](http://homer.ucsd.edu/homer/download.html))
    * Merge peaks across all samples and create tabular file to aid in the filtering of the data ([`BEDTools`](https://github.com/arq5x/bedtools2/))
    * Count reads in merged peaks from replicate-level alignments ([`featureCounts`](http://bioinf.wehi.edu.au/featureCounts/))
    * Differential binding analysis, PCA and clustering ([`R`](https://www.r-project.org/), [`DESeq2`](https://bioconductor.org/packages/release/bioc/html/DESeq2.html))
7. Create IGV session file containing bigWig tracks, peaks and differential sites for data visualisation ([`IGV`](https://software.broadinstitute.org/software/igv/)).
8. Collect and present QC at the raw read, alignment and peak-level ([`MultiQC`](http://multiqc.info/), [`R`](https://www.r-project.org/))

## Documentation

The documentation for the pipeline can be found in the `docs/` directory:

1. [Installation](docs/install.md)
2. [Pipeline configuration](docs/config.md)
3. [Reference genome](docs/genome.md)
4. [Design file](docs/design.md)
5. [Running the pipeline](docs/usage.md)
6. [Output and interpretation of results](docs/output.md)
7. [Troubleshooting](docs/troubleshooting.md)

## Pipeline DAG

# ![BABS-ATACSeqPE directed acyclic graph](https://raw.githubusercontent.com/crickbabs/BABS-ATACSeqPE/master/docs/images/BABS-ATACSeqPE_dag.png)

## Credits

The pipeline was written by the [The Bioinformatics & Biostatistics Group](https://www.crick.ac.uk/research/science-technology-platforms/bioinformatics-and-biostatistics/) at [The Francis Crick Institute](https://www.crick.ac.uk/), London.

The pipeline was developed by [Harshil Patel](mailto:harshil.patel@crick.ac.uk), [Philip East](mailto:philip.east@crick.ac.uk) and [Nourdine Bah](mailto:nourdine.bah@crick.ac.uk).

The [NGI-RNAseq](https://github.com/SciLifeLab/NGI-RNAseq) pipeline developed by Phil Ewels was used a template for this pipeline. Many thanks to Phil and the team at SciLifeLab. The help, tips and tricks provided by Paolo Di Tommaso were also invaluable. Thank you!

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details.
