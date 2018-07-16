# PIPELINE & DOCUMENTATION STILL UNDER CONSTRUCTION! ETA AUGUST 2018.   

# ![BABS-ATACSeqPE](https://raw.githubusercontent.com/crickbabs/BABS-ATACSeqPE/master/docs/images/BABS-ATACSeqPE_logo.png)

## Introduction

A [Nextflow](https://www.nextflow.io/) pipeline for processing paired-end Illumina ATACSeq sequencing data.

The pipeline was written by [The Bioinformatics & Biostatistics Group](https://www.crick.ac.uk/research/science-technology-platforms/bioinformatics-and-biostatistics/) at [The Francis Crick Institute](https://www.crick.ac.uk/), London.

## Pipeline summary

1. Raw read QC ([FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/),[Fastq Screen](https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/))
2. Adapter trimming ([cutadapt](http://cutadapt.readthedocs.io/en/stable/installation.html))
3. Alignment ([BWA](https://sourceforge.net/projects/bio-bwa/files/))
4. Mark duplicates ([picard](https://broadinstitute.github.io/picard/))
5. Filtering to remove
    * reads mapping to mitochondrial DNA ([SAMtools](https://sourceforge.net/projects/samtools/files/samtools/))
    * reads mapping to blacklisted regions ([SAMtools](https://sourceforge.net/projects/samtools/files/samtools/))
    * reads that are unmapped ([SAMtools](https://sourceforge.net/projects/samtools/files/samtools/))
    * reads that map to multiple locations ([SAMtools](https://sourceforge.net/projects/samtools/files/samtools/))
    * reads that map to different chromosomes ([Pysam](http://pysam.readthedocs.io/en/latest/installation.html))
    * reads that arent marked as properly paired ([SAMtools](https://sourceforge.net/projects/samtools/files/samtools/))
    * reads that arent marked as primary alignments ([SAMtools](https://sourceforge.net/projects/samtools/files/samtools/))
    * reads containing > 3 mismatches in either read of the pair ([BAMTools](https://github.com/pezmaster31/bamtools))
    * reads that have an insert size > 2kb ([BAMTools](https://github.com/pezmaster31/bamtools))
    * reads that are soft-clipped ([BAMTools](https://github.com/pezmaster31/bamtools))
    * reads where only one read of the pair fails the above criteria ([Pysam](http://pysam.readthedocs.io/en/latest/installation.html))
6. Merge alignments at replicate and sample level ([picard](https://broadinstitute.github.io/picard/))
    * Remove duplicate reads ([SAMtools](https://sourceforge.net/projects/samtools/files/samtools/))
    * Create normalised bigwig files scaled to 1 million mapped read pairs ([BEDTools](https://github.com/arq5x/bedtools2/), [wigToBigWig](http://hgdownload.soe.ucsc.edu/admin/exe/))
    * Call peaks ([MACS2](https://github.com/taoliu/MACS))
    * Annotate peaks ([HOMER](http://homer.ucsd.edu/homer/download.html))
    * Merge peaks across all samples and create tabular file to aid in the filtering of the data ([BEDTools](https://github.com/arq5x/bedtools2/))
    * Count reads in merged peaks at replicate level ([featureCounts](http://bioinf.wehi.edu.au/featureCounts/))
    * Differential binding analysis, PCA and clustering ([R](https://www.r-project.org/),[DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html))
7. Create IGV session file containing bigWig tracks, peaks and differential sites for data visualisation ([IGV](https://software.broadinstitute.org/software/igv/)).
8. Collect and present QC at the raw read, alignment and peak-level in a single report ([MultiQC](http://multiqc.info/),[R](https://www.r-project.org/))

## Documentation

The documentation for the pipeline can be found in the `docs/` directory:

1. [Software requirements](docs/software.md)
2. [Reference genome](docs/genome.md)
3. [Design file](docs/design.md)
4. [Configuration](docs/config.md)
5. [Running the pipeline](docs/usage.md)
6. [Output and interpretation of results](docs/output.md)
7. [Troubleshooting](docs/troubleshooting.md)

## Analysis Pipeline

# ![BABS-ATACSeqPE directed acyclic graph](https://raw.githubusercontent.com/crickbabs/BABS-ATACSeqPE/master/docs/images/BABS-ATACSeqPE_dag.png)

## Credits

The pipeline was written by the [The Bioinformatics & Biostatistics Group](https://www.crick.ac.uk/research/science-technology-platforms/bioinformatics-and-biostatistics/) at [The Francis Crick Institute](https://www.crick.ac.uk/), London.

The pipeline was developed by [Harshil Patel](mailto:harshil.patel@crick.ac.uk), [Philip East](mailto:philip.east@crick.ac.uk) and [Nourdine Bah](mailto:nourdine.bah@crick.ac.uk).

The [NGI-RNAseq](https://github.com/SciLifeLab/NGI-RNAseq) pipeline developed by Phil Ewels was used a template for this pipeline. Many thanks to Phil and the team at SciLifeLab.

<!---
Update dag plot

Given a design file containing sample name to FastQ file mappings the pipeline pre-processes the raw reads ([Cutadapt](http://cutadapt.readthedocs.io/en/stable/guide.html)), aligns the reads to a user-specified genome ([BWA](http://bio-bwa.sourceforge.net/)), marks duplicates ([picard](https://broadinstitute.github.io/picard/)), filters mitochondrial and spurious alignments ([SAMtools](http://samtools.sourceforge.net/), [BamTools](https://github.com/pezmaster31/bamtools)), merges alignments for multiple runs of the same library at the replicate-level ([picard](https://broadinstitute.github.io/picard/)), merges the replicate-level alignments at the sample-level ([picard](https://broadinstitute.github.io/picard/)), generates library-size normalised bigWig tracks ([BEDTools](http://bedtools.readthedocs.io/en/latest/), [wigToBigWig](https://www.encodeproject.org/software/wigtobigwig/)), calls and annotates both broad and narrow peaks ([MACS2](https://github.com/taoliu/MACS), [HOMER](http://homer.ucsd.edu/homer/)), merges peaks across all samples to aid filtering ([BEDTools](http://bedtools.readthedocs.io/en/latest/)), performs differential binding analysis ([featureCounts](http://bioinf.wehi.edu.au/featureCounts/), [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)), creates an IGV session file for visualisation ([IGV](https://software.broadinstitute.org/software/igv/)), and presents various quality-control measures at the raw read, alignment and peak-level ([FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/), [FastQ Screen](https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/), [R](https://www.r-project.org/), [MultiQC](http://multiqc.info/)).
-->
