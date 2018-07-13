
## Creating a sample design file

You will need to create a design file with information about the samples in your experiment before running the pipeline. It has to be a comma-separated file with 5 columns, and a header row as shown in the example below:

```bash
sample,replicate,run,fastq_1,fastq_2
control,1,1,control_R1_L1_1.fastq.gz,control_R1_L1_2.fastq.gz
control,2,1,control_R2_L1_1.fastq.gz,control_R2_L1_2.fastq.gz
treatment,1,1,treatment_R1_L1_1.fastq.gz,treatment_R1_L1_2.fastq.gz
treatment,2,1,treatment_R2_L1_1.fastq.gz,treatment_R2_L1_2.fastq.gz
treatment,2,2,treatment_R2_L2_1.fastq.gz,treatment_R2_L2_2.fastq.gz
```

| Column      | Description                                                                                                                               |
|-------------|-------------------------------------------------------------------------------------------------------------------------------------------|
| `sample`    | Group identifier for sample                                                                                                               |
| `replicate` | Integer representing replicate number                                                                                                     |
| `run`       | Integer representing the number of times the same library has been sequenced. This will be used later for merging at the replicate-level. |
| `fastq_1`   | Path to FastQ file for read 1. File has to be zipped and have the extension ".fastq.gz".                                                  |
| `fastq_2`   | Path to FastQ file for read 2. File has to be zipped and have the extension ".fastq.gz".                                                  |
