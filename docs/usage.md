
## Run pipeline

When you have [installed all the necessary software](https://github.com/crickbabs/BABS-ATACSeqPE/blob/master/docs/software.md), [prepared the reference genome](https://github.com/crickbabs/BABS-ATACSeqPE/blob/master/docs/genome.md), [created a design file](https://github.com/crickbabs/BABS-ATACSeqPE/blob/master/docs/design.md), and [configured the pipeline](https://github.com/crickbabs/BABS-ATACSeqPE/blob/master/docs/configuration/local.md) you can run it with the command below:

```bash
nextflow run main.nf --design <DESIGN_FILE> --genome <GENOME_NAME> -profile <PROFILE_NAME>
```

You can also specify parameters at the command-line to override those provided in the config files. A list of these can be obtained with the command below:

```bash
nextflow run main.nf --help
```

<!---
Add troubleshooting section here?
-->
