#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
==============================================================================
    B A B S - A T A C S E Q  P A I R E D - E N D   B E S T    P R A C T I C E
==============================================================================
 ATACSeq Analysis Pipeline For Paired-End Illumina Samples.
 Started 11th May 2018.
 #### Homepage / Documentation
 https://github.com/crickbabs/BABS-ATACSeqPE
 #### Authors
 Harshil Patel <harshil.patel@crick.ac.uk>
 Philip East   <philip.east@crick.ac.uk>
 Nourdine Bah  <nourdine.bah@crick.ac.uk>
-------------------------------------------------------------------------------
*/

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                             PARAMETERS                              -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

def helpMessage() {
    log.info"""
    ======================================================
    BABS ATACSeq Paired-End Pipeline v${params.version}
    ======================================================

    Usage:
    The typical command for running the pipeline is as follows:

    nextflow run main.nf --design design.csv --genome hg19 -profile babs

    Mandatory arguments:
      --design                      Comma separted file containing information about the samples (see README)
      --genome                      Genome shortname (e.g. hg19)
      -profile                      Hardware config to use e.g. babs

    References:                     If not specified in the configuration file or you wish to overwrite any of the references
      --fasta                       Path to Fasta reference file containing all chromosomes/contigs
      --chrom_size_file             Path to tab-delimited file with two columns i.e. chrom\tchrom_size
      --gtf                         Path to GTF file
      --bwa_index                   Path to BWA index
      --mito_name                   Name of Mitochondrial chomosome in genome fasta (e.g. chrM)
      --genome_mask                 BED file specifying target regions for analysis. This excludes problematic genomic loci
      --macs_genome_size            Effective genome size parameter required by MACS2

    Trimming options:
      --adapter                     3' adapter sequence trimmed by cutadapt (default: CTGTCTCTTATA)

    Output options:
      --outdir                      The output directory where the results will be saved (default: './results')

    Other options:
      -name                         Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.
    """.stripIndent()
}

// SHOW HELP MESSAGE
params.help = false
if (params.help){
    helpMessage()
    exit 0
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                SET UP CONFIGURATION VARIABLES                       -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////
/* --         DEFAULT PARAMETER VALUES         -- */
////////////////////////////////////////////////////

params.design = false
params.genome = false
params.profile = false
params.fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
params.chrom_size_file = params.genome ? params.genomes[ params.genome ].chrom_size_file ?: false : false
params.gtf = params.genome ? params.genomes[ params.genome ].gtf ?: false : false
params.mito_name = params.genome ? params.genomes[ params.genome ].mito_name ?: false : false
params.bwa_index = params.genome ? params.genomes[ params.genome ].bwa_index ?: false : false
params.genome_mask = params.genome ? params.genomes[ params.genome ].genome_mask ?: false : false
params.macs_genome_size = params.genome ? params.genomes[ params.genome ].macs_genome_size ?: false : false
params.adapter = false
params.outdir = './results'
params.outdir_abspath = new File(params.outdir).getCanonicalPath().toString()
params.name = false
params.project = false
params.multiqc_config = "$baseDir/conf/multiqc_config.yaml"
params.fastqscreen_config = "$baseDir/conf/fastq_screen.conf.txt"
params.bamtools_filter_pe_config = "$baseDir/conf/bamtools_filter_pe.json"

// PRESET ADAPTER TRIMMING OPTION
adapter_seq = 'CTGTCTCTTATA'
if (params.adapter){
    adapter_seq = params.adapter
}

// SPECIFY THE RUN NAME
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
    custom_runName = workflow.runName
}

output_docs = file("$baseDir/docs/output.md")

////////////////////////////////////////////////////
/* --          CONFIG FILES                    -- */
////////////////////////////////////////////////////

multiqc_config = file(params.multiqc_config)
if( !multiqc_config.exists() ) exit 1, "MultiQC config file not found: ${params.multiqc_config}"
bamtools_filter_pe_config = file(params.bamtools_filter_pe_config)
if( !bamtools_filter_pe_config.exists() ) exit 1, "BamTools filter config file not found: ${params.bamtools_filter_pe_config}"
fastqscreen_config = file(params.fastqscreen_config)

////////////////////////////////////////////////////
/* --          CHECK INPUT FILES               -- */
////////////////////////////////////////////////////

if( params.fasta && file(params.fasta).exists() ){
    Channel
      .fromPath(params.fasta)
      .into { fasta_index_ch;
              fasta_markdup_collectmetrics_ch;
              fasta_merge_replicate_macs2_homer_ch;
              fasta_merge_replicate_macs2_merge_peaks_homer_ch;
              fasta_merge_sample_macs2_homer_ch;
              fasta_merge_sample_macs2_merge_peaks_homer_ch;
              fasta_igv_ch }
} else {
    exit 1, "Genome fasta file not found: ${params.fasta}"
}

if( params.gtf && file(params.gtf).exists() ){
    Channel
      .fromPath(params.gtf)
      .into { gtf_merge_replicate_macs2_homer_ch;
              gtf_merge_replicate_macs2_merge_peaks_homer_ch;
              gtf_merge_sample_macs2_homer_ch;
              gtf_merge_sample_macs2_merge_peaks_homer_ch;
              gtf_igv_ch }
} else {
    exit 1, "GTF file not found: ${params.gtf}"
}

if( params.bwa_index && file("${params.bwa_index}.amb").exists() ) {
    Channel
        .fromPath("${params.bwa_index}*.{amb,ann,bwt,pac,sa}")
        .into { bwa_index_bwa_aln_r1_ch;
                bwa_index_bwa_aln_r2_ch;
                bwa_index_bwa_sampe_ch }
}  else {
    exit 1, "BWA indices not found: ${params.bwa_index}"
}

if( params.design && file(params.design).exists() ){
    Channel
      .fromPath(params.design)
      .splitCsv(header:true, sep:',')
      .map { row -> [ [row.sample,"R"+row.replicate,"L"+row.run].join("_"),
                       row.sample,
                       row.replicate,
                       row.run,
                       file(row.fastq_1),
                       file(row.fastq_2) ] }
      .into { design_raw_fastqc_ch;
              design_raw_fastqscreen_ch;
              design_cutadapt_ch }
} else {
    exit 1, "Design file not found: ${params.design}"
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                       HEADER LOG INFO                               -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

log.info "======================================================"
log.info "BABS ATACSeq Paired-End Pipeline v${params.version}"
log.info "======================================================"
def summary = [:]
summary['Run name']                   = custom_runName ?: workflow.runName
summary['Design file']                = params.design
summary['Genome version']             = params.genome
summary['Genome fasta file']          = params.fasta
summary['GTF annotation file']        = params.gtf
summary['Mitochondrial contig']       = params.mito_name
summary['BWA index']                  = params.bwa_index
summary['Genome masked BED file']     = params.genome_mask
summary['MACS2 genome size']          = params.macs_genome_size
summary['Adapter sequence']           = adapter_seq
summary['Max memory']                 = params.max_memory
summary['Max CPUs']                   = params.max_cpus
summary['Max time']                   = params.max_time
summary['Output directory']           = params.outdir
summary['Output directory path']      = params.outdir_abspath
summary['Working directory']          = workflow.workDir
summary['Current home']               = "$HOME"
summary['Current user']               = "$USER"
summary['Current path']               = "$PWD"
summary['Script directory']           = workflow.projectDir
summary['Config profile']             = workflow.profile
summary['Container']                  = workflow.container
if(workflow.revision) summary['Pipeline release'] = workflow.revision
if(params.project) summary['BABS project'] = params.project
log.info summary.collect { k,v -> "${k.padRight(26)}: $v" }.join("\n")
log.info "======================================================"

// CHECK THAT NEXTFLOW VERSION IS UP TO DATE ENOUGH
try {
    if( ! nextflow.version.matches(">= $params.nf_required_version") ){
        throw GroovyException('Nextflow version too old')
    }
} catch (all) {
    log.error "============================================================\n" +
              "  Nextflow version $params.nf_required_version required! You are running v$workflow.nextflow.version.\n" +
              "  Pipeline execution will continue, but things may break.\n" +
              "  Please run `nextflow self-update` to update Nextflow.\n" +
              "============================================================"
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                     PREPARE ANNOTATION FILES                        -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

// CREATE FAIDX FOR REFERENCE GENOME
// CREATE CHROMOSOME SIZES FILE FOR BEDTOOLS
// CREATE BED FILE WITHOUT MITOCHONDRIAL CONTIG FOR SAMTOOLS FILTERING
process prep_genome {

    tag "$fasta"

    publishDir "${params.outdir}/genome", mode: 'copy'

    input:
    file fasta from fasta_index_ch

    output:
    file "*.fai" into prep_genome_index_ch
    file "*.sizes" into prep_genome_sizes_replicate_bedgraph_ch,
                        prep_genome_sizes_sample_bedgraph_ch,
                        prep_genome_sizes_replicate_bigwig_ch,
                        prep_genome_sizes_sample_bigwig_ch
    file "*.rmMito.bed" into prep_genome_bed_filter_bam_ch

    script:
        prefix = fasta.toString().take(fasta.toString().lastIndexOf('.'))
        """
        samtools faidx ${fasta}
        cut -f 1,2 ${fasta}.fai > ${fasta}.sizes
        awk 'BEGIN{OFS="\t"}{print \$1, '0' , \$2}' ${fasta}.sizes > ${fasta}.bed
        awk '\$1 != "${params.mito_name}"' ${fasta}.bed > ${fasta}.rmMito.bed
        """
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                        FASTQ QC                                     -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

// RUN FASTQC ON RAW FASTQ FILES
process raw_fastqc {

    tag "$sampleid"

    label 'mediumcpu'

    publishDir "${params.outdir}/qc/fastqc/raw", mode: 'copy',
               saveAs: {filename -> filename.endsWith(".zip") ? "zip/$filename" : "$filename"}

    input:
    set val(sampleid), val(sample), val(replicate), val(run), file(fastq_1), file(fastq_2) from design_raw_fastqc_ch

    output:
    set val(sampleid), file("*.{zip,html}") into raw_fastqc_ch

    script:
        """
        ln -s ${fastq_1} ${sampleid}_1.fastq.gz
        ln -s ${fastq_2} ${sampleid}_2.fastq.gz
        fastqc --outdir ./ --threads ${task.cpus} ${sampleid}_1.fastq.gz
        fastqc --outdir ./ --threads ${task.cpus} ${sampleid}_2.fastq.gz
        """
}

// RUN FASTQSCREEN ON RAW FASTQ FILES
process raw_fastqscreen {

    tag "$sampleid"

    label 'mediumcpu'

    publishDir "${params.outdir}/qc/fastq_screen", mode: 'copy',
               saveAs: {filename -> filename.endsWith(".txt") ? "txt/$filename" : "$filename"}

    input:
    set val(sampleid), val(sample), val(replicate), val(run), file(fastq_1), file(fastq_2) from design_raw_fastqscreen_ch

    output:
    set val(sampleid), file("*.{txt,html}") into raw_fastqscreen_ch

    when:
    params.fastqscreen_config

    script:
        """
        ln -s ${fastq_1} ${sampleid}_1.fastq.gz
        ln -s ${fastq_2} ${sampleid}_2.fastq.gz
        fastq_screen --outdir ./ \\
                     --subset 200000 \\
                     --aligner bowtie2 \\
                     --conf ${fastqscreen_config} \\
                     --threads ${task.cpus} \\
                     ${sampleid}_1.fastq.gz
        fastq_screen --outdir ./ \\
                     --subset 200000 \\
                     --aligner bowtie2 \\
                     --conf ${fastqscreen_config} \\
                     --threads ${task.cpus} \\
                     ${sampleid}_2.fastq.gz
        """
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                        ADAPTER TRIMMING                             -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

// RUN CUTADAPT ON RAW FASTQ FILES
process cutadapt {

    tag "$sampleid"

    publishDir "${params.outdir}/qc/cutadapt", mode: 'copy',
                saveAs: {filename ->
                            if (filename.endsWith(".log")) "$filename"
                            else null
                        }

    input:
    set val(sampleid), val(sample), val(replicate), val(run), file(fastq_1), file(fastq_2) from design_cutadapt_ch

    output:
    set val(sampleid), file('*.trim.fastq.gz') into cutadapt_fastqc_ch,
                                                    cutadapt_bwa_aln_r1_ch,
                                                    cutadapt_bwa_aln_r2_ch,
                                                    cutadapt_bwa_sai_to_sam_ch
    set val(sampleid), file('*.log') into cutadapt_log_ch

    script:
        """
        ln -s ${fastq_1} ${sampleid}_1.fastq.gz
        ln -s ${fastq_2} ${sampleid}_2.fastq.gz
        cutadapt -a ${adapter_seq} \\
                 -A ${adapter_seq} \\
                 --minimum-length=25 \\
                 --quality-cutoff=20 \\
                 -o ${sampleid}_1.trim.fastq.gz \\
                 -p ${sampleid}_2.trim.fastq.gz \\
                 ${sampleid}_1.fastq.gz \\
                 ${sampleid}_2.fastq.gz > ${sampleid}.cutadapt.log
        """
}

// RUN FASTQC ON CUTADAPT TRIMMED FASTQ FILES
process trim_fastqc {

    tag "$sampleid"

    label 'mediumcpu'

    publishDir "${params.outdir}/qc/fastqc/trim", mode: 'copy',
               saveAs: {filename -> filename.endsWith(".zip") ? "zip/$filename" : "$filename"}

    input:
    set val(sampleid), file(fastqs) from cutadapt_fastqc_ch

    output:
    set val(sampleid), file("*.{zip,html}") into trim_fastqc_ch

    script:
        """
        fastqc --outdir ./ --threads ${task.cpus} ${fastqs[0]}
        fastqc --outdir ./ --threads ${task.cpus} ${fastqs[1]}
        """
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                        ALIGN                                        -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

// RUN BWA ALN ON CUTADAPT TRIMMED FASTQ FILE FOR READ 1
process bwa_aln_r1 {

    tag "$sampleid"

    label 'highcpu'

    publishDir "${params.outdir}/align", mode: 'copy',
                saveAs: {filename ->
                            if (filename.endsWith(".sysout")) "sysout/$filename"
                            else null
                        }

    input:
    set val(sampleid), file(fastqs) from cutadapt_bwa_aln_r1_ch
    file index from bwa_index_bwa_aln_r1_ch.collect()

    output:
    set val(sampleid), file("*.sai") into bwa_aln_r1_ch
    set val(sampleid), file("*.sysout") into bwa_aln_r1_sysout_ch

    script:
        """
        bwa aln -t ${task.cpus} ${params.bwa_index} ${fastqs[0]} > ${sampleid}_1.sai 2> ${sampleid}_bwa_aln_r1.sysout
        """
}

// RUN BWA ALN ON CUTADAPT TRIMMED FASTQ FILE FOR READ 2
process bwa_aln_r2 {

    tag "$sampleid"

    label 'highcpu'

    publishDir "${params.outdir}/align", mode: 'copy',
                saveAs: {filename ->
                            if (filename.endsWith(".sysout")) "sysout/$filename"
                            else null
                        }

    input:
    set val(sampleid), file(fastqs) from cutadapt_bwa_aln_r2_ch
    file index from bwa_index_bwa_aln_r2_ch.collect()

    output:
    set val(sampleid), file("*.sai") into bwa_aln_r2_ch
    set val(sampleid), file("*.sysout") into bwa_aln_r2_sysout_ch

    script:
        """
        bwa aln -t ${task.cpus} ${params.bwa_index} ${fastqs[1]} > ${sampleid}_2.sai 2> ${sampleid}_bwa_aln_r2.sysout
        """
}

// RUN BWA SAMPE FOR SAI TO SAM CONVERSION
process bwa_sampe {

    tag "$sampleid"

    label "lowcpu"

    publishDir "${params.outdir}/align", mode: 'copy',
                saveAs: {filename ->
                            if (filename.endsWith(".sysout")) "sysout/$filename"
                            else null
                        }

    input:
    set val(sampleid), file(fastqs), file(sai1), file(sai2) from cutadapt_bwa_sai_to_sam_ch.join(bwa_aln_r1_ch, by: [0]).join(bwa_aln_r2_ch, by: [0])
    file index from bwa_index_bwa_sampe_ch.collect()

    output:
    set val(sampleid), file("*.sam") into bwa_sampe_ch
    set val(sampleid), file("*.sysout") into bwa_sampe_sysout_ch

    script:
        rg="\'@RG\\tID:${sampleid}\\tSM:${sampleid.toString().subSequence(0, sampleid.length() - 3)}\\tPL:illumina\\tLB:1\\tPU:1\'"
        """
        bwa sampe -r $rg ${params.bwa_index} ${sai1} ${sai2} ${fastqs} > ${sampleid}.sam 2> ${sampleid}_bwa_sampe.sysout
        """
}

// CONVERT SAM TO COORDINATE SORTED BAM
process bwa_bam {

    tag "$sampleid"

    label 'lowcpu'

    publishDir "${params.outdir}/align", mode: 'copy',
                saveAs: {filename ->
                            if (filename.endsWith(".flagstat")) "flagstat/$filename"
                            else null
                        }

    input:
    set val(sampleid), file(sam) from bwa_sampe_ch

    output:
    set val(sampleid), file("*.sorted.{bam,bam.bai}") into bwa_bam_ch
    set val(sampleid), file("*.flagstat") into bwa_bam_flagstat_ch

    script:
        out_prefix="${sampleid}"
        """
        samtools view -b -h -O BAM -@ ${task.cpus} -o ${out_prefix}.bam ${sam}
        samtools sort -@ ${task.cpus} -o ${out_prefix}.sorted.bam -T ${out_prefix} ${out_prefix}.bam
        samtools index ${out_prefix}.sorted.bam
        samtools flagstat ${out_prefix}.sorted.bam > ${out_prefix}.sorted.bam.flagstat
        """
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                    BAM POST-PROCESSING                              -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

// RUN PICARD MARK DUPLICATES ON COORDINATE SORTED BAM
process markdup {

    tag "$sampleid"

    label 'lowcpu'

    publishDir "${params.outdir}/align", mode: 'copy',
                saveAs: {filename ->
                            if (filename.endsWith(".flagstat")) "flagstat/$filename"
                            else if (filename.endsWith(".idxstats")) "idxstats/$filename"
                            else if (filename.endsWith(".sysout")) "sysout/$filename"
                            else if (filename.endsWith(".metrics.txt")) "picard_metrics/$filename"
                            else null
                        }

    input:
    set val(sampleid), file(bam) from bwa_bam_ch

    output:
    set val(sampleid), file("*.{bam,bam.bai}") into markdup_filter_bam_ch, markdup_collectmetrics_in_ch
    set val(sampleid), file("*.flagstat") into markdup_flagstat_ch
    set val(sampleid), file("*.idxstats") into markdup_idxstats_ch
    set val(sampleid), file("*.metrics.txt") into markdup_metrics_ch
    set val(sampleid), file("*.sysout") into markdup_sysout_ch

    script:
        out_prefix="${sampleid}.mkD"
        """
        java -Xmx${task.memory.toString().split(" ")[0]}g -jar \${EBROOTPICARD}/picard.jar MarkDuplicates \\
             VALIDATION_STRINGENCY=LENIENT \\
             REMOVE_DUPLICATES=false \\
             ASSUME_SORTED=true \\
             TMP_DIR=tmp \\
             INPUT=${bam[0]} \\
             OUTPUT=${out_prefix}.sorted.bam \\
             METRICS_FILE=${out_prefix}.MarkDuplicates.metrics.txt \\
             >> ${out_prefix}.MarkDuplicates.sysout 2>&1
        samtools index ${out_prefix}.sorted.bam
        samtools flagstat ${out_prefix}.sorted.bam > ${out_prefix}.sorted.bam.flagstat
        samtools idxstats ${out_prefix}.sorted.bam > ${out_prefix}.sorted.bam.idxstats
        """
}

// RUN PICARD COLLECTMULTIPLE METRICS ON COORDINATE SORTED BAM
process markdup_collectmetrics {

    tag "$sampleid"

    label 'lowcpu'

    publishDir "${params.outdir}/align", mode: 'copy',
                saveAs: {filename ->
                            if (filename.endsWith("_metrics")) "picard_metrics/$filename"
                            else if (filename.endsWith(".pdf")) "picard_metrics/pdf/$filename"
                            else if (filename.endsWith(".sysout")) "sysout/$filename"
                        }

    input:
    set val(sampleid), file(bam) from markdup_collectmetrics_in_ch
    file fasta from fasta_markdup_collectmetrics_ch.collect()

    output:
    set val(sampleid), file("*{_metrics,.pdf}") into markdup_collectmetrics_ch
    set val(sampleid), file("*.sysout") into markdup_collectmetrics_sysout_ch

    script:
        out_prefix="${sampleid}.mkD"
        """
        java -Xmx${task.memory.toString().split(" ")[0]}g -jar \${EBROOTPICARD}/picard.jar CollectMultipleMetrics \\
             VALIDATION_STRINGENCY=LENIENT \\
             TMP_DIR=tmp \\
             INPUT=${bam[0]} \\
             OUTPUT=${out_prefix}.CollectMultipleMetrics \\
             REFERENCE_SEQUENCE=${fasta} \\
             >> ${out_prefix}.CollectMultipleMetrics.sysout 2>&1
        """
}

// FILTER BAM FILE TO KEEP (1) UNIQUELY MAPPED, (2) PRIMARY ALIGNMENT, (3) PROPERLY-PAIRED, (4) NON-MITOCHONDRIAL
//                         (5) NON-SOFTCLIPPED (bamtools), (6) INSERT SIZE < 2KB (bamtools), (7) MISMATCH <= 3 (bamtools)
process filter_bam {

    tag "$sampleid"

    label 'lowcpu'

    publishDir "${params.outdir}/align", mode: 'copy',
                saveAs: {filename ->
                            if (filename.endsWith(".flagstat")) "flagstat/$filename"
                            else null
                        }

    input:
    set val(sampleid), file(bam) from markdup_filter_bam_ch
    file mito_bed from prep_genome_bed_filter_bam_ch.collect()

    output:
    set val(sampleid), file("*.mkD.clN.bam") into filter_bam_ch
    set val(sampleid), file("*.sorted.{bam,bam.bai}") into filter_bam_sort_ch
    set val(sampleid), file("*.flagstat") into filter_bam_flagstat_ch

    script:
        // 0x0001 = read paired
        // 0x0002 = read mapped in proper pair
        // 0x0004 = read unmapped
        // 0x0008 = mate unmapped
        // 0x0100 = not primary alignment
        // 0x0400 = read is PCR or optical duplicate
        out_prefix="${sampleid}.mkD.clN"
        """
        samtools view \\
                 -f 0x001 \\
                 -f 0x002 \\
                 -F 0x004 \\
                 -F 0x0008 \\
                 -F 0x0100 \\
                 -q 1 \\
                 -L ${mito_bed} \\
                 -b ${bam[0]} \\
                 | bamtools filter \\
                            -out ${out_prefix}.sorted.bam \\
                            -script ${bamtools_filter_pe_config}
        samtools index ${out_prefix}.sorted.bam
        samtools flagstat ${out_prefix}.sorted.bam > ${out_prefix}.sorted.bam.flagstat
        samtools sort -n -@ ${task.cpus} -o ${out_prefix}.bam -T ${out_prefix} ${out_prefix}.sorted.bam
        """
}

// FILTER BAM FILE TO REMOVE READS WHERE ONE MATE HAS BEEN PHYSICALLY FILTERED OUT OF BAM FILE FROM PREVIOUS STEP.
process rm_orphan {

    tag "$sampleid"

    input:
    set val(sampleid), file(bam) from filter_bam_ch

    output:
    set val(sampleid), file("*.bam") into rm_orphan_ch

    script:
        out_prefix="${sampleid}.mkD.flT"
        """
        python $baseDir/bin/rm_orphan_from_bampe.py ${bam} ${out_prefix}.bam --excl_diff_chrom
        """
}

// SORT ORPHAN BAM FILE BY COORDINATE
process  rm_orphan_sort_bam {

    tag "$sampleid"

    label 'lowcpu'

    publishDir "${params.outdir}/align", mode: 'copy',
                saveAs: {filename ->
                            if (filename.endsWith(".flagstat")) "flagstat/$filename"
                            else if (filename.endsWith(".idxstats")) "idxstats/$filename"
                            else if (filename.endsWith(".bam")) "$filename"
                            else if (filename.endsWith(".bai")) "$filename"
                        }

    input:
    set val(sampleid), file(filter_bam), file(orphan_bam) from filter_bam_sort_ch.join(rm_orphan_ch, by: [0])

    output:
    set val(sampleid), file("*.{bam,bam.bai,flagstat}") into rm_orphan_sort_bam_replicate_ch,
                                                             rm_orphan_sort_bam_sample_ch,
                                                             rm_orphan_sort_bam_replicate_rmdup_ch,
                                                             rm_orphan_sort_bam_sample_rmdup_ch
    set val(sampleid), file("*.idxstats") into rm_orphan_sort_bam_idxstats_ch

    script:
        out_prefix="${sampleid}.mkD.flT"
        """
        samtools sort -@ ${task.cpus} -o ${out_prefix}.sorted.bam -T ${out_prefix} ${out_prefix}.bam
        samtools index ${out_prefix}.sorted.bam
        samtools flagstat ${out_prefix}.sorted.bam > ${out_prefix}.sorted.bam.flagstat
        samtools idxstats ${out_prefix}.sorted.bam > ${out_prefix}.sorted.bam.idxstats
        """
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                    MERGE REPLICATE BAM                              -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

// CREATE CHANNEL TO MERGE AT REPLICATE LEVEL
rm_orphan_sort_bam_replicate_ch.map { it -> [ it[0].toString().subSequence(0, it[0].length() - 3), it[1] ] }
                               .groupTuple(by: [0])
                               .map { it ->  [ it[0], it[1].flatten() ] }
                               .into { rm_orphan_sort_bam_replicate_merge_ch;
                                       rm_orphan_sort_bam_replicate_markdup_ch;
                                       rm_orphan_sort_bam_replicate_rmdup_ch }

// MERGE FILTERED BAM FILES AT REPLICATE LEVEL USING PICARD MERGESAMFILES
process merge_replicate {

    tag "$sampleid"

    label 'lowcpu'

    publishDir "${params.outdir}/align/mergeReplicate", mode: 'copy',
                saveAs: {filename ->
                            if (filename.endsWith(".flagstat")) "flagstat/$filename"
                            else if (filename.endsWith(".sysout")) "sysout/$filename"
                            else null
                        }

    input:
    set val(sampleid), file(bams) from rm_orphan_sort_bam_replicate_merge_ch

    output:
    set val(sampleid), file("*.{bam,bam.bai}") into merge_replicate_ch
    set val(sampleid), file("*.flagstat") into merge_replicate_flagstat_ch
    set val(sampleid), file("*.sysout") into merge_replicate_sysout_ch

    script:
        out_prefix="${sampleid}.mRp"
        bam_files = bams.findAll { it.toString().endsWith('.bam') }.sort()
        flagstat_files = bams.findAll { it.toString().endsWith('.flagstat') }.sort()
        if (bam_files.size() > 1) {
            """
            java -Xmx${task.memory.toString().split(" ")[0]}g -jar \${EBROOTPICARD}/picard.jar MergeSamFiles \\
                 VALIDATION_STRINGENCY=LENIENT \\
                 SORT_ORDER=coordinate \\
                 TMP_DIR=tmp \\
                 ${'INPUT='+bam_files.join(' INPUT=')} \\
                 OUTPUT=${out_prefix}.sorted.bam \\
                 >> ${out_prefix}.MergeSamFiles.sysout 2>&1
            samtools index ${out_prefix}.sorted.bam
            samtools flagstat ${out_prefix}.sorted.bam > ${out_prefix}.sorted.bam.flagstat
            """
        } else {
            """
            touch ${out_prefix}.sorted.bam
            touch ${out_prefix}.sorted.bam.bai
            touch ${out_prefix}.MergeSamFiles.sysout
            cp ${flagstat_files[0]} ${out_prefix}.sorted.bam.flagstat
            """
        }
}

// RUN PICARD MARK DUPLICATES
process merge_replicate_markdup {

    tag "$sampleid"

    label 'lowcpu'

    publishDir "${params.outdir}/align/mergeReplicate", mode: 'copy',
                saveAs: {filename ->
                            if (filename.endsWith(".flagstat")) "flagstat/$filename"
                            else if (filename.endsWith(".sysout")) "sysout/$filename"
                            else if (filename.endsWith(".metrics.txt")) "picard_metrics/$filename"
                            else null
                        }

    input:
    set val(sampleid), file(orphan_bams), file(merged_bam) from rm_orphan_sort_bam_replicate_markdup_ch.join(merge_replicate_ch, by: [0])

    output:
    set val(sampleid), file("*.{bam,bam.bai}") into merge_replicate_markdup_ch
    set val(sampleid), file("*.flagstat") into merge_replicate_markdup_flagstat_ch
    set val(sampleid), file("*.metrics.txt") into merge_replicate_markdup_metrics_ch
    set val(sampleid), file("*.sysout") into merge_replicate_markdup_sysout_ch

    script:
        out_prefix="${sampleid}.mRp.mkD"
        bam_files = orphan_bams.findAll { it.toString().endsWith('.bam') }.sort()
        flagstat_files = orphan_bams.findAll { it.toString().endsWith('.flagstat') }.sort()
        if (bam_files.size() > 1) {
            """
            java -Xmx${task.memory.toString().split(" ")[0]}g -jar \${EBROOTPICARD}/picard.jar MarkDuplicates \\
                 VALIDATION_STRINGENCY=LENIENT \\
                 REMOVE_DUPLICATES=false \\
                 ASSUME_SORTED=true \\
                 TMP_DIR=tmp \\
                 INPUT=${merged_bam[0]} \\
                 OUTPUT=${out_prefix}.sorted.bam \\
                 METRICS_FILE=${out_prefix}.MarkDuplicates.metrics.txt \\
                 >> ${out_prefix}.MarkDuplicates.sysout 2>&1
            samtools index ${out_prefix}.sorted.bam
            samtools flagstat ${out_prefix}.sorted.bam > ${out_prefix}.sorted.bam.flagstat
            """
        } else {
            """
            touch ${out_prefix}.sorted.bam
            touch ${out_prefix}.sorted.bam.bai
            touch ${out_prefix}.MarkDuplicates.sysout
            touch ${out_prefix}.MarkDuplicates.metrics.txt
            cp ${flagstat_files[0]} ${out_prefix}.sorted.bam.flagstat
            """
        }
}

// REMOVE DUPLICATES USING SAMTOOLS
process merge_replicate_rmdup {

    tag "$sampleid"

    publishDir "${params.outdir}/align/mergeReplicate", mode: 'copy',
                saveAs: {filename ->
                            if (filename.endsWith(".flagstat")) "flagstat/$filename"
                            else if (filename.endsWith(".bam")) "$filename"
                            else if (filename.endsWith(".bai")) "$filename"
                        }

    input:
    set val(sampleid), file(orphan_bams), file(markdup_bam) from rm_orphan_sort_bam_replicate_rmdup_ch.join(merge_replicate_markdup_ch, by: [0])

    output:
    set val(sampleid), file("*.{bam,bam.bai}") into merge_replicate_rmdup_name_bam_ch,
                                                    merge_replicate_rmdup_bedgraph_ch,
                                                    merge_replicate_rmdup_macs2_ch,
                                                    merge_replicate_rmdup_macs2_frip_ch
    set val(sampleid), file("*.flagstat") into merge_replicate_rmdup_flagstat_ch,
                                               merge_replicate_rmdup_flagstat_bedgraph_ch,
                                               merge_replicate_rmdup_flagstat_macs2_frip_ch

    script:
        out_prefix="${sampleid}.mRp.rmD"
        bam_files = orphan_bams.findAll { it.toString().endsWith('.bam') }.sort()
        if (bam_files.size() > 1) {
            """
            samtools view -bF 0x400 ${markdup_bam[0]} > ${out_prefix}.sorted.bam
            samtools index ${out_prefix}.sorted.bam
            samtools flagstat ${out_prefix}.sorted.bam > ${out_prefix}.sorted.bam.flagstat
            """
        } else {
            """
            samtools view -bF 0x400 ${bam_files[0]} > ${out_prefix}.sorted.bam
            samtools index ${out_prefix}.sorted.bam
            samtools flagstat ${out_prefix}.sorted.bam > ${out_prefix}.sorted.bam.flagstat
            """
        }
}

// SORT BAM FILE BY NAME
process merge_replicate_name_bam {

    tag "$sampleid"

    label 'lowcpu'

    input:
    set val(sampleid), file(bam) from merge_replicate_rmdup_name_bam_ch

    output:
    set val(sampleid), file("*.bam") into merge_replicate_name_bam_merge_replicate_featurecounts_ch,
                                         merge_replicate_name_bam_merge_sample_featurecounts_ch

    script:
        out_prefix="${sampleid}.mRp.rmD"
        """
        samtools sort -n -@ ${task.cpus} -o ${out_prefix}.bam -T ${out_prefix} ${out_prefix}.sorted.bam
        """
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                 MERGE REPLICATE BAM POST-ANALYSIS                   -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////
/* --          COVERAGE TRACKS                 -- */
////////////////////////////////////////////////////

// CREATE NORMALISED BEDGRAPH FILE USING BEDTOOLS GENOMECOVERAGEBED
// CALCULATE SCALE-FACTOR FROM FLAGSTAT FILE
process merge_replicate_bedgraph {

    tag "$sampleid"

    publishDir "${params.outdir}/align/mergeReplicate/bigwig", mode: 'copy',
                saveAs: {filename ->
                            if (filename.endsWith(".txt")) "scale_factor/$filename"
                            else null
                        }

    input:
    set val(sampleid), file(bam), file(flagstat) from merge_replicate_rmdup_bedgraph_ch.join(merge_replicate_rmdup_flagstat_bedgraph_ch, by: [0])
    file chrom_sizes from prep_genome_sizes_replicate_bedgraph_ch.collect()

    output:
    set val(sampleid), file("*.bg") into merge_replicate_bedgraph_ch
    set val(sampleid), file("*.txt") into merge_replicate_bedgraph_scale_ch

    script:
        out_prefix="${sampleid}.mRp.rmD"
        """
        SCALE_FACTOR=\$(grep 'read1' ${flagstat} | awk '{print 1000000/\$1}')
        echo \$SCALE_FACTOR > ${out_prefix}.scale_factor.txt
        genomeCoverageBed -ibam ${bam[0]} -bg -trackline -scale \$SCALE_FACTOR -pc -g ${chrom_sizes} >  ${out_prefix}.bg
        """
}

// CONVERT BEDGRAPH TO BIGWIG USING KENTOOLS
process merge_replicate_bigwig {

    tag "$sampleid"

    publishDir "${params.outdir}/align/mergeReplicate/bigwig", mode: 'copy'

    input:
    set val(sampleid), file(bedgraph) from merge_replicate_bedgraph_ch
    file chrom_sizes from prep_genome_sizes_replicate_bigwig_ch.collect()

    output:
    set val(sampleid), file("*.bigWig") into merge_replicate_bigwig_ch

    script:
        out_prefix="${sampleid}.mRp.rmD"
        """

        wigToBigWig -clip ${bedgraph}  ${chrom_sizes} ${out_prefix}.bigWig
        """
}

////////////////////////////////////////////////////
/* --          MACS2 BROAD                     -- */
////////////////////////////////////////////////////

// CALL BROAD PEAKS USING MACS2
process merge_replicate_macs2 {

    tag "$sampleid"

    publishDir "${params.outdir}/align/mergeReplicate/macs2", mode: 'copy'

    input:
    set val(sampleid), file(bam) from merge_replicate_rmdup_macs2_ch

    output:
    set val(sampleid), file("*/*.broadPeak") into merge_replicate_macs2_homer_in_ch,
                                                  merge_replicate_macs2_frip_in_ch,
                                                  merge_replicate_macs2_merge_peaks_in_ch
    set val(sampleid), file("*/*.{gappedPeak,xls}") into merge_replicate_macs2_output_ch
    set val(sampleid), file("*/*.sysout") into merge_replicate_macs2_sysout_ch

    script:
        """
        mkdir -p ${sampleid}
        macs2 callpeak \\
              --verbose=2 \\
              -t ${bam[0]} \\
              --gsize=${params.macs_genome_size} \\
              --outdir=${sampleid} \\
              --name=${sampleid} \\
              --format BAMPE \\
              --keep-dup all \\
              --nomodel \\
              --broad \\
              >> ${sampleid}/${sampleid}.macs2.sysout 2>&1
        """
}

// ANNOTATE PEAKS USING HOMER ANNOTATEPEAKS
process merge_replicate_macs2_homer {

    tag "$sampleid"

    publishDir "${params.outdir}/align/mergeReplicate/macs2", mode: 'copy'

    input:
    set val(sampleid), file(peak) from merge_replicate_macs2_homer_in_ch
    file fasta from fasta_merge_replicate_macs2_homer_ch.collect()
    file gtf from gtf_merge_replicate_macs2_homer_ch.collect()

    output:
    set val(sampleid), file("*/*.annotatePeaks.txt") into merge_replicate_macs2_homer_ch
    set val(sampleid), file("*/*.annotatePeaks.sysout") into merge_replicate_macs2_homer_sysout_ch

    script:
        """
        mkdir -p ${sampleid}
        annotatePeaks.pl ${peak} \\
                         ${fasta} \\
                         -gid \\
                         -gtf ${gtf} \\
                         > ${sampleid}/${sampleid}_peaks.homer.annotatePeaks.txt \\
                         2> ${sampleid}/${sampleid}_peaks.homer.annotatePeaks.sysout
        """
}

// CALCULATE FRIP SCORE USING BEDTOOLS
process merge_replicate_macs2_frip {

    tag "$sampleid"

    publishDir "${params.outdir}/align/mergeReplicate/macs2", mode: 'copy'

    input:
    set val(sampleid), file(peak), file(bam), file(flagstat) from merge_replicate_macs2_frip_in_ch.join(merge_replicate_rmdup_macs2_frip_ch, by: [0]).join(merge_replicate_rmdup_flagstat_macs2_frip_ch, by: [0])

    output:
    set val(sampleid), file("*/*.frip.txt") into merge_replicate_macs2_frip_ch

    script:
        """
        mkdir -p ${sampleid}
        READS_IN_PEAKS=\$(intersectBed -a ${bam[0]} -b ${peak} -bed -c -f 0.20 | awk -F '\t' '{sum += \$NF} END {print sum}')
        grep 'mapped (' ${flagstat} | awk -v a="\$READS_IN_PEAKS" '{print a/\$1}' > ${sampleid}/${sampleid}_peaks.frip.txt
        """
}

// GET LIST OF FRIP FILES ACROSS ALL SAMPLES FOR QC PLOT
merge_replicate_macs2_frip_ch.map { it -> [ "macs2_frip", [it[1]] ] }
                             .groupTuple(by: [0])
                             .map { it ->  [ it[0],
                                             it[1].flatten().sort().collect{ it.getName().replace("_peaks.frip.txt","") },
                                             it[1].flatten().sort() ] }
                             .set { merge_replicate_macs2_frip_collate_ch }

// GENERATE PLOTS FOR FRIP SCORES ACROSS ALL SAMPLES
process merge_replicate_macs2_fripqc {

   tag "$name"

   publishDir "${params.outdir}/align/mergeReplicate/macs2/qc", mode: 'copy'

   input:
   set val(name), val(sampleids), file(frips) from merge_replicate_macs2_frip_collate_ch

   output:
   file ("*.pdf") into merge_replicate_macs2_fripqc_ch

   script:
       """
       Rscript $baseDir/bin/plot_frip.R -i ${frips.join(',')} -s ${sampleids.join(',')} -o ./ -p ${name}
       """
}

// GET LIST OF HOMER ANNOTATED FILES ACROSS ALL SAMPLES FOR QC PLOTS
merge_replicate_macs2_homer_ch.map { it -> [ "macs2_homer_annotation", [it[1]] ] }
                              .groupTuple(by: [0])
                              .map { it ->  [ it[0],
                                              it[1].flatten().sort().collect{ it.getName().replace("_peaks.homer.annotatePeaks.txt","") },
                                              it[1].flatten().sort() ] }
                              .set { merge_replicate_macs2_homer_collate_ch }

// GENERATE PLOTS FOR VARIOUS ASPECTS OF HOMER ANNOTATION ACROSS ALL SAMPLES
process merge_replicate_macs2_homerqc {

   tag "$name"

   publishDir "${params.outdir}/align/mergeReplicate/macs2/qc", mode: 'copy'

   input:
   set val(name), val(sampleids), file(homers) from merge_replicate_macs2_homer_collate_ch

   output:
   file ("*.pdf") into merge_replicate_macs2_homerqc_ch

   script:
       """
       Rscript $baseDir/bin/homer_plot_annotation.R -i ${homers.join(',')} -s ${sampleids.join(',')} -o ./ -p ${name}
       """
}

// GET LIST OF PEAK ACROSS ALL SAMPLES FOR MERGING
merge_replicate_macs2_merge_peaks_in_ch.map { it -> [ "merged_peaks", [it[1]] ] }
                                       .groupTuple(by: [0])
                                       .map { it ->  [ it[0],
                                                       it[1].flatten().sort().collect{ it.getName().replace("_peaks.broadPeak","") },
                                                       it[1].flatten().sort() ] }
                                       .into { merge_replicate_macs2_merge_peaks_in_ch;
                                               merge_replicate_macs2_peakqc_in_ch }

// GENERATE PLOTS FOR VARIOUS ASPECTS OF PEAKS ACROSS ALL SAMPLES
process merge_replicate_macs2_peakqc {

   tag "$name"

   publishDir "${params.outdir}/align/mergeReplicate/macs2/qc", mode: 'copy'

   input:
   set val(name), val(sampleids), file(peaks) from merge_replicate_macs2_peakqc_in_ch

   output:
   file ("*.{txt,pdf}") into merge_replicate_macs2_peakqc_ch

   script:
       """
       Rscript $baseDir/bin/macs2_peakqc.R -i ${peaks.join(',')} -s ${sampleids.join(',')} -o ./ -p macs2_peakqc
       """
}

// MERGE PEAKS ACROSS ALL SAMPLES AND CREATE A BOOLEAN OUTPUT FILE TO AID FILTERING
process merge_replicate_macs2_merge_peaks {

   tag "$name"

   publishDir "${params.outdir}/align/mergeReplicate/macs2/merged_peaks", mode: 'copy'

   input:
   set val(name), val(sampleids), file(peaks) from merge_replicate_macs2_merge_peaks_in_ch

   output:
   set val(name), file("*.boolean.txt") into merge_replicate_macs2_merge_peaks_bool_ch
   set val(name), file("*.bed") into merge_replicate_macs2_merge_peaks_bed_ch
   set val(name), file("*.saf") into merge_replicate_macs2_merge_peaks_saf_ch
   set val(name), file("*.log") into merge_replicate_macs2_merge_peaks_log_ch

   script:
       """
       sort -k1,1 -k2,2n ${peaks.join(' ')} \\
            | mergeBed -c 2,3,4,5,6,7,8,9 -o collapse,collapse,collapse,collapse,collapse,collapse,collapse,collapse > ${name}.txt
       python $baseDir/bin/expand_merged_macs.py ${name}.txt ${sampleids.join(',')} ${name}.boolean.txt --min_samples 1
       awk -v FS='\t' -v OFS='\t' 'FNR > 1 { print \$1, \$2, \$3, \$4, "0", "+" }' ${name}.boolean.txt > ${name}.bed
       echo -e "GeneID\tChr\tStart\tEnd\tStrand" > ${name}.saf
       awk -v FS='\t' -v OFS='\t' 'FNR > 1 { print \$4, \$1, \$2, \$3,  "+" }' ${name}.boolean.txt >> ${name}.saf
       """
}

// ANNOTATE PEAKS USING HOMER ANNOTATE AND ADD ANNOTATION TO BOOLEAN OUTPUT FILE
process merge_replicate_macs2_merge_peaks_homer {

    tag "${name}"

    publishDir "${params.outdir}/align/mergeReplicate/macs2/merged_peaks", mode: 'copy'

    input:
    set val(name), file(bed) from merge_replicate_macs2_merge_peaks_bed_ch
    set val(name), file(bool) from merge_replicate_macs2_merge_peaks_bool_ch
    file fasta from fasta_merge_replicate_macs2_merge_peaks_homer_ch.collect()
    file gtf from gtf_merge_replicate_macs2_merge_peaks_homer_ch.collect()

    output:
    set val(name), file("*.annotatePeaks.txt") into merge_replicate_macs2_merge_peaks_homer_ch
    set val(name), file("*.annotatePeaks.sysout") into merge_replicate_macs2_merge_peaks_homer_sysout_ch

    script:
        """
        annotatePeaks.pl ${bed} \\
                         ${fasta} \\
                         -gid \\
                         -gtf \\
                         ${gtf} \\
                         > ${name}.homer.annotatePeaks.txt \\
                         2> ${name}.homer.annotatePeaks.sysout
        cut -f2- ${name}.homer.annotatePeaks.txt | awk 'NR==1; NR > 1 {print \$0 | "sort -k1,1 -k2,2n"}' | cut -f6- > tmp.txt
        paste ${bool} tmp.txt > ${bool.toString().replace(".txt","")}.homer.annotatePeaks.txt
        """
}

// GET LIST OF REPLICATE LEVEL BAM FILES
merge_replicate_name_bam_merge_replicate_featurecounts_ch.map { it -> [ "merged_bam", [it[1]] ] }
                                                         .groupTuple(by: [0])
                                                         .map { it ->  [ it[1].flatten().sort().collect{ it.getName().replace(".mRp.rmD.bam","") },
                                                                         it[1].flatten().sort() ] }
                                                         .set { merge_replicate_name_bam_merge_replicate_featurecounts_ch }

process merge_replicate_macs2_merge_peaks_featurecounts {

    tag "${name}"

    label 'highcpu'

    publishDir "${params.outdir}/align/mergeReplicate/macs2/merged_peaks", mode: 'copy'

    input:
    set val (sampleids), file(bams) from merge_replicate_name_bam_merge_replicate_featurecounts_ch
    set val(name), file(saf) from merge_replicate_macs2_merge_peaks_saf_ch.collect()

    output:
    set val(name), file("*.featureCounts.txt") into merge_replicate_macs2_merge_peaks_featurecounts_ch
    set val(name), file("*.featureCounts.sysout") into merge_replicate_macs2_merge_peaks_featurecounts_sysout_ch

    script:
        """
        featureCounts -F SAF \\
                      -O \\
                      -T ${task.cpus} \\
                      --fracOverlap 0.5 \\
                      --primary \\
                      --ignoreDup \\
                      -p \\
                      -B \\
                      -C \\
                      --donotsort \\
                      -a ${saf} \\
                      -o ${name}.replicate.featureCounts.txt \\
                      ${bams.join(' ')} \\
                      >> ${name}.replicate.featureCounts.sysout 2>&1
        """
}

process merge_replicate_macs2_merge_peaks_differential {

    tag "${name}"

    publishDir "${params.outdir}/align/mergeReplicate/macs2/merged_peaks", mode: 'copy'

    input:
    set val(name), file(counts) from merge_replicate_macs2_merge_peaks_featurecounts_ch

    output:
    set val(name), file("deseq2/*") into merge_replicate_macs2_merge_peaks_differential_ch
    file "*.txt" into merge_replicate_macs2_merge_peaks_differential_complete_ch

    script:
        """
        Rscript $baseDir/bin/featurecounts_deseq2.R -i ${counts} -b '.mRp.rmD.bam' -o ./deseq2 -p ${name}
        touch merge_replicate_macs2_merge_peaks_differential.complete.txt
        """
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                    MERGE SAMPLE BAM                                 -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

// CREATE CHANNEL TO MERGE AT SAMPLE LEVEL
rm_orphan_sort_bam_sample_ch.map { it -> [ it[0].toString().subSequence(0, it[0].length() - 6), it[1] ] }
                            .groupTuple(by: [0])
                            .map { it ->  [ it[0],
                                            it[1].flatten() ] }
                            .into { rm_orphan_sort_bam_sample_merge_ch;
                                    rm_orphan_sort_bam_sample_markdup_ch;
                                    rm_orphan_sort_bam_sample_rmdup_ch }

// MERGE FILTERED BAM FILES AT SAMPLE LEVEL USING PICARD MERGESAMFILES
process merge_sample {

    tag "$sampleid"

    label 'lowcpu'

    publishDir "${params.outdir}/align/mergeSample", mode: 'copy',
                saveAs: {filename ->
                            if (filename.endsWith(".flagstat")) "flagstat/$filename"
                            else if (filename.endsWith(".sysout")) "sysout/$filename"
                            else null
                        }

    input:
    set val(sampleid), file(bams) from rm_orphan_sort_bam_sample_merge_ch

    output:
    set val(sampleid), file("*.{bam,bam.bai}") into merge_sample_ch
    set val(sampleid), file("*.flagstat") into merge_sample_flagstat_ch
    set val(sampleid), file("*.sysout") into merge_sample_sysout_ch

    script:
        out_prefix="${sampleid}.mSm"
        bam_files = bams.findAll { it.toString().endsWith('.bam') }.sort()
        flagstat_files = bams.findAll { it.toString().endsWith('.flagstat') }.sort()
        if (bam_files.size() > 1) {
            """
            java -Xmx${task.memory.toString().split(" ")[0]}g -jar \${EBROOTPICARD}/picard.jar MergeSamFiles \\
                 VALIDATION_STRINGENCY=LENIENT \\
                 SORT_ORDER=coordinate \\
                 TMP_DIR=tmp \\
                 ${'INPUT='+bam_files.join(' INPUT=')} \\
                 OUTPUT=${out_prefix}.sorted.bam \\
                 >> ${out_prefix}.MergeSamFiles.sysout 2>&1
            samtools index ${out_prefix}.sorted.bam
            samtools flagstat ${out_prefix}.sorted.bam > ${out_prefix}.sorted.bam.flagstat
            """
        } else {
            """
            touch ${out_prefix}.sorted.bam
            touch ${out_prefix}.sorted.bam.bai
            touch ${out_prefix}.MergeSamFiles.sysout
            cp ${flagstat_files[0]} ${out_prefix}.sorted.bam.flagstat
            """
        }
}

// RUN PICARD MARK DUPLICATES
process merge_sample_markdup {

    tag "$sampleid"

    label 'lowcpu'

    publishDir "${params.outdir}/align/mergeSample", mode: 'copy',
                saveAs: {filename ->
                            if (filename.endsWith(".flagstat")) "flagstat/$filename"
                            else if (filename.endsWith(".sysout")) "sysout/$filename"
                            else if (filename.endsWith(".metrics.txt")) "picard_metrics/$filename"
                            else null
                        }

    input:
    set val(sampleid), file(orphan_bams), file(merged_bam) from rm_orphan_sort_bam_sample_markdup_ch.join(merge_sample_ch, by: [0])

    output:
    set val(sampleid), file("*.{bam,bam.bai}") into merge_sample_markdup_ch
    set val(sampleid), file("*.flagstat") into merge_sample_markdup_flagstat_ch
    set val(sampleid), file("*.metrics.txt") into merge_sample_markdup_metrics_ch
    set val(sampleid), file("*.sysout") into merge_sample_markdup_sysout_ch

    script:
        out_prefix="${sampleid}.mSm.mkD"
        bam_files = orphan_bams.findAll { it.toString().endsWith('.bam') }.sort()
        flagstat_files = orphan_bams.findAll { it.toString().endsWith('.flagstat') }.sort()
        if (bam_files.size() > 1) {
            """
            java -Xmx${task.memory.toString().split(" ")[0]}g -jar \${EBROOTPICARD}/picard.jar MarkDuplicates \\
                 VALIDATION_STRINGENCY=LENIENT \\
                 REMOVE_DUPLICATES=false \\
                 ASSUME_SORTED=true \\
                 TMP_DIR=tmp \\
                 INPUT=${merged_bam[0]} \\
                 OUTPUT=${out_prefix}.sorted.bam \\
                 METRICS_FILE=${out_prefix}.MarkDuplicates.metrics.txt \\
                 >> ${out_prefix}.MarkDuplicates.sysout 2>&1
            samtools index ${out_prefix}.sorted.bam
            samtools flagstat ${out_prefix}.sorted.bam > ${out_prefix}.sorted.bam.flagstat
            """
        } else {
            """
            touch ${out_prefix}.sorted.bam
            touch ${out_prefix}.sorted.bam.bai
            touch ${out_prefix}.MarkDuplicates.sysout
            touch ${out_prefix}.MarkDuplicates.metrics.txt
            cp ${flagstat_files[0]} ${out_prefix}.sorted.bam.flagstat
            """
        }
}

// REMOVE DUPLICATES USING SAMTOOLS
process merge_sample_rmdup {

    tag "$sampleid"

    publishDir "${params.outdir}/align/mergeSample", mode: 'copy',
                saveAs: {filename ->
                            if (filename.endsWith(".flagstat")) "flagstat/$filename"
                            else if (filename.endsWith(".bam")) "$filename"
                            else if (filename.endsWith(".bai")) "$filename"
                        }

    input:
    set val(sampleid), file(orphan_bams), file(markdup_bam) from rm_orphan_sort_bam_sample_rmdup_ch.join(merge_sample_markdup_ch, by: [0])

    output:
    set val(sampleid), file("*.{bam,bam.bai}") into merge_sample_rmdup_bedgraph_ch,
                                                    merge_sample_rmdup_macs2_ch,
                                                    merge_sample_rmdup_macs2_frip_ch
    set val(sampleid), file("*.flagstat") into merge_sample_rmdup_flagstat_ch,
                                               merge_sample_rmdup_flagstat_bedgraph_ch,
                                               merge_sample_rmdup_flagstat_macs2_frip_ch

    script:
        out_prefix="${sampleid}.mSm.rmD"
        bam_files = orphan_bams.findAll { it.toString().endsWith('.bam') }.sort()
        if (bam_files.size() > 1) {
            """
            samtools view -bF 0x400 ${markdup_bam[0]} > ${out_prefix}.sorted.bam
            samtools index ${out_prefix}.sorted.bam
            samtools flagstat ${out_prefix}.sorted.bam > ${out_prefix}.sorted.bam.flagstat
            """
        } else {
            """
            samtools view -bF 0x400 ${bam_files[0]} > ${out_prefix}.sorted.bam
            samtools index ${out_prefix}.sorted.bam
            samtools flagstat ${out_prefix}.sorted.bam > ${out_prefix}.sorted.bam.flagstat
            """
        }
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                 MERGE SAMPLE BAM POST-ANALYSIS                      -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////
/* --          COVERAGE TRACKS                 -- */
////////////////////////////////////////////////////

// CREATE NORMALISED BEDGRAPH FILE USING BEDTOOLS GENOMECOVERAGEBED
// CALCULATE SCALE-FACTOR FROM FLAGSTAT FILE
process merge_sample_bedgraph {

    tag "$sampleid"

    publishDir "${params.outdir}/align/mergeSample/bigwig", mode: 'copy',
                saveAs: {filename ->
                            if (filename.endsWith(".txt")) "scale_factor/$filename"
                            else null
                        }

    input:
    set val(sampleid), file(bam), file(flagstat) from merge_sample_rmdup_bedgraph_ch.join(merge_sample_rmdup_flagstat_bedgraph_ch, by: [0])
    file chrom_sizes from prep_genome_sizes_sample_bedgraph_ch.collect()

    output:
    set val(sampleid), file("*.bg") into merge_sample_bedgraph_ch
    set val(sampleid), file("*.txt") into merge_sample_bedgraph_scale_ch

    script:
        out_prefix="${sampleid}.mSm.rmD"
        """
        SCALE_FACTOR=\$(grep 'read1' ${flagstat} | awk '{print 1000000/\$1}')
        echo \$SCALE_FACTOR > ${out_prefix}.scale_factor.txt
        genomeCoverageBed -ibam ${bam[0]} -bg -trackline -scale \$SCALE_FACTOR -pc -g ${chrom_sizes} >  ${out_prefix}.bg
        """
}

// CONVERT BEDGRAPH TO BIGWIG USING KENTOOLS
process merge_sample_bigwig {

    tag "$sampleid"

    publishDir "${params.outdir}/align/mergeSample/bigwig", mode: 'copy'

    input:
    set val(sampleid), file(bedgraph) from merge_sample_bedgraph_ch
    file chrom_sizes from prep_genome_sizes_sample_bigwig_ch.collect()

    output:
    set val(sampleid), file("*.bigWig") into merge_sample_bigwig_ch

    script:
        out_prefix="${sampleid}.mSm.rmD"
        """
        wigToBigWig -clip ${bedgraph}  ${chrom_sizes} ${out_prefix}.bigWig
        """
}

////////////////////////////////////////////////////
/* --          MACS2 BROAD                     -- */
////////////////////////////////////////////////////

// CALL BROAD PEAKS USING MACS2
process merge_sample_macs2 {

    tag "$sampleid"

    publishDir "${params.outdir}/align/mergeSample/macs2", mode: 'copy'

    input:
    set val(sampleid), file(bam) from merge_sample_rmdup_macs2_ch

    output:
    set val(sampleid), file("*/*.broadPeak") into merge_sample_macs2_homer_in_ch,
                                                  merge_sample_macs2_frip_in_ch,
                                                  merge_sample_macs2_merge_peaks_in_ch
    set val(sampleid), file("*/*.{gappedPeak,xls}") into merge_sample_macs2_output_ch
    set val(sampleid), file("*/*.sysout") into merge_sample_macs2_sysout_ch

    script:
        """
        mkdir -p ${sampleid}
        macs2 callpeak \\
              --verbose=2 \\
              -t ${bam[0]} \\
              --gsize=${params.macs_genome_size} \\
              --outdir=${sampleid} --name=${sampleid} \\
              --format BAMPE \\
              --keep-dup all \\
              --nomodel \\
              --broad \\
              >> ${sampleid}/${sampleid}.macs2.sysout 2>&1
        """
}

// ANNOTATE PEAKS USING HOMER ANNOTATEPEAKS
process merge_sample_macs2_homer {

    tag "$sampleid"

    publishDir "${params.outdir}/align/mergeSample/macs2", mode: 'copy'

    input:
    set val(sampleid), file(peak) from merge_sample_macs2_homer_in_ch
    file fasta from fasta_merge_sample_macs2_homer_ch.collect()
    file gtf from gtf_merge_sample_macs2_homer_ch.collect()

    output:
    set val(sampleid), file("*/*.annotatePeaks.txt") into merge_sample_macs2_homer_ch
    set val(sampleid), file("*/*.annotatePeaks.sysout") into merge_sample_macs2_homer_sysout_ch

    script:
        """
        mkdir -p ${sampleid}
        annotatePeaks.pl ${peak} \\
                         ${fasta} \\
                         -gid \\
                         -gtf ${gtf} \\
                         > ${sampleid}/${sampleid}_peaks.homer.annotatePeaks.txt \\
                         2> ${sampleid}/${sampleid}_peaks.homer.annotatePeaks.sysout
        """
}

// CALCULATE FRIP SCORE USING BEDTOOLS
process merge_sample_macs2_frip {

    tag "$sampleid"

    publishDir "${params.outdir}/align/mergeSample/macs2", mode: 'copy'

    input:
    set val(sampleid), file(peak), file(bam), file(flagstat) from merge_sample_macs2_frip_in_ch.join(merge_sample_rmdup_macs2_frip_ch, by: [0]).join(merge_sample_rmdup_flagstat_macs2_frip_ch, by: [0])

    output:
    set val(sampleid), file("*/*.frip.txt") into merge_sample_macs2_frip_ch

    script:
        """
        mkdir -p ${sampleid}
        READS_IN_PEAKS=\$(intersectBed -a ${bam[0]} -b ${peak} -bed -c -f 0.20 | awk -F '\t' '{sum += \$NF} END {print sum}')
        grep 'mapped (' ${flagstat} | awk -v a="\$READS_IN_PEAKS" '{print a/\$1}' > ${sampleid}/${sampleid}_peaks.frip.txt
        """
}

// GET LIST OF FRIP FILES ACROSS ALL SAMPLES FOR QC PLOT
merge_sample_macs2_frip_ch.map { it -> [ "macs2_frip", [it[1]] ] }
                          .groupTuple(by: [0])
                          .map { it ->  [ it[0],
                                          it[1].flatten().sort().collect{ it.getName().replace("_peaks.frip.txt","") },
                                          it[1].flatten().sort() ] }
                          .set { merge_sample_macs2_frip_collate_ch }

// GENERATE PLOTS FOR FRIP SCORES ACROSS ALL SAMPLES
process merge_sample_macs2_fripqc {

   tag "$name"

   publishDir "${params.outdir}/align/mergeSample/macs2/qc", mode: 'copy'

   input:
   set val(name), val(sampleids), file(frips) from merge_sample_macs2_frip_collate_ch

   output:
   file ("*.pdf") into merge_sample_macs2_fripqc_ch

   script:
       """
       Rscript $baseDir/bin/plot_frip.R -i ${frips.join(',')} -s ${sampleids.join(',')} -o ./ -p ${name}
       """
}

// GET LIST OF HOMER ANNOTATED FILES ACROSS ALL SAMPLES FOR QC PLOTS
merge_sample_macs2_homer_ch.map { it -> [ "macs2_homer_annotation", [it[1]] ] }
                           .groupTuple(by: [0])
                           .map { it ->  [ it[0],
                                           it[1].flatten().sort().collect{ it.getName().replace("_peaks.homer.annotatePeaks.txt","") },
                                           it[1].flatten().sort() ] }
                           .set { merge_sample_macs2_homer_collate_ch }

// GENERATE PLOTS FOR VARIOUS ASPECTS OF HOMER ANNOTATION ACROSS ALL SAMPLES
process merge_sample_macs2_homerqc {

   tag "$name"

   publishDir "${params.outdir}/align/mergeSample/macs2/qc", mode: 'copy'

   input:
   set val(name), val(sampleids), file(homers) from merge_sample_macs2_homer_collate_ch

   output:
   file ("*.pdf") into merge_sample_macs2_homerqc_ch

   script:
       """
       Rscript $baseDir/bin/homer_plot_annotation.R -i ${homers.join(',')} -s ${sampleids.join(',')} -o ./ -p ${name}
       """
}

// GET LIST OF PEAK ACROSS ALL SAMPLES FOR MERGING
merge_sample_macs2_merge_peaks_in_ch.map { it -> [ "merged_peaks", [ it[1]] ] }
                                    .groupTuple(by: [0])
                                    .map { it ->  [ it[0],
                                                    it[1].flatten().sort().collect{ it.getName().replace("_peaks.broadPeak","") },
                                                    it[1].flatten().sort() ] }
                                    .into { merge_sample_macs2_merge_peaks_in_ch;
                                            merge_sample_macs2_peakqc_in_ch }

// GENERATE PLOTS FOR VARIOUS ASPECTS OF PEAKS ACROSS ALL SAMPLES
process merge_sample_macs2_peakqc {

   tag "$name"

   publishDir "${params.outdir}/align/mergeSample/macs2/qc", mode: 'copy'

   input:
   set val(name), val(sampleids), file(peaks) from merge_sample_macs2_peakqc_in_ch

   output:
   file ("*.{txt,pdf}") into merge_sample_macs2_peakqc_ch

   script:
       """
       Rscript $baseDir/bin/macs2_peakqc.R -i ${peaks.join(',')} -s ${sampleids.join(',')} -o ./ -p macs2_peakqc
       """
}

// MERGE PEAKS ACROSS ALL SAMPLES AND CREATE A BOOLEAN OUTPUT FILE TO AID FILTERING
process merge_sample_macs2_merge_peaks {

   tag "$name"

   publishDir "${params.outdir}/align/mergeSample/macs2/merged_peaks", mode: 'copy'

   input:
   set val(name), val(sampleids), file(peaks) from merge_sample_macs2_merge_peaks_in_ch

   output:
   set val(name), file("*.boolean.txt") into merge_sample_macs2_merge_peaks_bool_ch
   set val(name), file("*.bed") into merge_sample_macs2_merge_peaks_bed_ch
   set val(name), file("*.saf") into merge_sample_macs2_merge_peaks_saf_ch
   set val(name), file("*.log") into merge_sample_macs2_merge_peaks_log_ch

   script:
       """
       sort -k1,1 -k2,2n ${peaks.join(' ')} \\
            | mergeBed -c 2,3,4,5,6,7,8,9 -o collapse,collapse,collapse,collapse,collapse,collapse,collapse,collapse > ${name}.txt
       python $baseDir/bin/expand_merged_macs.py ${name}.txt ${sampleids.join(',')} ${name}.boolean.txt --min_samples 1
       awk -v FS='\t' -v OFS='\t' 'FNR > 1 { print \$1, \$2, \$3, \$4, "0", "+" }' ${name}.boolean.txt > ${name}.bed
       echo -e "GeneID\tChr\tStart\tEnd\tStrand" > ${name}.saf
       awk -v FS='\t' -v OFS='\t' 'FNR > 1 { print \$4, \$1, \$2, \$3,  "+" }' ${name}.boolean.txt >> ${name}.saf
       """
}

// ANNOTATE PEAKS USING HOMER ANNOTATE AND ADD ANNOTATION TO BOOLEAN OUTPUT FILE
process merge_sample_macs2_merge_peaks_homer {

    tag "${name}"

    publishDir "${params.outdir}/align/mergeSample/macs2/merged_peaks", mode: 'copy'

    input:
    set val(name), file(bed) from merge_sample_macs2_merge_peaks_bed_ch
    set val(name), file(bool) from merge_sample_macs2_merge_peaks_bool_ch
    file fasta from fasta_merge_sample_macs2_merge_peaks_homer_ch.collect()
    file gtf from gtf_merge_sample_macs2_merge_peaks_homer_ch.collect()

    output:
    set val(name), file("*.annotatePeaks.txt") into merge_sample_macs2_merge_peaks_homer_ch
    set val(name), file("*.annotatePeaks.sysout") into merge_sample_macs2_merge_peaks_homer_sysout_ch

    script:
        """
        annotatePeaks.pl ${bed} \\
                         ${fasta} \\
                         -gid \\
                         -gtf ${gtf} \\
                         > ${name}.homer.annotatePeaks.txt \\
                         2> ${name}.homer.annotatePeaks.sysout
        cut -f2- ${name}.homer.annotatePeaks.txt | awk 'NR==1; NR > 1 {print \$0 | "sort -k1,1 -k2,2n"}' | cut -f6- > tmp.txt
        paste ${bool} tmp.txt > ${bool.toString().replace(".txt","")}.homer.annotatePeaks.txt
        """
}

// GET LIST OF REPLICATE LEVEL BAM FILES
merge_replicate_name_bam_merge_sample_featurecounts_ch.map { it -> [ "merged_bam", [ it[1]] ]}
                                                      .groupTuple(by: [0])
                                                      .map { it ->  [ it[1].flatten().sort().collect{ it.getName().replace(".mRp.rmD.bam","") },
                                                                      it[1].flatten().sort() ] }
                                                      .set { merge_replicate_name_bam_merge_sample_featurecounts_ch }

process merge_sample_macs2_merge_peaks_featurecounts {

    tag "${name}"

    label 'highcpu'

    publishDir "${params.outdir}/align/mergeSample/macs2/merged_peaks", mode: 'copy'

    input:
    set val (sampleids), file(bams) from merge_replicate_name_bam_merge_sample_featurecounts_ch
    set val(name), file(saf) from merge_sample_macs2_merge_peaks_saf_ch.collect()

    output:
    set val(name), file("*.featureCounts.txt") into merge_sample_macs2_merge_peaks_featurecounts_ch
    set val(name), file("*.featureCounts.sysout") into merge_sample_macs2_merge_peaks_featurecounts_sysout_ch

    script:
        """
        featureCounts -F SAF \\
                      -O \\
                      -T ${task.cpus} \\
                      --fracOverlap 0.5 \\
                      --primary \\
                      --ignoreDup \\
                      -p \\
                      -B \\
                      -C \\
                      --donotsort \\
                      -a ${saf} \\
                      -o ${name}.replicate.featureCounts.txt \\
                      ${bams.join(' ')} \\
                      >> ${name}.replicate.featureCounts.sysout 2>&1
        """
}

process merge_sample_macs2_merge_peaks_differential {

    tag "${name}"

    publishDir "${params.outdir}/align/mergeSample/macs2/merged_peaks", mode: 'copy'

    input:
    set val(name), file(counts) from merge_sample_macs2_merge_peaks_featurecounts_ch

    output:
    set val(name), file("deseq2/*") into merge_sample_macs2_merge_peaks_differential_ch
    file "*.txt" into merge_sample_macs2_merge_peaks_differential_complete_ch

    script:
        """
        Rscript $baseDir/bin/featurecounts_deseq2.R -i ${counts} -b '.mRp.rmD.bam' -o ./deseq2 -p ${name}
        touch merge_sample_macs2_merge_peaks_differential.complete.txt
        """
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                             IGV                                     -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

// MERGE CHANNELS FOR REPLICATE AND SAMPLE LEVEL FILES TO AVOID RUNNING PROCESS PER SAMPLE
merge_replicate_bigwig_ch.map { it -> [ 'bigwig', it[1] ] }
                         .groupTuple()
                         .map { it -> it[1].sort() }
                         .set { merge_replicate_bigwig_ch }

merge_sample_bigwig_ch.map { it -> [ 'bigwig', it[1] ] }
                      .groupTuple()
                      .map { it -> it[1].sort() }
                      .set { merge_sample_bigwig_ch }

// HAD TO WRITE CUSTOM SCRIPT TO REORDER TRACKS IN SESSION FILE MORE SENSIBLY.
// PROBABLY POSSIBLE WITH NEXTFLOW BUT MUCH EASIER IN PYTHON!
process igv_session {

    tag "igv_session"

    publishDir "${params.outdir}/igv", mode: 'copy'

    input:
    file bigwigs from merge_replicate_bigwig_ch.merge(merge_sample_bigwig_ch)
    file fasta from fasta_igv_ch.collect()
    file gtf from gtf_igv_ch.collect()
    file replicate_diff from merge_replicate_macs2_merge_peaks_differential_complete_ch.collect()
    file sample_diff from merge_sample_macs2_merge_peaks_differential_complete_ch.collect()

    output:
    file "*.{xml,txt}" into igv_session_ch

    script:
        """
        [ ! -f ${params.outdir_abspath}/genome/${fasta.getName()} ] && ln -s ${params.fasta} ${params.outdir_abspath}/genome/${fasta.getName()}
        [ ! -f ${params.outdir_abspath}/genome/${gtf.getName()} ] && ln -s ${params.gtf} ${params.outdir_abspath}/genome/${gtf.getName()}
        python $baseDir/bin/get_files_for_igv.py ${params.outdir_abspath} igv_files.txt
        python $baseDir/bin/files_to_igv_session.py igv_session.xml igv_files.txt ${params.outdir_abspath}/genome/${fasta.getName()}
        """
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                          MULTIQC                                    -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

// RUN THIS AFTER IGV BECAUSE ITS AT THE END OF THE PIPELINE
process multiqc {

    tag "multiqc"

    publishDir "${params.outdir}/qc/multiqc", mode: 'copy'

    input:
    file igvs from igv_session_ch

    output:
    file "*" into multiqc_ch

    script:
        """
        multiqc ${params.outdir_abspath} \\
                -f \\
                --config ${params.multiqc_config} \\
                --filename BABS-ATACSeqPE_multiqc_report.html \\
                -m fastqc \\
                -m fastq_screen \\
                -m cutadapt \\
                -m samtools \\
                -m picard
        """
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                        END OF PIPELINE                              -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
