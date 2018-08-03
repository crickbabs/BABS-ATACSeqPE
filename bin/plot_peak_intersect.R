#!/usr/bin/env Rscript

################################################
################################################
## LOAD LIBRARIES                             ##
################################################
################################################

library(optparse)
library(UpSetR)

################################################
################################################
## PARSE COMMAND-LINE PARAMETERS              ##
################################################
################################################

option_list <- list(make_option(c("-i", "--frip_files"), type="character", default=NULL, help="Comma-separated list of FRiP score files. Each file should contain one line with FRiP score. ", metavar="frip_files"),
										make_option(c("-s", "--sample_ids"), type="character", default=NULL, help="Comma-separated list of sample ids associated with FRiP files. Must be unique and in same order as FRiP files input.", metavar="sampleids"),
									  make_option(c("-o", "--outdir"), type="character", default='./', help="Output directory", metavar="path"),
									  make_option(c("-p", "--outprefix"), type="character", default='plot_frip', help="Output prefix", metavar="character"))


option_list <- list(make_option(c("-i", "--input_file"), type="character", default=NULL, help="Path to tab-delimited file containing two columns i.e sample1&sample2&sample3 indicating intersect between samples <TAB> set size.", metavar="infile"),
										make_option(c("-o", "--output_file"), type="character", default=NULL, help="Path to output file with '.pdf' extension.", metavar="outfile"))

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$infile)){
		print_help(opt_parser)
		stop("Input file must be supplied.", call.=FALSE)
}
if (is.null(opt$outfile)){
		print_help(opt_parser)
		stop("Output pdf file must be supplied.", call.=FALSE)
}

OutDir <- dirname(opt$outfile)
if (file.exists(OutDir) == FALSE) {
    dir.create(OutDir,recursive=TRUE)
}

################################################
################################################
## PLOT DATA                                  ##
################################################
################################################

comb.dat <- read.table(opt$infile,sep="\t",header=FALSE)
comb.vec <- comb.dat[,2]
comb.vec <- setNames(comb.vec,comb.dat[,1])

pdf(opt$outfile,height=7,width=12)

upset(fromExpression(comb.vec),
			sets.bar.color = "#56B4E9",
			point.size = 5,
			line.size = 2,
			order.by = "freq",
      text.scale = c(1.7, 1.5, 1.7, 1.5, 1.7, 1.7))

dev.off()

################################################
################################################
################################################
################################################
