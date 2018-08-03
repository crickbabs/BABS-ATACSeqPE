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

option_list <- list(make_option(c("-i", "--input_file"), type="character", default=NULL, help="Path to tab-delimited file containing two columns i.e sample1&sample2&sample3 indicating intersect between samples <TAB> set size.", metavar="input_file"),
										make_option(c("-o", "--output_file"), type="character", default=NULL, help="Path to output file with '.pdf' extension.", metavar="output_file"))

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input_file)){
		print_help(opt_parser)
		stop("Input file must be supplied.", call.=FALSE)
}
if (is.null(opt$output_file)){
		print_help(opt_parser)
		stop("Output pdf file must be supplied.", call.=FALSE)
}

OutDir <- dirname(opt$output_file)
if (file.exists(OutDir) == FALSE) {
    dir.create(OutDir,recursive=TRUE)
}

################################################
################################################
## PLOT DATA                                  ##
################################################
################################################

comb.dat <- read.table(opt$input_file,sep="\t",header=FALSE)
comb.vec <- comb.dat[,2]
comb.vec <- setNames(comb.vec,comb.dat[,1])

pdf(opt$output_file,onefile=F,height=10,width=14)

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
