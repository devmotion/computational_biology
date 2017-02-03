#!/usr/bin/env Rscript

# Authors: Rebekka Mueller, David Widmann
# 2016-12-27
# v0.5

"Run AneuFinder with variable width bins and blacklist of euploid reference.

Usage:
  aneufinder.R [--nocluster --ncpu <ncpu> --ntrials <ntrials> --mchrom <mchrom> --binsize <binsize>...] <output> [<input>]
  aneufinder.R -h | --help
  aneufinder.R --version

Arguments:
  <output>  Output folder.
  <input>   Input folder of BAM files [default: .].

Options:
  -h, --help                         Show this screen.
  --version                          Show version.
  --nocluster                        Do not cluster plots by similarity. If the input directory contains only one sample, this option will be set automatically.
  -c <ncpu>, --ncpu <ncpu>           Number of CPUs [default: 2].
  -t <ntrials>, --ntrials <ntrials>  Number of trials to find a fit where state '2-somy' is most frequent [default: 15].
  -m <mchrom>, --mchrom <mchrom>     Maximum number of chromosomes appearing in the Hidden Markov Model CNV states 'c('zero-inflation','0-somy','1-somy','2-somy','3-somy', ...)'. It has to be greater or equal than 2 since state '2-somy' is expected to be most frequent [default: 10].
  -b <binsize>, --binsize <binsize>  Size of bins (in bp) [default: 1e6].
" -> doc

suppressPackageStartupMessages(library(docopt))

opt <- docopt(doc, version="0.5")

# Default options for arguments not working in docopt (I guess)
if (is.null(opt$input))
  opt$input <- "."

cat("Options:\n")
cat("\toutput =", opt$output, "\n")
cat("\tinput =", opt$input, "\n")
cat("\tnocluster =", opt$nocluster, "\n")
cat("\tncpu =", opt$ncpu, "\n")
cat("\tntrials = ", opt$ntrials, "\n")
cat("\tmchrom = ", opt$mchrom, "\n")
cat("\tbinsize =", opt$binsize, "\n\n")

# Input check
mchrom <- as.integer(opt$mchrom)
if (mchrom < 2)
  stop("Please specify a maximum number of chromosomes greater or equal than 2.")

# List input files
cat("Sample files in input directory:\n")
files <- list.files(opt$input, pattern = '\\.(bam|bed(\\.gz)?)$', ignore.case = TRUE)
print(files)

# Do not cluster plots for only one sample
if (length(files) < 2 & !opt$nocluster) {
  opt$nocluster <- TRUE
  warning("Do not cluster plots by similarity since only one sample is found in input directory.")
}

# Load blacklist of euploid reference from AneuFinderData package
blacklist <- system.file("extdata", "blacklist-hg19.bed.gz", package="AneuFinderData")

# Load mappability correction reference from AneuFinderData package
var.width.ref <- system.file("extdata", "hg19_diploid.bam.bed.gz", package="AneuFinderData")

# Load full genome sequence for GC correction
suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg19))

# Run AneuFinder
suppressPackageStartupMessages(library(AneuFinder))
cnv.states <- c('zero-inflation', paste0(0:mchrom, '-somy'))
Aneufinder(inputfolder = opt$input, outputfolder = opt$output,
           assembly = 'hg19', numCPU = opt$ncpu,
           binsizes = as.numeric(opt$binsize),
           variable.width.reference = var.width.ref,
           chromosomes = c(1:22,'X'), blacklist = blacklist,
           states = cnv.states, correction.method = 'GC',
           GC.BSgenome = BSgenome.Hsapiens.UCSC.hg19,
           use.bamsignals = TRUE, max.time = -1,
           num.trials = as.integer(opt$ntrials), cluster.plots = !opt$nocluster)
