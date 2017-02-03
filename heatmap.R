#!/usr/bin/env Rscript

# Authors: Rebekka Mueller, David Widmann
# 2017-01-10
# v0.1

"Plot heatmap of AneuFinder models.

Usage:
  heatmap.R [options] <model>...
  heatmap.R -h | --help
  heatmap.R --version

Arguments:
  <model>                AneuFinder model (dnacopy or HMM).

Options:
  -h, --help             Show this screen.
  --version              Show version.
  -f, --force            Force to overwrite the output files if present.
  -o, --output <output>  Output file [default: ./genomeHeatmap.pdf].
" -> doc

suppressPackageStartupMessages(library(docopt))

opt <- docopt(doc, version="0.1")

cat("Options:\n")
cat("\tmodel =", opt$model, "\n")
cat("\tforce =", opt$force, "\n")
cat("\toutput =", opt$output, "\n\n")

# Set up output directory
if (!dir.exists(dirname(opt$output)))
  dir.create(dirname(opt$output), recursive=TRUE)

if (!opt$force && file.exists(opt$output))
  stop("Output file already exists. Specify '-f' or '--force' on the commandline to overwrite it.")

# Load AneuFinder
cat("Loading AneuFinder ...\n")
suppressPackageStartupMessages(library(AneuFinder))

# Plot genomewide heatmap without clustering
heatmapGenomewide(rev(opt$model), cluster=FALSE, file=opt$output)
