#!/usr/bin/env Rscript

# Authors: Rebekka Mueller, David Widmann
# 2017-01-09
# v0.1

"Print segment of smallest size in AneuFinder results.

Usage:
  shortestsegment.R [--bins] <model>...
  shortestsegment.R -h | --help
  shortestsegment.R --version

Options:
  -h, --help  Show this screen.
  --version  Show version.
  -b, --bins  Print segment with smallest number of overlapped bins.

Arguments:
  <model>  AneuFinder model (dnacopy or HMM).
" -> doc

suppressPackageStartupMessages(library(docopt))

opt <- docopt(doc, version="0.1")

suppressPackageStartupMessages(library(GenomicRanges))

cat("width", "bins", "chrom", "state", "model", sep="\t")
cat("\n")

for (model in opt$model) {
  # Check if file exists
  if (!file.exists(model)) {
    warning("File ", model, " does not exist.")
    next
  }

  tryCatch ({
    # Load AneuFinder model in new temporary environment
    tmpenv <- new.env()
    load(model, envir=tmpenv)

    # Calculate and output minimal segment
    if (opt$bins) {
      overlaps <- countOverlaps(tmpenv$model$segments, tmpenv$model$bins)
      min.overlaps <- tmpenv$model$segments[overlaps==min(overlaps)]
      min.segment <- min.overlaps[which.min(width(min.overlaps))]
    } else {
      widths <- width(tmpenv$model$segments)
      min.widths <- tmpenv$model$segments[widths==min(widths)]
      min.segment <- min.widths[which.min(countOverlaps(tmpenv$model$segments, min.widths))]
    }
    cat(width(min.segment), countOverlaps(min.segment, tmpenv$model$bins), as.character(seqnames(min.segment)), as.character(min.segment$state), model, sep="\t")
    cat("\n")
  }, error = function(err) {
    warning("Could not load file ", model, ". Error: ", err)
  }, finally = {
    # Delete temporary environment
    rm(tmpenv)
  })
}
