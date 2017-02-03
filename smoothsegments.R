#!/usr/bin/env Rscript

# Authors: Rebekka Mueller, David Widmann
# 2017-01-09
# v0.1

"Smooth segments in AneuFinder results.

Usage:
  smoothsegments.R [options] <model> [<output>]
  smoothsegments.R -h | --help
  smoothsegments.R --version

Arguments:
  <model>  AneuFinder model (dnacopy or HMM).
  <output>  Output directory for smoothed AneuFinder model [default: .].

Options:
  -h, --help  Show this screen.
  --version  Show version.
  -c, --cutoff <cutoff>  Segments consisting of a number of bins below or equal to this threshold are merged with neighbouring segment, if possible [default: 10].
  -d, --dry  Dry run. Do not apply or save changes.
" -> doc

suppressPackageStartupMessages(library(docopt))

opt <- docopt(doc, version="0.1")

# Default options for arguments not working in docopt (I guess)
if (is.null(opt$output))
  opt$output <- "."

cat("Options:\n")
cat("\tmodel = ", opt$model, "\n")
cat("\toutput =", opt$output, "\n")
cat("\tcutoff =", opt$cutoff, "\n")
cat("\tdry =", opt$dry, "\n\n")

# Check if file exists
if (!file.exists(opt$model))
  stop("File ", opt$model, " does not exist.")

# Convert cutoff to number
cutoff <- as.numeric(opt$cutoff)

# Try to load AneuFinder model
tryCatch ({
  cat("Loading AneuFinder model ...\n")
  suppressPackageStartupMessages(library(AneuFinder))
  model <- loadFromFiles(opt$model, check.class="aneuHMM")[[1]]
}, error = function(err) {
  stop("Could not load model ", opt$model, ". Error: ", err)
})

# Check that model was successfully loaded
if (is.null(model$bins) | class(model$bins) != 'GRanges' | length(model$bins) == 0)
    stop("Model does not contain any bins.")


cat("Computing segmentation ...\n")

# Split states into groups of chromosomes
state.per.chrom <- lapply(split(as.character(model$bins$state), as.factor(seqnames(model$bins))),
                          as, 'Rle')

# Merge segments (runs) below threshold for every chromosome
for (chrom.name in names(state.per.chrom)) {
  chrom <- state.per.chrom[[chrom.name]]
  chrom.nrun <- nrun(chrom)
  chrom.length <- runLength(chrom)

  cat("\nSegmentation for chromosome ", chrom.name, ":\n", sep="")
  print(chrom)

  min.index <- which.min(chrom.length)
  updated <- FALSE
  while (chrom.nrun > 1 & chrom.length[min.index] <= cutoff) {
    cat("Merging segment of state", runValue(chrom)[min.index],
        "consisting of", chrom.length[min.index], "bins ...\n")
    updated <- TRUE

    # Calculate index of neighbour with maximal length
    if (min.index == 1) {
      max.index <- 2
    } else if (min.index == chrom.nrun) {
      max.index <- chrom.nrun-1
    } else {
      max.index <- min.index - 3 + 2*which.max(chrom.length[c(-1,1) + min.index])
    }

    # Update length of runs
    chrom.length[max.index] <- chrom.length[max.index] + chrom.length[min.index]
    chrom.length[min.index] <- 0

    # Update Rle object
    runLength(chrom) <- chrom.length

    # Update local variables
    chrom.nrun <- nrun(chrom)
    chrom.length <- runLength(chrom)
    min.index <- which.min(chrom.length)
  }

  # Update list of chromosomes
  if (updated) {
    if (!opt$dry) state.per.chrom[[chrom.name]] <- chrom
    cat("New segmentation for chromosome ", chrom.name, ":\n", sep="")
    print(chrom)
  }
}

if (opt$dry) {
  cat("\nDry run completed.\n")
} else {
  # Write new state to bins
  cat("\nUpdating bins ...\n")
  model$bins$state <- factor(unsplit(lapply(state.per.chrom, as.vector),
                                     as.factor(seqnames(model$bins))))
  model$bins$copy.number <- as.numeric(sub('^(\\d+)-somy$', '\\1', model$bins$state))
  names(model$bins$copy.number) <- model$bins$state

  # Update segments
  cat("Updating segments ...\n")
  model$segments <- as(collapseBins(as.data.frame(model$bins), column2collapseBy='state', columns2drop='width', columns2average=c('counts','mcounts','pcounts')), 'GRanges')
  seqlevels(model$segments) <- seqlevels(model$bins) # correct order from as()
  seqlengths(model$segments) <- seqlengths(model$bins)[names(seqlengths(model$segments))]

  # Update ID
  cat("Updating ID ...\n")
  model$ID <- sub("^([^\\.]*)\\.(.*)$", paste0("\\1", "_cutoff-", opt$cutoff, ".\\2"), model$ID)

  # Update quality information
  cat("Updating quality information ...\n")
  model$qualityInfo <- as.list(getQC(model))

  # Save model to output directory
  cat("Saving smoothed model to output directory ...\n")
  if (!dir.exists(opt$output))
    dir.create(opt$output, recursive=TRUE)

  model.method <- match.call(definition = findCNVs <- function(method, ...) {},
                             attr(model, "call"))$method
  savename.model <- file.path(opt$output,
                              paste0("smoothed_", model.method, "-",
                                     sub('^(.*)\\.bam_binsize_.*$', '\\1', basename(opt$model), ignore.case=TRUE),
                                     "_cutoff-", opt$cutoff, ".RData"))

  if (file.exists(savename.model))
    warning(savename.model, " already exists and is overwritten.")

  save(model, file=savename.model)

  # Save profile plot to output directory
  cat("Saving profile plot to output directory ...\n")
  suppressPackageStartupMessages(library(ggplot2))

  savename.plot <- file.path(opt$output,
                             paste0("profile_smoothed_", model.method, "-",
                                    sub('^(.*)\\.bam_binsize_.*$', '\\1', basename(opt$model), ignore.case=TRUE),
                                    "_cutoff-", opt$cutoff, ".pdf"))

  if (file.exists(savename.plot))
    warning(savename.plot, " already exists and is overwritten.")

  pdf(file=savename.plot, width=20, height=10)
  print(plot(model, type='profile'))
  invisible(dev.off())

  cat("Done.\n")
}
