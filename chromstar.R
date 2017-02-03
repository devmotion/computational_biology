#!/usr/bin/env Rscript

# Authors: Rebekka Mueller, David Widmann
# 2017-01-25
# v0.2

"Run chromstaR on copy number segments for patient SBT15003 or SBT15008 obtained from AneuFinder.

Usage:
  chromstar.R [options] <model> <chipseq>
  chromstar.R -h | --help
  chromstar.R --version

Arguments:
  <model>                  AneuFinder model (dnacopy or HMM).
  <chipseq>                Folder with ChIP-seq files from MS160621 run.

Options:
  -h, --help               Show this screen.
  --version                Show version.
  -b, --binsize <binsize>  Approximate size of fixed-width bins [default: 1000].
  --exact                  Truncate ends of segments to obtain bins of exactly the same size.
  -c, --ncpu <ncpu>        Number of CPUs [default: 2].
  -o, --output <output>    Output directory for chromstaR runs [default: .].
" -> doc

suppressPackageStartupMessages(library(docopt))

opt <- docopt(doc, version="0.2")

cat("Options:\n")
cat("\tmodel =", opt$model, "\n")
cat("\tchipseq =", opt$chipseq, "\n")
cat("\tbinsize =", opt$binsize, "\n")
cat("\texact =", opt$exact, "\n")
cat("\tncpu =", opt$ncpu, "\n")
cat("\toutput =", opt$output, "\n\n")

# Check if file exists
if (!file.exists(opt$model))
  stop("Model ", opt$model, " does not exist.")

if (!dir.exists(opt$chipseq))
  stop("ChIP-seq folder ", opt$chipseq, " does not exist.")

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

# Check that model belongs to one of the conditions 'SBT003' or 'SBT008'
model.cond <- c('SBT15003','SBT15008')[startsWith(model$ID, c('BB150521_I_', 'BB151124_I_'))]
if (length(model.cond) != 1)
  stop("Model does not belong to one of the conditions 'SBT15003' or 'SBT15008'.")

# Change copy numbers > 4 to state '>4-somy'
levels(model$bins$state)[!(levels(model$bins$state) %in% paste0(0:4, '-somy'))] <- '>4-somy'
model$segments <- as(collapseBins(as.data.frame(model$bins), column2collapseBy='state',
                                  columns2drop='width', columns2average=c('counts','mcounts','pcounts')),
                     'GRanges')
seqlevels(model$segments) <- seqlevels(model$bins) # correct order from as()
seqlengths(model$segments) <- seqlengths(model$bins)[names(seqlengths(model$segments))]

# Create bins
cat("Creating fixed-width bins of segments for each state ...\n")
binsize <- as.integer(opt$binsize)
fixedwidth.bins <- lapply(split(model$segments, as.factor(model$segments$state)),
                          function(x) {
                            if (opt$exact)
                              end(x) <- start(x)-1+floor(width(x)/binsize)*binsize
                            unlist(tile(x, width=binsize))
                          })

# Load chromstaR
cat("\nLoading chromstaR ...\n")
suppressPackageStartupMessages(library(chromstaR))
suppressPackageStartupMessages(library(rtracklayer))

# Load lengths of human chromosomes
suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg19))
seqlengths.hg19 <- seqlengths(Hsapiens)

# Set up color mapping
colors <- c('#89a651', '#965ea7', '#bc5842')
names(colors) <- c('[H3K9me3]', '[H3K27me3]', '[H3K27me3+H3K9me3]')

# Set up experiment table
exp <- data.frame(file=c('MS160621_002.bam','MS160621_004.bam','MS160621_007.bam','MS160621_009.bam'),
                  mark=rep(c('H3K27me3','H3K9me3'), 2),
                  condition=rep(c('SBT15003','SBT15008'), each=2),
                  replicate=1, pairedEndReads=FALSE,
                  controlFiles=rep(c('MS160621_001.bam','MS160621_006.bam'), each=2))
exp <- subset(exp, condition==model.cond)

# Define export functions

# Export browser files of read counts and combinations
export.read_combinations <- function(directory) {
  cat("Exporting browser files of read counts and combinations ...\n")

  # Obtain all combined models and define output prefixes
  models <- list.files(path=file.path(directory, 'combined'), pattern="\\.RData$", full.names=TRUE)
  browserfiles <- sub(file.path('^(.*)', 'combined', '(.*)\\.RData$'),
                      file.path('\\1', 'BROWSERFILES', '\\2'),
                      models)

  for (i in seq_along(models)) {
    # Load combined model
    model <- loadHmmsFromFiles(models[i], check.class='combinedMultiHMM')[[1]]
    model.bins <- granges(model$bins)
    model.counts <- model$bins$counts
    model.IDs <- model$info$ID
    model.conds <- names(model$segments.per.condition)

    # Convert chromosomes to UCSC format
    seqlevelsStyle(model.bins) <- 'UCSC'

    # Update lengths of chromosomes
    genome(model.bins) <- 'hg19'
    seqlengths(model.bins) <- seqlengths.hg19[names(seqlengths.hg19) %in% seqlevels(model.bins)]

    # Export counts of every model
    for (j in seq_along(model.IDs)) {
      ID <- model.IDs[j]

      model.bins$score <- model.counts[,j]
      export.bw(object=model.bins, con=paste0(browserfiles[i], '_counts_', ID, '.bw'))
    }

    # Export combinations of every condition
    for (j in seq_along(model.conds)) {
      cond <- model.conds[j]

      model.combinations <- model$segments.per.condition[[j]]

      # Remove empty combinations
      model.combinations <- subset(model.combinations, combination != '[]')

      # Convert chromosomes to UCSC format
      seqlevelsStyle(model.combinations) <- 'UCSC'

      # Update genome
      genome(model.combinations) <- 'hg19'

      # Rename metadata column of combinations
      model.combinations$Name <- model.combinations$combination
      model.combinations$combination <- NULL

      # Add type and colour column
      model.combinations$type <- 'ChIP_seq_region'
      model.combinations$color <- colors[as.character(model.combinations$Name)]

      # Export combinations
      output <- paste0(browserfiles[i], '_combinations_', cond, '.gff3')
      export.gff3(object=model.combinations, con=output, source='chromstaR')
      system(command=paste('gzip', output))
    }
  }
}

# Export browser files of univariate peaks
export.peaks <- function(directory) {
  cat("Exporting browser files of univariate peaks ...\n")

  # Obtain all univariate models and define output prefixes
  models <- list.files(path=file.path(directory , 'univariate'), pattern="\\.RData$", full.names=TRUE)
  browserfiles <- sub(file.path('^(.*)', 'univariate', '(.*)\\.RData$'),
                      file.path('\\1', 'BROWSERFILES', paste0('univariate_peaks_', '\\2', '.gff3')),
                      models)

  for (i in seq_along(models)) {
    # Load univariate model
    model <- loadHmmsFromFiles(models[i], check.class='uniHMM')[[1]]
    model.peaks <- model$peaks

    # Convert chromosomes to UCSC format
    seqlevelsStyle(model.peaks) <- 'UCSC'

    # Update genome
    genome(model.peaks) <- 'hg19'

    # Rename metadata column of peakscores
    model.peaks$score <- model.peaks$peakScores
    model.peaks$peakScores <- NULL

    # Add type column
    model.peaks$type <- 'ChIP_seq_region'
    model.peaks$ID <- paste0('peak_', seq_along(model.peaks))

    # Export peaks
    export.gff3(object=model.peaks, con=browserfiles[i], source='chromstaR')
    system(command=paste('gzip', browserfiles[i]))
  }
}

# Set up output directory
if (!dir.exists(opt$output))
  dir.create(opt$output, recursive=TRUE)

# Run chromstaR for every copy number state
for (state in names(fixedwidth.bins)) {
  cat("\nRunning chromstaR for state", state, "...\n")

  bins <- list()
  bins[[as.character(binsize)]] <- fixedwidth.bins[[state]]
  output <- file.path(opt$output, state)

  Chromstar(inputfolder=opt$chipseq, experiment.table=exp,
            outputfolder=output, mode='full',
            numCPU=as.integer(opt$ncpu), pre.bins=bins)


  # Export read counts and combinations
  export.read_combinations(output)

  # Export univariate peaks
  export.peaks(output)
}
