#!/usr/bin/env Rscript

# Authors: Rebekka Mueller, David Widmann
# 2017-01-14
# v0.1

"Plot frequencies of combinations for each state resulting from chromstaR output.

Usage:
  plotfrequency.R  <inputfolder> [<output>]
  plotfrequency.R -h | --help
  plotfrequency.R --version

Arguments:
  <inputfolder>            Folder with chromstaR output.
  <output>                 Output directory for [default: .].

Options:
  -h, --help               Show this screen.
  --version                Show version.
" -> doc

suppressPackageStartupMessages(library(docopt))

opt <- docopt(doc, version="0.1")

# Default options for arguments not working in docopt (I guess)
if (is.null(opt$output))
  opt$output <- "."

cat("Options:\n")
cat("\tinputfolder =", opt$inputfolder, "\n")
cat("\toutput =", opt$output, "\n\n")

# Check if file exists
if (!dir.exists(opt$inputfolder))
  stop("Input folder ", opt$inputfolder, " does not exist.")

cat("Loading chromstaR ...\n")
suppressPackageStartupMessages(library(chromstaR))

# Local variables for states, combinations, conditions, and frequencies
states <- character()
combinations <- character()
conditions <- character()
frequencies <- numeric()

for (state in list.dirs(path=opt$inputfolder, full.names=FALSE, recursive=FALSE)) {
  cat("Loading combined model for", state, "...\n")

  # Load combined model
  model.file <- file.path(opt$inputfolder, state, 'combined', 'combined_mode-full.RData')
  if (!file.exists(model.file))
    next
  combined.model <- loadHmmsFromFiles(model.file, check.class = "combinedMultiHMM")[[1]]

  state.with.bins <- paste0(state, '\n(', length(combined.model$bins), ')')

  for (condition in levels(factor(combined.model$info$condition))) {
    states <- c(states, rep(state.with.bins, nrow(combined.model$frequencies)))
    combinations <- c(combinations, as.character(combined.model$frequencies[[paste0('combination.', condition)]]))
    conditions <- c(conditions, rep(condition, nrow(combined.model$frequencies)))
    frequencies <- c(frequencies, combined.model$frequencies$frequency)
  }
}

# Create data.frame of states, combinations, conditions, and frequencies
freq <- data.frame(state=states, combination=combinations, condition=conditions, frequency=frequencies)

# Correct order of states
poly.state <- startsWith(levels(freq$state), '>')
freq$state <- factor(freq$state, levels=c(sort(levels(freq$state)[!poly.state]), sort(levels(freq$state)[poly.state])))

# Split table per condition
freq.per.condition <- split(freq, freq$condition)

# For every condition
for (condition in names(freq.per.condition)) {
  # Select frequencies of condition
  data <- freq.per.condition[[condition]]

  # Output
  savename <- file.path(opt$output, paste0('frequency_per_state_', condition, '.pdf'))

  cat("Saving plot in", savename, "...\n")

  if (file.exists(savename))
    warning(savename, " already exists and is overwritten.")

  # Visualize the distributions for each state
  ggplt <- ggplot(data=data, aes(x=state, y=frequency, fill=combination)) +
    geom_bar(stat="identity") + xlab("CNV state") + theme_bw()

  ggsave(savename, plot=ggplt, width=20, height=15, limitsize=FALSE, units='cm')
}
