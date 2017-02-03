#!/usr/bin/env Rscript

# Authors: Rebekka Mueller, David Widmann
# 2016-12-31
# v0.2

"Plot ChIP-seq read counts of chromstaR runs for different copy number states.

Usage:
  plotreads.R [options] <chromstar>
  plotreads.R -h | --help
  plotreads.R --version

Arguments:
  <chromstar>            Folder of chromstaR runs for different copy number states.

Options:
  -h, --help             Show this screen.
  --version              Show version.
  --outliers             Plot outliers.
  --nostats              Do not plot statistics.
  --raw                  Plot raw read counts. By default normalized read counts are plotted.
  -o, --output <output>  Output folder [default: .].
" -> doc

suppressPackageStartupMessages(library(docopt))

opt <- docopt(doc, version="0.2")

cat("Arguments:\n")
cat("\tchromstar =", opt$chromstar, "\n")
cat("\toutliers =", opt$outliers, "\n")
cat("\tnostats =", opt$nostats, "\n")
cat("\traw =", opt$raw, "\n")
cat("\toutput =", opt$output, "\n\n")

# Check if specified file or directory exists
if (!dir.exists(opt$chromstar))
  stop("Please specify a directory of chromstaR runs for different copy number states.")

# Load chromstaR models
cat("Loading chromstaR models ...\n")
suppressPackageStartupMessages(library(chromstaR))

# Local variables for states, experiment IDs, and counts
states <- character()
IDs <- character()
counts <- integer()

# Load counts for every state
for (state in list.dirs(path=opt$chromstar, full.names=FALSE, recursive=FALSE)) {
  cat("Loading counts for state", state, "...\n")

  if (opt$raw) {
    # Load binned data
    binned.files <- list.files(path=file.path(opt$chromstar, state, 'binned'), pattern='\\.RData$', full.names=TRUE)

    for (binned.file in binned.files) {
      temp.env <- new.env()
      bins <- get(load(binned.file, envir=temp.env), envir=temp.env)

      states <- c(states, rep(state, length(bins)))
      IDs <- c(IDs, rep(sub('^(.*)_binsize.*\\.RData$', '\\1', basename(binned.file)), length(bins)))
      counts <- c(counts, bins$counts)
    }
  } else {
    # Load combined model
    model.file <- file.path(opt$chromstar, state, 'combined', 'combined_mode-full.RData')
    if (!file.exists(model.file))
      next
    model <- loadHmmsFromFiles(model.file, check.class='combinedMultiHMM')[[1]]

    model.IDs <- model$info$ID

    for (i in seq_along(model.IDs)) {
      states <- c(states, rep(state, length(model$bins)))
      IDs <- c(IDs, rep(model.IDs[i], length(model$bins)))
      counts <- c(counts, model$bins$counts[,i])
    }
  }
}

# Create data.frame of states, experiment IDs, and counts
read.counts <- data.frame(state=states, ID=IDs, counts=counts)

# Correct order of states
poly.state <- startsWith(levels(read.counts$state), '>')
read.counts$state <- factor(read.counts$state, levels=c(sort(levels(read.counts$state)[!poly.state]), sort(levels(read.counts$state)[poly.state])))

# Create plots
library(ggplot2)

if (!dir.exists(opt$output))
  dir.create(opt$output, recursive=TRUE)

# Split table per experiment
read.counts.per.ID <- split(read.counts, read.counts$ID)

# For every ChIP-seq experiment
for (ID in names(read.counts.per.ID)) {
  cat("\nAnalysing ChIP-seq experiment", ID, "...\n")

  # Select counts of experiment
  df <- read.counts.per.ID[[ID]]

  if (!opt$outliers) {
    # Calculate lower and upper whiskers for all states (see ?ggplot2::geom_boxplot)
    ylims <- sapply(levels(df$state),
                    function(x) {
                      state.counts <- subset(df, state==x, select=counts)$counts
                      hinges <- quantile(state.counts, c(1,3)/4)
                      extremeWhiskers <- hinges + 1.5*diff(hinges)*c(-1,1)
                      return(c(min(state.counts[state.counts >= extremeWhiskers[1]]),
                               max(state.counts[state.counts <= extremeWhiskers[2]])))
                    })

    # Get extremal values of read counts
    ylim <- c(min(ylims[1,]), max(ylims[2,]))
  }

  # Create boxplot
  savename.box <- file.path(opt$output,
                            paste0("boxplot_chromstar-", basename(opt$chromstar),
                                   if (opt$raw) "_rawcounts-" else "_normalizedcounts-",
                                   ID,
                                   if (opt$outliers) "_outliers" else "",
                                   if (opt$nostats) "_nostats" else "",
                                   ".pdf"))

  cat("Saving boxplot in", savename.box, "...\n")
  if (file.exists(savename.box))
    warning(savename.box, " already exists and is overwritten.")

  ggplt.box <- ggplot(df, aes(state, counts)) +
    xlab("CNV state") +
    ylab(paste0(if(opt$raw) 'raw' else 'normalized',
                ' read counts (',
                unlist(strsplit(ID, '-'))[[1]], ')')) +
    theme_bw()

  # Show/hide outliers
  if (opt$outliers) {
    ggplt.box <- ggplt.box + geom_boxplot()
  } else {
    ggplt.box <- ggplt.box + geom_boxplot(outlier.colour = NA) +
      coord_cartesian(ylim = ylim + c(-0.01, 0.01) * diff(ylim))
  }

  # Add number of samples to plot
  if (!opt$nostats) {
    .give.n <- function(x){
      return(c(y = mean(quantile(x, probs=c(0.5, 0.75))), label = length(x)))
    }
    ggplt.box <- ggplt.box + stat_summary(fun.data = .give.n, size = 2, geom = "text")
  }

  ggsave(savename.box, plot=ggplt.box, width=10, height=10, limitsize=FALSE, units='cm')

  # Perform statistical tests
#  if (length(levels(df$state)) == 1) {
#    cat("Can not perform statistical testes since only one copy number state found.\n")
#  } else {
#    # Test normality with QQ plots
#    savename.qq <- file.path(opt$output,
#                             paste0("qqplot_chromstar-", basename(opt$chromstar),
#                                    if (opt$raw) "_rawcounts-" else "_normalizedcounts-",
#                                    ID, ".pdf"))
#
#    cat("Plotting quantile-quantile plots to", savename.qq, "...\n")
#
#    if (file.exists(savename.qq))
#      warning(savename.qq, "already exists and is overwritten.")
#
#    pdf(savename.qq, onefile=TRUE, width=8, height=8)
#    for (state in levels(df$state)) {
#      df.qq <- subset(df, state==state, select=counts)
#
#      # Calculations for quantile line (taken from https://gist.github.com/meren/4485081)
#      y <- quantile(df.qq$counts, c(0.25, 0.75))
#      x <- qnorm(c(0.25, 0.75))
#      slope <- diff(y)/diff(x)
#      int <- y[1L] - slope * x[1L]
#
#      ggplt.qq <- ggplot(data=df.qq, aes(sample=counts)) +
#        stat_qq() +
#        geom_abline(slope = slope, intercept = int) +
#        ggtitle(paste0("Normal Q-Q Plot of CNV state ", state)) +
#        theme_bw()
#
#      print(ggplt.qq)
#    }
#    invisible(dev.off())
#
#    if (length(levels(df$state)) == 2) {
#      # Perform t-test if only two copy number states
#      cat("Performing t-test since only two copy number states found ...\n")
#
#      tt <- t.test(counts ~ state, data=df)
#      print(tt)
#
#      ## Save results of t-test to file
#      savename.t <- file.path(opt$output,
#                              paste0("t-test_chromstar-", basename(opt$chromstar),
#                                     if (opt$raw) "_rawcounts-" else "_normalizedcounts-",
#                                     ID, ".txt"))
#
#      cat("Saving results of t-test to", savename.t, "...\n")
#
#      if (file.exists(savename.t))
#        warning(savename.t, " already exists and is overwritten.")
#
#      capture.output(tt, file=savename)
#    } else  {
#      # Perform ANOVA + post-hoc TukeyHSD if more than two different copy number states
#      cat("Performing ANOVA since more than two copy number states found ...\n")
#
#      model <- aov(data=df, counts ~ state)
#      aa <- anova(model)
#      print(aa)
#
#      # Save results of ANOVA to file
#      savename.anova <- file.path(opt$output,
#                                  paste0("anova_chromstar-", basename(opt$chromstar),
#                                         if (opt$raw) "_rawcounts-" else "_normalizedcounts-",
#                                         ID, ".txt"))
#
#      cat("Saving results of ANOVA to", savename.anova, "...\n")
#
#      if (file.exists(savename.anova))
#        warning(savename.anova, " already exists and is overwritten.")
#
#      capture.output(aa, file=savename.anova)
#
#      # Obtain p value of ANOVA
#      pvalue <- aa$"Pr(>F)"[1]
#
#      if (pvalue > 0.05) {
#        # Perform no post-hoc tests if no significant results in ANOVA
#        cat("Performing no post-hoc tests since ANOVA did not find any significant differences.\n")
#      } else {
#        # Perform post-hoc tests if significant results in ANOVA
#        cat("Performing post-hoc test since ANOVA did find significant differences ...\n")
#
#        tk <- TukeyHSD(model)
#        print(tk)
#
#        # Save results of TukeyHSD to file
#        savename.tukey <- file.path(opt$output,
#                                    paste0("tukey-hsd_chromstar-", basename(opt$chromstar),
#                                           if (opt$raw) "_rawcounts-" else "_normalizedcounts-",
#                                           ID, ".txt"))
#
#        cat("Saving results of TukeyHSD to", savename.tukey, "...\n")
#
#        if (file.exists(savename.tukey))
#          warning(savename.tukey, " already exists and is overwritten.")
#
#        capture.output(tk, file=savename.tukey)
#
#        # Plot results of TukeyHSD to file (http://stackoverflow.com/questions/33644034/how-to-visualize-pairwise-comparisons-with-ggplot2)
#        savename.tukeyplot <- file.path(opt$output,
#                                        paste0("plot_tukey-hsd_chromstar-", basename(opt$chromstar),
#                                               if (opt$raw) "_rawcounts-" else "_normalizedcounts-",
#                                               ID, ".pdf"))
#
#        cat("Plotting results of TukeyHSD to", savename.tukeyplot, "...")
#
#        if (file.exists(savename.tukeyplot))
#          warning(savename.tukeyplot, " already exists and is overwritten.")
#
#        tk.state = as.data.frame(tk$state)
#        tk.state$pair = rownames(tk.state)
#
#        print(tk.state)
#
#        ggplt.tk <- ggplot(tk.state, aes(colour=cut(`p adj`, c(-0.0001, 0.001, 0.01, 0.05, 1),
#                                                    label=c("p < 0.001", "p < 0.01","p < 0.05","not signficant")))) +
#geom_hline(yintercept=0, lty="11", colour="grey30") +
#geom_errorbar(aes(pair, ymin=lwr, ymax=upr), width=0.2) +
#geom_point(aes(pair, diff)) +
#labs(colour="") +
#ggtitle("Results of TukeyHSD") +
#theme_bw() +
#theme(axis.text.x=element_text(angle=90, hjust=1))
#
#        width <- max(15, 1.5*(1+nrow(tk.state))) + 1
#height <- max(3, sapply(tk.state$pair, nchar)) / 10 + 15
#
#ggsave(savename.tukeyplot, plot=ggplt.tk, width=width, height=height, limitsize=FALSE, units='cm')
#      }
#    }
#  }
}

cat("Done.\n")
