if (F){
  "
    This script implements count ~ group plotting integration.
    It requires Stage I integration outputs dots.csv, samples.csv.
  "
}
if (T){
  instruction.samples <- "
    The user is required to mofidy the samples.csv file to include the following image-specific information:
      column    values      description
      include   T or F      should R include this image for plotting
      num_cells integer     how many cells are there in this image
  "
  instruction.plot.csv <- "
    The user is also required to modify the plot.csv file to include the following sample-specific information:
      column    values      description
      group     character   what x-label should the sample take (will merge samples with the same group)
      batch     character   what batch is the sample (allows batch splitting of samples from the same group, QC only)
  "
  instruction.plot.yaml <- "
    If plot.yaml exists in the same folder as samples.csv, the script will use it for plot customization:
      --- only put in top level entries ---
      entry       ggplot2_param     default
      xlab        xlab()            'group'
      ylab        ylab()            'counts/cell'
  "
}

library(tidyverse)
library(yaml)

# Determine context
tryCatch(snakemake, error = function(e) stop("Workflow script is only callable via snakemake."))

# Read in samples and plot integration data
tryCatch({
  samples <- read_csv(snakemake@input[['samples']])
  plot <- read_csv(snakemake@input[['plot']])
}, error = function(e) stop("Could not read samples and/or plot integration data."))

# Read in optional plot.yaml and throw a message if not present
tryCatch({
  yaml.path <- file.path(dirname(snakemake@input[['samples']]), "plot.yaml")
  plot.config <- yaml.load_file(yaml.path)
}, error = function(e) {message(paste("To further configure your plots:", instruction.plot.yaml, sep = "\n"))})

# Check completeness of samples data (all columns present?)

# Check completeness of plot data (all columns present?)

# Get the plot