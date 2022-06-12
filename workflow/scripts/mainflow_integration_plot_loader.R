# This script contains helper functions and a data loading flow
# Will load summarized data from Stage I integration, together with user input to plot.csv, samples.csv, and plot.yaml

library(tidyverse)
library(yaml)

# ------ HELPER FUNCTION DEFINITIONS ------
verify.samples <- function(samples){
  # Input: a samples df
  valid <-
    all(c("sample", "image", "mask_path", "fishdot_path",
        "physicalSizeX", "physicalSizeY", "physicalSizeZ",
        "include", "num_cells") %in% names(samples))
  if (! valid){
    stop("Please modify samples.csv to include required columns")
  }
}

verify.plot.csv <- function(plot){
  # Input: a plot df
  valid <- all(c("sample", "group", "batch") %in% names(plot))
  if (! valid){
    stop("Please modify plot.csv to include required columns")
  }
}

verify.context <- function(){
  # Verifies whether snakemake object is present
  if (!exists("snakemake")){
  stop("Workflow script is only callable via snakemake.")
}
}

read.sum <- function(snakemake){
  # Reads in summarized data returning a list
  tryCatch({
    samples <- read_csv(snakemake@input[['samples']], show_col_types = F)
    plot <- read_csv(snakemake@input[['plot']], show_col_types = F)
    dots <- read_csv(snakemake@input[['dots']], show_col_types = F)
  }, error = function(e) stop("Could not read samples and/or plot and/or dots integration data."))
  list(samples = samples,
       plot = plot,
       dots = dots)
}

read.plot.config <- function(snakemake){

  # Reads in plot.yaml returning a list
  tryCatch({
    yaml.path <- file.path(dirname(snakemake@input[['samples']]), "plot.yaml")
    plot.config <- yaml.load_file(yaml.path)
  }, error = function(e) {
    stop("Please use plot.yaml to customize plotting. See documentation.")
  })

  # Processes plot parameter defaults
  plot.xlab <- plot.config$xlab
  plot.ylab <- plot.config$ylab
  plot.basefs <- plot.config$base_fs
  plot.xlab <- ifelse(is.null(plot.xlab), "Group", plot.xlab)
  plot.ylab <- ifelse(is.null(plot.ylab), "dots/cell", plot.ylab)
  plot.basefs <- ifelse(is.null(plot.basefs), 24, as.integer(plot.basefs))

  # Returns parsed config entries
  list(
    # If any of these are not given, NULL will be the value.
    plot.xlab = plot.xlab,
    plot.ylab = plot.ylab,
    plot.basefs = plot.basefs,
    group.ordered = plot.config$group_ordered,
    plot.ymin = plot.config$ymin,
    plot.ymax = plot.config$ymax,
  )
}

# ------ DATA LOADING ------

# Determine context
verify.context()


# Read in samples and plot integration data
#   adds: dots, plot, samples
data <- read.sum(snakemake)
list2env(data, .GlobalEnv)


# Read in optional plot.yaml and throw a message if not present
#   adds: plot config entries
plot.config <- read.plot.config(snakemake)
list2env(plot.config, .GlobalEnv)


# Check completeness of samples and plot data (all columns present?)
verify.samples(samples)
verify.plot.csv(plot)


# Get the plot
theme_set(theme_classic(base_size = plot.basefs))
dots %>%
  inner_join(plot, by = "sample") %>%
  inner_join(samples, by = c("sample", "image")) %>%
  group_by(group, image, sample, batch) %>%
  summarize(dot.count = n(), num.cells = min(num_cells), .groups = "drop") %>%
  mutate(group = ordered(group, levels = group.ordered)) -> dots

# Verify that groups are all defined by the YAML spec. Otherwise raise error.
if (any(is.na(dots$group)))
  stop("Please make sure that ALL groups are defined in plot.yaml.")
