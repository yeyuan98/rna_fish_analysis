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
      entry         ggplot2_param     default
      xlab          xlab()            'group'
      ylab          ylab()            'counts/cell'
      base_fs       theme(base_size)  24
      group_ordered NA                NA
      ymin          scale_y           NA
      ymax          scale_y           NA
    group_ordered should be a YAML sequence (dash lines) representing plotting order (left->right).
  "
}

library(tidyverse)
library(yaml)
source("workflow/scripts/mainflow_integration_countPlot_helper.R")

# Determine context
if (!exists("snakemake")){
  stop("Workflow script is only callable via snakemake.")
}

# Read in samples and plot integration data
tryCatch({
  samples <- read_csv(snakemake@input[['samples']], show_col_types = F)
  plot <- read_csv(snakemake@input[['plot']], show_col_types = F)
  dots <- read_csv(snakemake@input[['dots']], show_col_types = F)

  # Special - group parameter must be integer for plotting x-axis, and we check this
  tryCatch(plot$group <- as.integer(plot$group),
          warning = function(e) stop("Group must be integer for plotting"))
}, error = function(e) stop("Could not read samples and/or plot and/or dots integration data."))

# Read in optional plot.yaml and throw a message if not present
tryCatch({
  yaml.path <- file.path(dirname(snakemake@input[['samples']]), "plot.yaml")
  plot.config <- yaml.load_file(yaml.path)
  plot.xlab <- plot.config$xlab
  plot.ylab <- plot.config$ylab  # If any of these are not given, NULL will be the value.
  plot.basefs <- plot.config$base_fs
  group.ordered <- plot.config$group_ordered
  plot.ymin <- plot.config$ymin
  plot.ymax <- plot.config$ymax
}, error = function(e) {stop("Please use plot.yaml to customize plotting. See documentation.")})

# Check completeness of samples data (all columns present?)
verify.samples(samples)

# Check completeness of plot data (all columns present?)
verify.plot.csv(plot)

# Process plot parameter defaults
plot.xlab <- ifelse(is.null(plot.xlab), "Group", plot.xlab)
plot.ylab <- ifelse(is.null(plot.ylab), "dots/cell", plot.ylab)
plot.basefs <- ifelse(is.null(plot.basefs), 24, as.integer(plot.basefs))

# Get the plot
theme_set(theme_classic(base_size = plot.basefs))
dots %>%
  inner_join(plot, by = "sample") %>%
  inner_join(samples, by = c("sample", "image")) %>%
  group_by(group, image, sample) %>%
  summarize(dot.count = n(), num.cells = min(num_cells), .groups = "drop") %>%
  mutate(group = ordered(group, levels = group.ordered)) -> dots

# Verify that groups are all defined by the YAML spec. Otherwise raise error.
if (any(is.na(dots$group)))
  stop("Please make sure that ALL groups are defined in plot.yaml.")

dots %>%
  mutate(dots.per.cell = dot.count / num.cells) %>%
  ggplot(aes(x = group, y = dots.per.cell))+
  geom_point()+
  scale_y_continuous(expand = c(0,0), limits = c(plot.ymin, plot.ymax))+
  xlab(plot.xlab)+
  ylab(plot.ylab) -> plot.countPlot

out.path <- snakemake@output[["plot"]]
ggsave(filename = basename(out.path),
       path = dirname(out.path),
       dpi = "retina")
