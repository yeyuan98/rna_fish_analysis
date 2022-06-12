if (F){
  "
    This script implements count ~ group plotting integration.
    It requires Stage I integration outputs dots.csv, samples.csv.
    Two possible output files: merged.pdf and batch.qc.pdf
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


source("workflow/scripts/mainflow_integration_plot_loader.R")


dots %>%
  mutate(dots.per.cell = dot.count / num.cells) %>%
  ggplot(aes(x = group, y = dots.per.cell))+
  geom_point()+
  scale_y_continuous(expand = c(0.05,0.05), limits = c(plot.ymin, plot.ymax))+
  xlab(plot.xlab)+
  ylab(plot.ylab) -> plot.countPlot

# Check plotting type and add batch facet if QC plot is requested.
plot.type <- snakemake@wildcards[["plot_type"]]
switch(plot.type,
       merged=message("Generating merged count plot"),
       batch.qc={
         plot.countPlot <- plot.countPlot + facet_wrap(vars(batch))
         message("Generating batch faceted plot")
       },
       stop("Unsupported plotting type"))

out.path <- snakemake@output[["plot"]]
ggsave(filename = basename(out.path),
       path = dirname(out.path),
       dpi = "retina")
