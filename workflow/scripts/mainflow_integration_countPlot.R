if (F){
  "
    This script implements count ~ group plotting integration.
    It requires Stage I integration outputs dots.csv, samples.csv.
    Two possible output files: merged.pdf and batch.qc.pdf
  "
}


library(tidyverse)


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
