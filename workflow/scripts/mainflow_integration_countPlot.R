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
  mutate(dots.per.cell = dot.count / num.cells) -> dots
dots %>%
  ggplot(aes(x = group, y = dots.per.cell))+
  geom_point()+
  scale_y_continuous(expand = c(0.05,0.05), limits = c(plot.ymin, plot.ymax))+
  xlab(plot.xlab)+
  ylab(plot.ylab)+
  custom.theme -> plot.countPlot


# Check plotting type and add batch facet if QC plot is requested.
plot.type <- snakemake@wildcards[["plot_type"]]
switch(plot.type,
       merged=message("Generating merged count plot"),
       batch.qc={
         plot.countPlot <- plot.countPlot + facet_wrap(vars(batch))
         message("Generating batch faceted plot")
       },
       replot={
         gmm.path <- file.path(dirname(snakemake@input[["dots"]]), "qcPlots", "gmm.fit.csv")
         if (!exists(gmm.path)) stop("Replotting requires GMM fit data. Please performed int_qc rule first.")
         gmm.fit <- read_csv(gmm.path)
         intensity.threshold <- gmm.fit$mean[1]
         plot.countPlot <- plot.countPlot %+% subset(dots, integratedIntensity >= intensity.threshold)
       },
       stop("Unsupported plotting type"))


out.path <- snakemake@output[["plot"]]

width.scale.factor <- ifelse(plot.type == "batch.qc", 2, 1)

ggsave(filename = basename(out.path),
       path = dirname(out.path),
       dpi = "retina",
       width = width.scale.factor * plot.width.in,
       height = plot.height.in,
       units = "in")
