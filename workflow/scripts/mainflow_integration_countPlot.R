if (F){
  "
    This script implements count ~ group plotting integration.
    It requires Stage I integration outputs dots.csv, samples.csv.
    Two possible output files: merged.pdf and batch.qc.pdf
  "
}


library(tidyverse)


source("workflow/scripts/mainflow_integration_plot_loader.R")


# Define the base plot without data
base.countPlot <-
  ggplot(data = NULL, aes(x = group, y = dots.per.cell))+
  geom_point()+
  stat_summary(fun.data = "mean_cl_normal", fun.args = list(mult=1),
               geom="errorbar", color="red", width=0.2)+  # This calculates s.e.m.
  stat_summary(fun = mean, geom="point", color="red", size=2, alpha=0.5)+  # This calculates mean.
  scale_y_continuous(expand = c(0.05,0.05), limits = c(plot.ymin, plot.ymax))+
  xlab(plot.xlab)+
  ylab(plot.ylab)+
  custom.theme


# Check plotting type and add batch facet if QC plot is requested.
plot.type <- snakemake@wildcards[["plot_type"]]
probe <- snakemake@wildcards[["probe"]]
switch(plot.type,
       merged={
         message(paste("Generating merged count plot for", probe))
         plot.data <- dots %>% mutate(dots.per.cell = dot.count / num.cells)

         print(plot.data)

         plot.countPlot <- base.countPlot %+% plot.data
       },
       batch.qc={
         plot.data <- dots %>% mutate(dots.per.cell = dot.count / num.cells)
         plot.countPlot <- base.countPlot %+% plot.data
         plot.countPlot <- plot.countPlot + facet_wrap(vars(batch))
         message(paste("Generating batch faceted count plot for", probe))
       },
       replot={
         # read in GMM fit
         gmm.path <- file.path(dirname(snakemake@input[["dots"]]), "qcPlots", "gmm.fit.csv")
         if (!file.exists(gmm.path)) stop("Replotting requires GMM fit data. Please performed int_qc rule first.")
         gmm.fit <- read_csv(gmm.path)
         intensity.threshold <- exp(gmm.fit$mean[1]) - 1
         # filter data
         dots.full %>%
           filter(integratedIntensity >= intensity.threshold) -> dots.filtered
         message(paste("Filtered dots by intensity for ", probe, "before filter nrows=", nrow(dots.full)))
         print(paste("Filtered dots by intensity for ", probe, "after filter nrows=", nrow(dot.filtered)))
         dots.filtered %>%
           group_by(group, image, sample, batch) %>%
           summarize(dot.count = n(), num.cells = min(num_cells), .groups = "drop") -> dots.filtered
         # add on data to plot
         plot.countPlot <- base.countPlot %+% (dots.filtered %>% mutate(dots.per.cell = dot.count / num.cells))
         message(paste("Generating filtered count plot for", probe))
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
