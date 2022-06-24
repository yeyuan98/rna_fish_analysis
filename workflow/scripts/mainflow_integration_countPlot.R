if (F){
  "
    This script implements count ~ group plotting integration.
    It requires Stage I integration outputs dots.csv, samples.csv.
    Two possible output files: merged.pdf and batch.qc.pdf
  "
}


library(tidyverse)
library(ggrepel)


source("workflow/scripts/mainflow_integration_plot_loader.R")

# A custom function for error bar calculation
sem.error <- function(x){
  sem <- sd(x) / sqrt(length(x))
  m <- mean(x)
  data.frame(ymin=m-sem, y=m, ymax=m+sem)
}

# Define the base plot without data
base.countPlot <-
  ggplot(data = NULL, aes(x = group, y = dots.per.cell))+
  geom_point()+
  stat_summary(fun.data = sem.error,
               geom="errorbar", color="red", width=0.1)+  # This calculates s.e.m.
  stat_summary(fun = mean, geom="point", color="red", shape=3, size=3)+  # This calculates mean.
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
         plot.data <- dots %>% mutate(dots.per.cell = dot.count / num.cells) %>% filter(include)
         plot.countPlot <- base.countPlot %+% plot.data
       },
       batch.qc={
         plot.data <- dots %>%
                        mutate(dots.per.cell = dot.count / num.cells) %>%
                        mutate(image = ifelse(include, "", image))
         plot.countPlot <- base.countPlot %+% plot.data
         plot.countPlot <- plot.countPlot +
                           facet_wrap(vars(batch)) +
                           geom_point(color=ifelse(plot.data$include, "black", "red"),
                                      label=image) +
                           geom_text_repel(max.overlaps = 30)
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
