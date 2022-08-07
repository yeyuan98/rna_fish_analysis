if (F){
  "
    This script implements all QC plots except the batch faceted count plot.
    It requires Stage I integration outputs dots.csv, samples.csv.
    The following plots will be generated:
      1. per sample basis, intensity of the recognized & overlapped FISH dots
      2. per sample basis, residual of the recognized & overlapped FISH dots
      3. per sample basis, segmentation object volume normalized by num_cells
      4. per sample basis, working distance evaluation plots
  "
}


suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggpmisc))
suppressPackageStartupMessages(library(mclust))


source("workflow/scripts/mainflow_integration_plot_loader.R")
dots <- dots.full  # QC Plot needs raw dots.csv integrated data

# A custom function for error bar calculation. Note that this is copied from countPlot.R
# TODO: consider moving base.countPlot and sem.error to plot_loader.R
sem.error <- function(x){
  sem <- sd(x) / sqrt(length(x))
  m <- mean(x)
  data.frame(ymin=m-sem, y=m, ymax=m+sem)
}

# Check plotting type and add batch facet if QC plot is requested.
plot.type <- snakemake@wildcards[["plot_type"]]
switch(plot.type,
       volume={
         pixel3.volume <- with(dots, physicalSizeX * physicalSizeY * physicalSizeZ)
         dots %>%
           mutate(cell.volume = mask_pixel_volume * pixel3.volume / num_cells) %>%
           mutate(cell.equi.diameter = (cell.volume/pi*6)^(1/3)) -> dots
           dots %>%
             group_by(sample, image) %>%
             summarize(cell.equi.diameter = mean(cell.equi.diameter), .groups = "drop") %>%
             ggplot(aes(x=sample, y=cell.equi.diameter))+
             geom_point(alpha=0.3)+
             geom_jitter(alpha=0.3, width=0.25)+
             scale_y_continuous(expand = c(0.05,0.05))+
             xlab("Sample")+ylab("Object equivalent diameter (p.u.)")+  # p.u. = physical unit
             custom.theme
       },
       intensity={
         #  Intensity ~ sample
         dots %>%
           mutate(unique.label = paste(sample, image, sep = "\n")) %>%
           ggplot(aes(x=log1p(integratedIntensity)))+
           geom_histogram(bins=30)+
           scale_y_continuous(expand = c(0,0.01))+
           scale_x_continuous(expand=c(0,0))+
           xlab("log1p(Intensity) (a.u.)")+ylab("Count")+
           facet_wrap(vars(unique.label))+
           theme(strip.background = element_blank())+
           custom.theme
       },
       residual={
         #  Residual ~ sample
         dots %>%
           mutate(unique.label = paste(sample, image, sep = "\n")) %>%
           ggplot(aes(x=residuals))+
           geom_histogram(bins=30)+
           scale_y_continuous(expand = c(0.05,0.05))+
           scale_x_continuous(trans="log1p", expand=c(0,0))+
           xlab("Residuals (a.u.)")+ylab("Count")+
           facet_wrap(vars(unique.label))+
           theme(strip.background = element_blank())+
           custom.theme
       },
       workingDist={
         #  Intensity ~ Z, faceted by sample
         #    We align Z direction for each image, putting one end to be 0 with physical unit
         #    Simple linear model fitting is performed, and only the significant lines are plotted.
         dots %>%
           mutate(z.in.physical = # '+' direction increase with z_in_pix; otherwise decrease; both take the same range.
                    ifelse(z_direction == "+", (z_in_pix -1), (z_pixel_num - z_in_pix)) * physicalSizeZ) -> dots
         dots %>%
           mutate(unique.label = paste(sample, image, sep = "\n")) %>%
           ggplot(aes(x=z.in.physical, y=integratedIntensity))+
           geom_point(alpha=0.3)+
           scale_y_log10(expand = c(0.05, 0.05))+
           geom_smooth(method="lm", formula= y~x, se=F)+
           stat_fit_glance(method = 'lm', method.args = list(formula = y~x), geom = 'text_npc',
                       aes(label = paste0("P-value = ", signif(..p.value.., digits = 3),
                                          "\n","Adj. Rsq = ", signif(..adj.r.squared.., digits = 3)),
                           colour = ifelse(..p.value.. < 0.001, "#FF00FF","#000000")),
                       label.x = 1, label.y = 0.1, size = 6)+
           xlab("Z Position per dot (physical unit)")+
           ylab("Integrated intensity per dot (arbitrary unit)")+
           facet_wrap(vars(unique.label))+
           theme(strip.background = element_blank())+
           custom.theme
       },
       intensity_all={
         dots.no.overlap.full %>%
           mutate(integratedIntensity.log1p = log1p(integratedIntensity)) -> dots.no.overlap.full
         message(paste0("Performing GMM for ", snakemake@wildcards[["probe"]]))
         fit <- Mclust(dots.no.overlap.full$integratedIntensity.log1p,
                       G=2, modelNames = "V")
         fit.means <- fit$parameters$mean
         fit.sds <- sqrt(fit$parameters$variance$sigmasq)
         fit.pros <- fit$parameters$pro
         message(paste(snakemake@wildcards[["probe"]], "... component means=", fit.means[1], fit.means[2], sep="\t"))
         message(paste(snakemake@wildcards[["probe"]], "... component sds=", fit.sds[1], fit.sds[2], sep="\t"))
         message(paste(snakemake@wildcards[["probe"]], "... component proportions=",fit.pros[1], fit.pros[2], sep="\t"))
         #  We also save the GMM results to a file for replotting the results
         gmm.fit.path <- file.path(dirname(snakemake@output[["plot"]]), "gmm.fit.csv")
         write_csv(tibble(mean = fit.means, sd = fit.sds, proportions = fit.pros), file=gmm.fit.path)
         message(paste("GMM fit results written to", gmm.fit.path))
         dots.no.overlap.full %>%
           ggplot(aes(x=integratedIntensity.log1p))+
           geom_histogram(aes(y=..density..), bins=1500)+
           scale_y_continuous(expand = c(0,0.01))+
           scale_x_continuous(expand=c(0,0))+
           xlab("log1p(Intensity) (a.u.)")+ylab("Count (ALL dots in sample)")+
           stat_function(fun = function(x,...) dnorm(x,...) * fit.pros[1],
                         args = list(mean = fit.means[1], sd = fit.sds[1]), col="red", alpha=0.3, size=3)+
           stat_function(fun = function(x,...) dnorm(x,...) * fit.pros[2],
                         args = list(mean = fit.means[2], sd = fit.sds[2]), col="green", alpha = 0.3, size=3)+
           custom.theme
       },
       intensity_sample={
         # Plot dot intensity distribution for each sample. Each dot is a detected spot.
         dots %>%
           ggplot(aes(x=sample, y=log1p(integratedIntensity)))+
           geom_point()+
           stat_summary(fun.data = sem.error,
                        geom="errorbar", color="red", width=0.1)+  # This calculates s.e.m.
           stat_summary(fun = mean, geom="point", color="red", shape=3, size=3)+
           scale_y_continuous(expand = c(0,0.01))+
           xlab("Sample")+ylab("log1p(Intensity) (a.u.)")+
           custom.theme
       },
       stop("Unsupported QC plotting type"))


out.path <- snakemake@output[["plot"]]

plot.size.scale <- ifelse(plot.type %in% c("volume", "intensity_all"), 1, 5)

ggsave(filename = basename(out.path),
       path = dirname(out.path),
       dpi = "retina",
       #  Volume plot has the sample tags as X axis, give it wider plot.
       width = ifelse(plot.type == "volume" , 2, 1) * plot.size.scale * plot.width.in,
       height = plot.size.scale * plot.height.in,
       units = "in", limitsize = F)

warnings()