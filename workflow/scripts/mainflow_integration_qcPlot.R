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


library(tidyverse)
library(ggpmisc)
library(mclust)


source("workflow/scripts/mainflow_integration_plot_loader.R")
dots <- dots.full  # QC Plot needs raw dots.csv integrated data


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
           ggplot(aes(x=integratedIntensity))+
           geom_histogram(bins=30)+
           scale_y_continuous(expand = c(0.05,0.05))+
           scale_x_continuous(trans="log1p", expand=c(0,0))+
           xlab("Intensity (a.u.)")+ylab("Count")+
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
         message(paste("... component means=", fit.means[1], fit.means[2], sep="\t"))
         message(paste("... component sds=", fit.sds[1], fit.sds[2], sep="\t"))
         dots.no.overlap.full %>%
           ggplot(aes(x=integratedIntensity.log1p))+
           geom_density(bins=1500)+
           scale_y_continuous(expand = c(0.05,0.05))+
           scale_x_continuous(expand=c(0,0))+
           xlab("log1p(Intensity) (a.u.)")+ylab("Count (ALL dots in sample)")+
           stat_function(fun = dnorm, args = list(mean = fit.means[1], sd = fit.sds[1]), col="red", size=2)+
           stat_function(fun = dnorm, args = list(mean = fit.means[2], sd = fit.sds[2]), col="green", size=2)+
           custom.theme
       },
       stop("Unsupported QC plotting type"))


out.path <- snakemake@output[["plot"]]

ggsave(filename = basename(out.path),
       path = dirname(out.path),
       dpi = "retina",
       width = ifelse(plot.type != "volume", 5, 1) * plot.width.in,
       height = ifelse(plot.type != "volume", 5, 1) * plot.height.in,
       units = "in", limitsize = F)

warnings()