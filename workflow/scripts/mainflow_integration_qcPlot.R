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


source("workflow/scripts/mainflow_integration_plot_loader.R")


# Check plotting type and add batch facet if QC plot is requested.
plot.type <- snakemake@wildcards[["plot_type"]]
switch(plot.type,
       volume={
         pixel3.volume <- with(dots, physiscalSizeX * physicalSizeY * physicalSizeZ)
         dots %>%
           mutate(cell.volume = mask_pixel_volume * pixel3.volume / num_cells) %>%
           mutate(cell.equi.diameter = (cell.volume/pi*6)^(1/3)) %>%
           ggplot(aes(x=sample, y=cell.equi.diameter))+
           geom_point()+
           scale_y_continuous(expand = c(0.05,0.05))+
           xlab("Sample")+ylab("Equivalent diameter per object (physical unit)")
       },
       intensity={
         #  Intensity ~ sample
         dots %>%
           ggplot(aes(x=sample, y=integratedIntensity))+
           geom_point()+
           scale_y_log10(expand = c(0.05,0.05))+
           xlab("Sample")+ylab("Integrated intensity per dot (arbitrary unit)")
       },
       residual={
         #  Residual ~ sample
         dots %>%
           ggplot(aes(x=sample, y=residuals))+
           geom_point()+
           scale_y_log10(expand = c(0.05,0.05))+
           xlab("Sample")+ylab("Gaussian fit residuals per dot (arbitrary unit)")
       },
       workingDist={
         #  Intensity ~ Z, faceted by sample
         #    We align Z direction for each image, putting one end to be 0 with physical unit

       },
       stop("Unsupported QC plotting type"))


out.path <- snakemake@output[["plot"]]
ggsave(filename = basename(out.path),
       path = dirname(out.path),
       dpi = "retina")
