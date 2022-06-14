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
           ggplot(aes(x=sample, y=integratedIntensity))+
           geom_point(alpha=0.3)+
           geom_jitter(alpha=0.3)+
           scale_y_log10(expand = c(0.05,0.05))+
           xlab("Sample")+ylab("Integrated intensity per dot (a.u.)")+
           custom.theme
       },
       residual={
         #  Residual ~ sample
         dots %>%
           ggplot(aes(x=sample, y=residuals))+
           geom_point(alpha=0.3)+
           geom_jitter(alpha=0.3)+
           scale_y_log10(expand = c(0.05,0.05))+
           xlab("Sample")+ylab("Gaussian fit residuals per dot (a.u.)")+
           custom.theme
       },
       workingDist={
         #  Intensity ~ Z, faceted by sample
         #    We align Z direction for each image, putting one end to be 0 with physical unit
         dots %>%
           mutate(z.in.physical = # '+' direction increase with z_in_pix; otherwise decrease; both take the same range.
                    ifelse(z_direction == "+", (z_in_pix -1), (z_pixel_num - z_in_pix)) * physicalSizeZ) %>%
           ggplot(aes(x=z.in.physical, y=integratedIntensity, group=image, color=image))+
           #geom_point(alpha=0.3, show.legend = F)+
           geom_smooth(method="lm", formula= y~x, se=T, show.legend = F)+
           scale_y_log10(expand = c(0.05, 0.05))+
           xlab("Z Position per dot (physical unit)")+
           ylab("Integrated intensity per dot (arbitrary unit)")+
           facet_wrap(vars(sample))+
           custom.theme
       },
       stop("Unsupported QC plotting type"))


out.path <- snakemake@output[["plot"]]
ggsave(filename = basename(out.path),
       path = dirname(out.path),
       dpi = "retina",
       width = plot.width.in,
       height = plot.height.in,
       units = "in")

warnings()