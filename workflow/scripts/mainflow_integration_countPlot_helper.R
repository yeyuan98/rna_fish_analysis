# Helper functions for countPlot workflow.

verify.samples <- function(samples){
  # Input: a samples df
  valid <-
    all(c("sample", "image", "mask_path", "fishdot_path",
        "physicalSizeX", "physicalSizeY", "physicalSizeZ",
        "include", "num_cells") %in% names(samples))
  if (! valid){
    stop("Please modify samples.csv to include required columns")
  }
}

verify.plot.csv <- function(plot){
  # Input: a plot df
  valid <- all(c("sample", "group", "batch") %in% names(plot))
  if (! valid){
    stop("Please modify plot.csv to include required columns")
  }
}
