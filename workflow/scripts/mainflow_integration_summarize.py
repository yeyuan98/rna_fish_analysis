"""
    This script fetches the following information for each image of a given probe:
        mask path, fishdot path, sample name, pixel size x, pixel size y, pixel size z
"""

from helpers_and_sources.integration_visualize import visualize
from helpers_and_sources import image_base
from snake_helper import *
from helpers_and_sources.mainflow_helpers import *
import csv


def query_images_by_converted(sample_path, config):
    """
        Get a list of images (without extension) by looking at sample_convert workflow output
    :param sample_path: path to the sample result folder
    :param config: snakemake config
    :return: a list of image names (NOT path) without any extension
    """
    converted_path = join(sample_path, output_workflow("convert", config))  # get path to converted image dir
    images_converted = listdir(converted_path)  # listdir
    images_converted = [f for f in images_converted if not f.startswith(".")]  # remove hidden files
    images_converted = [f for f in images_converted if (f.find(".ome.tif") >= 0)]  # look for .ome.tif files
    images_converted = [f.replace(".ome.tif", "") for f in images_converted]  # remove extension .ome.tif
    return images_converted


def summarize(sample_paths, output_file_paths, config):
    """
        Generates summary entry for each image file of a given probe.
    :param sample_paths: list of paths to different sample dirs
    :param output_file_paths: named paths to the specified output files
        samples is for per image summary and plot is for per sample configuration of plotting
    :param config: snakemake config
    """
    # For each sample directory path
    #   get image file names WITHOUT extension
    #   For each image file
    #       get each data of interest, append to result lists
    # Write to csv with core csv module
    output_samples_path = output_file_paths['samples']
    output_plot_path = output_file_paths['plot']
    samples_header = ["sample", "image", "mask_path", "fishdot_path", "physicalSizeX", "physicalSizeY", "physicalSizeZ"]
    samples_rows = []
    plot_header = ["sample"]
    plot_rows = []
    for sample_path in sample_paths:
        image_names = query_images_by_converted(sample_path, config)
        sample_name = split(sample_path)[1]
        plot_rows.append([sample_name])
        print_current_time("Processing sample: " + sample_name)
        for image_name in image_names:
            # columns:
            mask_path = join(sample_path, output_workflow("segmentation", config), image_name + ".tif")
            fishdot_path = join(sample_path, output_workflow("fishdot", config), image_name + ".loc4")
            converted_path = join(sample_path, output_workflow("convert", config), image_name + ".ome.tif")
            print_current_time("      converted image path: " + converted_path)
            pixel_sizes = image_base.load_pixelSizes(converted_path)
            samples_rows.append([sample_name, image_name,
                                 mask_path, fishdot_path, pixel_sizes[0], pixel_sizes[1], pixel_sizes[2]])
    with open(output_samples_path, 'w') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(samples_header)
        writer.writerows(samples_rows)
    with open(output_plot_path, 'w') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(plot_header)
        writer.writerows(plot_rows)
    print_current_time("Done summary. Aggregating images for visualization...")
    output_vis_path = output_file_paths['visualization']
    visualize(output_samples_path, output_vis_path)


try:
    snakemake
except NameError:
    raise ReferenceError("Mainflow integration_summarize is only compatible with snakemake script directive.")
summarize(snakemake.input, snakemake.output, snakemake.config)
