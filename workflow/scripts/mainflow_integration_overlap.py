"""
    Workflow script to get fish dots in the given segmentation mask
"""
from helpers_and_sources.mainflow_helpers import *
from helpers_and_sources import gaussian_plot
import tifffile as tf
import numpy as np
import csv
from os.path import split, splitext, join
from os import makedirs


def integration_overlap(samples_csv, output_paths):
    """
        For each image file in samples_csv, reads in segmentation mask and dots list, and performs union.
    :param samples_csv: Path to the integrated per image information
    :param output_paths: Paths to the output dots csv file and visualization folder
    """
    result = []  # will be a list of dicts for csv.DictWriter.
    result_complete = []  # All dots without any overlap check
    dots_csv = output_paths["dots"]
    dots_csv_complete = output_paths["dots_complete"]
    viz_dir = output_paths["gaussian_fit"]
    makedirs(viz_dir, exist_ok=True)  # visualization folder
    with open(samples_csv, 'r') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            # Process each image (part of a sample)
            print_current_time("Integration overlap: processing image= " + row['mask_path'])
            sample = row['sample']
            mask = tf.TiffFile(row['mask_path']).asarray()  # shape=zyx
            mask = np.transpose(mask, axes=(2, 1, 0))  # shape=xyz
            psx, psy, psz = float(row['physicalSizeX']), float(row['physicalSizeY']), float(row['physicalSizeZ'])
            image_name = splitext(split(row['mask_path'])[1])[0]
            viz_dots_dict = {"x": [], "y": [], "z": [], "integratedIntensity": []}  # for gaussian_plot
            makedirs(join(viz_dir, sample), exist_ok=True)  # Make directory for the sample if needed
            viz_img_path = join(viz_dir, sample, image_name, ".tif")  # Save path for the visualization
            with open(row['fishdot_path'], 'r') as dotfile:
                dot_reader = csv.DictReader(dotfile, delimiter='\t')
                for dotrow in dot_reader:
                    # Process each dot
                    # WARNING ------ MATLAB's pixel order is flipped. Images with ImageJ shape (x, y) will be (y, x)!!!
                    # MATLAB uses (row, column) order, and row is basically y for other image processing routines.
                    x, y, z = float(dotrow['y_in_pix']), float(dotrow['x_in_pix']), float(dotrow['z_in_pix'])
                    x, y, z = round(x) - 1, round(y) - 1, round(z) - 1  # to shift 1px for 0-based array index
                    x, y, z = [max(t, 0) for t in (x, y, z)]  # make sure no negatives
                    dot_dict = dotrow.copy()
                    dot_dict.update({'sample': sample, 'image': image_name,
                                     'physicalSizeX': psx, 'physicalSizeY': psy, 'physicalSizeZ': psz})
                    if mask[x, y, z] > 0:  # dot in segmentation mask
                        result.append(dot_dict)
                        viz_dots_dict["x"].append(x)
                        viz_dots_dict["y"].append(y)
                        viz_dots_dict["z"].append(z)
                        viz_dots_dict["integratedIntensity"].append(round(dotrow["integratedIntensity"]))
                    result_complete.append(dot_dict)  # add to the complete list regardless
            #  Save visualization for current image
            sigma_xy = round(snakemake.config["fishdot"]["psf_xy"] / (1E3 * psx))
            sigma_z = round(snakemake.config["fishdot"]["psf_z"] / (1E3 * psz))
            viz_array = gaussian_plot.gaussian_plot(viz_dots_dict, mask.shape, sigma_xy, sigma_z)
            tf.imwrite(viz_img_path, np.transpose(viz_array, axes=(2, 1, 0)),
                       compression=tf.TIFF.COMPRESSION.DEFLATE, imagej=True,
                       metadata={'axes': 'ZYX'})
    # Next, write results
    with open(dots_csv, 'w') as csvfile:
        result_fields = result[0].keys()
        writer = csv.DictWriter(csvfile, fieldnames=result_fields)
        writer.writeheader()
        writer.writerows(result)
    with open(dots_csv_complete, 'w') as csvfile:
        result_fields = result_complete[0].keys()
        writer = csv.DictWriter(csvfile, fieldnames=result_fields)
        writer.writeheader()
        writer.writerows(result_complete)
        print_current_time("Integration overlap: done.")

try:
    snakemake
except NameError:
    raise ReferenceError("Mainflow integration is only compatible with snakemake script directive.")
integration_overlap(snakemake.input[0], snakemake.output)
