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


# The following functions decide what constitutes "overlap" given a dot (x,y,z) and the mask array (xyz)
def overlap_simple(x, y, z, mask):
    """
        Simple overlap check.
        Returns True if the dot is in the mask.
    :param x: x coordinate of the dot
    :param y: y coordinate of the dot
    :param z: z coordinate of the dot
    :param mask: mask array
    :return: True if the dot is in the mask
    """
    return mask[x, y, z] > 0


def overlap_manual(x, y, z, mask, physicalSizes, threshold):
    """
        Overlap check for the manual mask mode.
        Returns True if the dot is in vicinity of the mask.
        Exact distance is hard to compute and this function adopts a rough estimate.
    :param x: x coordinate of the dot
    :param y: y coordinate of the dot
    :param z: z coordinate of the dot
    :param mask: mask array
    :param physicalSizes: physical sizes of the image (x, y, z)
    :param threshold: physical distance threshold for overlap (unit same as physicalSizes, on all axes independently)
    :return: True if the dot is in the mask
    """
    # First, get the dot and physicalSizes
    center = np.array([x, y, z])  # shape will be (3,)
    physicalSizes = np.array(physicalSizes)  # shape will be (3,)
    # Next, get bounding box given its center and physicalSizes
    box_lower = center - threshold / physicalSizes / 2
    box_lower = np.maximum(np.int32(box_lower), 0)
    box_upper = center + threshold / physicalSizes / 2
    box_upper = np.minimum(np.int32(box_upper), np.array(mask.shape) - 1)
    # Finally, see if any of the mask is in the box
    return np.any(mask[box_lower[0]:box_upper[0], box_lower[1]:box_upper[1], box_lower[2]:box_upper[2]])


def integration_overlap(samples_csv, output_paths):
    """
        For each image file in samples_csv, reads in segmentation mask and dots list, and performs union.
        Supports manually annotated masks. Refer to the pipeline documentation on the manual overlap mode.
    :param samples_csv: Path to the integrated per image information
    :param output_paths: Paths to the output dots csv file and visualization folder
    """
    result = []  # will be a list of dicts for csv.DictWriter.
    result_complete = []  # All dots without any overlap check
    dots_csv = output_paths["dots"]
    dots_csv_complete = output_paths["dots_complete"]
    viz_dir = output_paths["gaussian_fit"]
    makedirs(viz_dir, exist_ok=True)  # visualization folder
    manual_mode = snakemake.config["segmentation"]["manual"]
    with open(samples_csv, 'r') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            # Process each image (part of a sample)
            if VERBOSE:
                print_current_time("Integration overlap: processing image= " + row['mask_path'])
            sample = row['sample']
            mask = tf.TiffFile(row['mask_path']).asarray()  # shape=zyx
            mask = np.transpose(mask, axes=(2, 1, 0))  # shape=xyz
            psx, psy, psz = float(row['physicalSizeX']), float(row['physicalSizeY']), float(row['physicalSizeZ'])
            image_name = splitext(split(row['mask_path'])[1])[0]
            viz_dots_dict = {"x": [], "y": [], "z": [], "integratedIntensity": []}  # for gaussian_plot
            makedirs(join(viz_dir, sample), exist_ok=True)  # Make directory for the sample if needed
            viz_img_path = join(viz_dir, sample, image_name + ".tif")  # Save path for the visualization
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
                    if manual_mode:
                        # TODO: physical unit may not be micrometer and current code may break for different setups.
                        dot_overlapped = overlap_manual(x, y, z, mask, [psx, psy, psz],
                                                        snakemake.config["segmentation"]["manual_threshold"] / 1E3)

                    else:
                        dot_overlapped = overlap_simple(x, y, z, mask)
                    if dot_overlapped:  # dot in segmentation mask
                        result.append(dot_dict)
                        viz_dots_dict["x"].append(x)
                        viz_dots_dict["y"].append(y)
                        viz_dots_dict["z"].append(z)
                        viz_dots_dict["integratedIntensity"].append(round(float(dotrow["integratedIntensity"])))
                    result_complete.append(dot_dict)  # add to the complete list regardless
            #  Save visualization for current image
            sigma_xy = round(snakemake.config["fishdot"]["psf_xy"] / (1E3 * psx))
            sigma_z = round(snakemake.config["fishdot"]["psf_z"] / (1E3 * psz))
            viz_array = gaussian_plot.gaussian_plot(viz_dots_dict, mask.shape, sigma_xy, sigma_z)
            tf.imwrite(viz_img_path, np.transpose(viz_array, axes=(2, 1, 0)),
                       compression=tf.TIFF.COMPRESSION.DEFLATE, imagej=True,
                       metadata={'axes': 'ZYX'})
    # Next, write results
    with open(dots_csv, 'w', newline='') as csvfile:
        result_fields = result[0].keys()
        writer = csv.DictWriter(csvfile, fieldnames=result_fields)
        writer.writeheader()
        writer.writerows(result)
    with open(dots_csv_complete, 'w', newline='') as csvfile:
        result_fields = result_complete[0].keys()
        writer = csv.DictWriter(csvfile, fieldnames=result_fields)
        writer.writeheader()
        writer.writerows(result_complete)
        if VERBOSE:
            print_current_time("Integration overlap: done.")

try:
    snakemake
except NameError:
    raise ReferenceError("Mainflow integration is only compatible with snakemake script directive.")
VERBOSE = snakemake.config["resources"]["verbose_log"]["integration"]
integration_overlap(snakemake.input[0], snakemake.output)
