"""
    Aggregates the following images together to form visualization of the results:
        Channel#       What                 Where
        1              marker raw data      singlechannel
        2              segmentation mask    segmentation
        3              fish raw data        singlechannel
        4              fishdot results      fishdot
"""

from os.path import join
from mainflow_helpers import *
import tifffile as tf  # image_base will not be that useful here
import numpy as np
from re import search
import csv
import os


def getChannelNameFromPath(path):
    """
        Extracts channel name from a segmentation or fishdot output path.
        :param path: path string to a segmentation or fishdot output
            This function looks for /*++ pattern and extracts the *
    """
    path = path[::-1]
    path = search(r"\+\+(.*?)/", path)
    return (path.groups()[0])[::-1]


def visualize(samples_csv_path, output_dir_base):
    """
        aggregates results together to form visualization images.
    :param samples_csv_path: path to the integration file samples.csv. Do NOT need user custom fields.
    :param output_dir_base: path to the output folder.
            Folders corresponding to different samples for this probe will be created.
            This function develops file paths based on samples.csv integration and snake_helper only.
            Assuming np.uint16 dtype
    :return: none.
    """
    with open(samples_csv_path, "r") as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            # Goal: get four paths to four single channel images.
            # Those images could be of different data type (uint8 or uint16, etc.)
            #   First, get mask and fishdot image paths
            mask_path = row['mask_path']
            fishdot_path = row['fishdot_path']
            fishdot_path = fishdot_path.replace(".loc4", "_spots.tif")
            #   Second, get channel names for getting the raw singlechannels
            segmentation_channel = getChannelNameFromPath(mask_path)
            fishdot_channel = getChannelNameFromPath(fishdot_path)
            #   Third, get base folder leading to the singlechannel
            sample = row['sample']
            image = row['image']
            probe_folder_path = search("(.*?)"+os.sep+sample, mask_path).groups()[0]
            sc_base_path = join(probe_folder_path, sample, "single_channel")
            marker_path = join(sc_base_path, segmentation_channel, image + ".tif")
            rawfish_path = join(sc_base_path, fishdot_channel, image + ".tif")
            #   Final, read in all images and put together and write
            #       Read in all images, each shape=zyx
            imgs = [tf.TiffFile(fp).asarray() for fp in [marker_path, mask_path, rawfish_path, fishdot_path]]
            img_shape = imgs[0].shape  # zyx for tifffile
            img_shape_z = img_shape[0]
            img_vis = np.empty(img_shape_z + (4,) + img_shape[1:], dtype=np.uint16)  # shape = zcyx
            for i in range(4):
                img_vis[:, i, :, :] = imgs[i]
            output_path = join(output_dir_base, sample)  # Output path part I - sample folder level
            os.makedirs(output_path, exist_ok=True)  # Make sure to create the sample folder
            output_path = join(output_path, row['image'] + ".tif")  # Full output path
            tf.imwrite(output_path, img_vis, photometric='minisblack')
            print_current_time("Visualization written for " + sample + "  " + image)
