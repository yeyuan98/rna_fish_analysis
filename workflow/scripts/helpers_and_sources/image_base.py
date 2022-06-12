"""
    Image basic IO capabilities
"""

import re
import numpy as np
from xml.etree import ElementTree as ET
import tifffile as tf


# OME-TIFF Metadata Processing

def get_namespace(element):
    """Get XML namespace string given an element.
    Use get_pixelsettings instead for accessing pixel attributes from TiffFile."""
    m = re.match('\{.*\}', element.tag)
    return m.group(0) if m else ''


def get_pixelsettings(tiff):
    """Get OMETIFF pixel attributes given an TiffFile.
     The original file must be OMETIFF."""
    # Input: TiffFile
    try:
        meta_root = ET.fromstring(tiff.ome_metadata)
        ome_schema = get_namespace(meta_root)
        ome_pixel = meta_root.find(f"{ome_schema}Image").find(f"{ome_schema}Pixels")
    except:
        raise ValueError("Tiff does not have OME metadata.")
    if ome_pixel is None:
        raise ValueError("OME-Tiff does not hold Pixel settings.")
    return ome_pixel.attrib


def get_z_direction(tiff):
    """
        Get Z-direction acquisition sequence of the image.
    :param tiff: TiffFile object of the image
    :return: either '+' or '-'.
        '+' means that the stack is taken such that absolute Z position is increasing with Z index
    """
    # Input: TiffFile
    try:
        meta_root = ET.fromstring(tiff.ome_metadata)
        ome_schema = get_namespace(meta_root)
        ome_pixel = meta_root.find(f"{ome_schema}Image").find(f"{ome_schema}Pixels")
    except:
        raise ValueError("Tiff does not have OME metadata.")
    # Get Z position values for the first and second planes, from OME-TIFF Plane metadata
    z_01 = [float(ome_pixel.find(f"{ome_schema}Plane[@TheZ='{i}']").attrib['PositionZ']) for i in range(2)]
    if z_01[0] < z_01[1]:
        return '+'
    else:
        return '-'


def get_channels(tiff, what="Name"):
    """
    Get channel information from the TiffFile object
    :param what: what entry to extract
    :param tiff: TiffFile
    :return: a list of channel information, only entry specified by what param.
    """
    try:
        meta_root = ET.fromstring(tiff.ome_metadata)
        ome_schema = get_namespace(meta_root)
        ome_channel = meta_root.find(f"{ome_schema}Image").find(f"{ome_schema}Pixels").findall(f"{ome_schema}Channel")
    except:
        raise ValueError("Tiff does not have OME metadata.")
    if ome_channel is None:
        raise ValueError("OME-Tiff does not hold Pixel settings.")
    return [a.attrib[what] for a in ome_channel]


# Image I/O and transformations


def load_image(path):
    """
    Couldn't be simpler equivalent to tf.TiffFile
    """
    return tf.TiffFile(path)


def load_channel(path, name):
    """
    Loads data from a given channel in the tiff file
    :param name: search string to find the channel
    :param path: path to the tif file
    :return: a numpy array, shape = zyx
    """
    img = load_image(path)
    img_chs = get_channels(img, what="Name")
    img_chs = [ch.find(name) for ch in img_chs]
    img_chs = [loc >= 0 for loc in img_chs]
    try:
        channel = img_chs.index(True)
    except:
        raise ValueError("Light train specified does not exist in file.")
    img = img.asarray()
    img = img[:, channel, :, :]
    return img


def save_channel(in_path, out_path, name):
    """
    Loads tiff image and extracts a given channel and saves the output.
    Also, will write shaped_metadata to the output tif file of the "PhysicalSize" entries from the input OME-TIFF.
        to access saved metadata, create TiffFile and use shaped_metadata[0] to access the dictionary.
        metadata keys & values will all be strings.
    :param in_path: Input ome.tif
    :param out_path: Output target path, will be raw tif wo/ metadata.
    :param name: search string to find the channel
    :return: none
    """
    img_ch = load_channel(in_path, name)  # Loads the specified channel
    img = tf.TiffFile(in_path)
    img = get_pixelsettings(img)
    img_PhysicalSize = dict()
    for (k, v) in img.items():
        if k.find("PhysicalSize") == 0:
            img_PhysicalSize[k] = v
    img_PhysicalSize['axes'] = 'ZYX'
    tf.imwrite(out_path, img_ch,
               metadata=img_PhysicalSize)  # Write image with PhysicalSize


def load_pixelSizes(path):
    """
        Reads an OME-TIFF file and fetch its physical sizes settings, returning a list for [X, Y, Z].
    :param path: path to an OME-TIFF file
    :return: a list of physical pixel sizes, in the order of [X, Y, Z]
    """
    f = load_image(path)
    pixels = get_pixelsettings(f)
    try:
        return [pixels["PhysicalSizeX"], pixels["PhysicalSizeY"], pixels["PhysicalSizeZ"]]
    except:
        raise ValueError("Target OME-TIFF file does not contain all PhysicalSizesXYZ info.")


def load_nonzero_count(path):
    """
        Reads a single channel ZYX TIFF file and return count of non-zero pixels.
    :param path: path to the zyx TIFF image
    :return: int of non-zero pixel count
    """
    f = load_image(path)
    return np.sum(f.asarray() != 0)


def load_z_direction(path):
    """
        Simple wrapper of get_z_direction to see whether larger Z-stack index means more positive stage Z location.
    :param path: path to the OME-TIFF image
    :return: '+' or '-', see get_z_direction
    """
    f = load_image(path)
    return get_z_direction(f)
