"""
    Functions for plotting a list of Gaussian dots
"""

import numpy as np


def gaussian_gen(sigma_xy, sigma_z):
    """
        Generates unit Gaussian numerical values given sigmas.
        Will take +- 3-times sigma.
    :param sigma_xy: Gaussian function. Unit: pixel.
    :param sigma_z: Gaussian function. Unit: pixel.
    """
    xy_dim_half = np.ceil(sigma_xy * 3).astype(np.int16)
    z_dim_half = np.ceil(sigma_z * 3).astype(np.int16)
    x = np.linspace(-xy_dim_half, xy_dim_half, 2 * xy_dim_half + 1)
    y = np.linspace(-xy_dim_half, xy_dim_half, 2 * xy_dim_half + 1)
    z = np.linspace(-z_dim_half, z_dim_half, 2 * z_dim_half + 1)
    xx, yy, zz = np.meshgrid(x, y, z)
    gaussian = np.exp(-xx**2 / (2*sigma_xy**2))
    gaussian *= np.exp(-yy**2 / (2*sigma_xy**2))
    gaussian *= np.exp(-zz**2 / (2*sigma_z**2))
    gaussian /= np.sum(gaussian)
    return gaussian


def bounding_box(x, y, z, sigma_xy, sigma_z, pixel_sizes):
    """
        Get bounding box for Gaussian value assignment, 3-sigma
    :param x: x pixel location, zero-based
    :param y: y pixel location, zero-based
    :param z: z pixel location, zero-based
    :param sigma_xy: Gaussian function. Unit: pixel
    :param sigma_z: Gaussian function. Unit: pixel
    :param pixel_sizes: A list of x, y, z dimensions of the image
    :return: dict with keys `image_loc`, `base_loc`.
        both locs are list [[x_min, x_max], [y_min, y_max], [z_min, z_max]]
        the min and max values are BOTH inclusive
    """
    sigma_xy_3 = np.ceil(sigma_xy * 3).astype(np.int16)
    sigma_z_3 = np.ceil(sigma_z * 3).astype(np.int16)
    image_loc = [[max(0, x - sigma_xy_3), min(pixel_sizes[0] - 1, x + sigma_xy_3)],
                 [max(0, y - sigma_xy_3), min(pixel_sizes[1] - 1, y + sigma_xy_3)],
                 [max(0, z - sigma_z_3), min(pixel_sizes[2] - 1, z + sigma_z_3)]]
    base_loc = [[max(0, sigma_xy_3 - x), min(2*sigma_xy_3, sigma_xy_3 + pixel_sizes[0] - 1 - x)],
                [max(0, sigma_xy_3 - y), min(2*sigma_xy_3, sigma_xy_3 + pixel_sizes[1] - 1 - y)],
                [max(0, sigma_z_3 - z), min(2*sigma_z_3, sigma_z_3 + pixel_sizes[2] - 1 - z)]]
    return {"image_loc": image_loc,
            "base_loc": base_loc}


def gaussian_plot(dots_dict, pixel_sizes, sigma_xy, sigma_z):
    """
        Plots a list of Gaussian dots
    :param dots_dict: A dict of dots to plot. Should have four keys `x`, `y`, `z`, `integratedIntensity`
        those should follow ImageJ dimension of xyz. Unit: pixel.
        assuming that all dict entry lists have the same length.
        x, y, z should be zero-based.
    :param pixel_sizes: A list of x, y, z dimensions of the image
        those should follow ImageJ dimension of xyz. Unit: pixel.
    :param sigma_xy: Gaussian function. Unit: pixel.
    :param sigma_z: Gaussian function. Unit: pixel.
    :return: An image of shape xyz
    """
    img = np.zeros(shape=pixel_sizes, dtype=np.uint16)
    xs, ys, zs = dots_dict['x'], dots_dict['y'], dots_dict['z']
    int_intensity = dots_dict['integratedIntensity']
    num_dots = len(xs)
    gaussian_base = gaussian_gen(sigma_xy, sigma_z)
    for i in range(num_dots):
        x, y, z = xs[i], ys[i], zs[i]
        b_box = bounding_box(x, y, z, sigma_xy, sigma_z, pixel_sizes)
        img_x = b_box["image_loc"][0]
        img_y = b_box["image_loc"][1]
        img_z = b_box["image_loc"][2]
        base_x = b_box["base_loc"][0]
        base_y = b_box["base_loc"][1]
        base_z = b_box["base_loc"][2]
        img[img_x[0]:(img_x[1]+1),
            img_y[0]:(img_y[1]+1),
            img_z[0]:(img_z[1]+1)] = gaussian_base[base_x[0]:(base_x[1]+1),
                                                   base_y[0]:(base_y[1]+1),
                                                   base_z[0]:(base_z[1]+1)] * int_intensity[i]
    return img
