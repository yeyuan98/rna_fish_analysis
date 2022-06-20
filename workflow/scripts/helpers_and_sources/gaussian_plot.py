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


def gaussian_plot(dots_dict, pixel_sizes, sigma_xy, sigma_z):
    """
        Plots a list of Gaussian dots
    :param dots_dict: A dict of dots to plot. Should have four keys `x`, `y`, `z`, `integratedIntensity`
        those should follow ImageJ dimension of xyz. Unit: pixel.
        assuming that all dict entry lists have the same length.
    :param pixel_sizes: A list of x, y, z dimensions of the image
        those should follow ImageJ dimension of xyz. Unit: pixel.
    :param sigma_xy: Gaussian function. Unit: pixel.
    :param sigma_z: Gaussian function. Unit: pixel.
    """
    img = np.zeros(shape=pixel_sizes, dtype=np.uint16)
    xs, ys, zs = dots_dict['x'], dots_dict['y'], dots_dict['z']
    int_intensity = dots_dict['integratedIntensity']
    num_dots = len(xs)
    gaussian_base = gaussian_gen(sigma_xy, sigma_z)
    for i in range(num_dots):
        x, y, z = xs[i], ys[i], zs[i]
        #  TODO
