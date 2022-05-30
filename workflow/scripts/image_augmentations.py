"""
Code conceptualized from model notebook.

Functions for image augmentation to improve model performance.

March, 2022

Ye Yuan (yeyu@umich.edu)

ALL AUGMENTATION FUNCTIONS MUST TAKE IN XYCZ INPUT IMAGE AND OUTPUTS THE SAME SHAPE AND DTYPE IMAGE.
ALL AUGMENTATION FUNCTIONS MUST SUPPORT AUGMENTATION OF MASK (dtype=np.uint8) AND IMAGE

"""

import numpy as np
import cv2 as cv
from image_operators import imStackAutoAdjust


# Helper functions for augmentations
#       Note - at a minimum, our images needs to be resized such that a single CNN model suffices.


def augmentations_resize(array, target_pixels):
    """
    This function takes in shape = xycz images and scale its x&y both to target_pixels in pixel count.
    """
    zcount = array.shape[3]
    result = np.empty((target_pixels, target_pixels, array.shape[2], array.shape[3]), dtype=array.dtype)
    for z in range(zcount):
        currentSlice = array[:, :, :, z]
        currentSlice = cv.resize(currentSlice, dsize=(target_pixels, target_pixels), interpolation=cv.INTER_AREA)
        # Note - the slicing method above will lead to dimension drop if any of the dimension is 1.
        # Therefore, for multichannel shape = xyc; for singlechannel (mask) shape = xy
        # For our case, mask input will drop its channel dimension. We need to add it back.
        if len(currentSlice.shape) == 2:
            currentSlice = np.expand_dims(currentSlice, axis=-1)
        result[:, :, :, z] = currentSlice
    return result


def augmentations_enhanceThenResize(array, target_pixels, **kwargs):
    """
    This function enhances the image contrast and then resize the image
    Contrast enhancement is done via image_operators.imStackAutoAdjust
    """
    # First, check whether we are looking at mask images, which should be uint8
    if array.dtype == np.uint8:
        # For mask image, we don't need contrast enhancement
        return augmentations_resize(array, target_pixels)
    else:
        # For real data, we perform contrast enhancement
        ccount = array.shape[2]
        enhanced = np.empty(array.shape, dtype=array.dtype)
        for c in range(ccount):
            currentStack = array[:, :, c, :]
            enhanced[:, :, c, :] = imStackAutoAdjust(currentStack, **kwargs)
        return augmentations_resize(array, target_pixels)
