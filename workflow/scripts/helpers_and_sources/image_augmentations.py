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
# TODO: use package relative import for more pythonic code.
# noinspection PyUnresolvedReferences
# Search path includes scripts folder only.
from helpers_and_sources.image_operators import imStackAutoAdjust, imPadding, imPaddingRestore


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


def augmentations_enhanceThenResize(array, target_pixels, padding=False, **kwargs):
    """
    This function enhances the image contrast and then resize the image
    Contrast enhancement is done via image_operators.imStackAutoAdjust
    :param array: Input image, either mask or stack with shape xycz
    :param target_pixels: Target size of the image after resizing. x=y=target_pixels
    :param padding: bool, whether to perform image padding before resizing. Otherwise, rectangular images can distort
    :return: augmented image w/ contrast enhancement, optional padding, and resizing. shape xycz is maintained
        Note - if padding is specified, the image will be padded before resizing
    """
    # First, check whether we are looking at mask images, which should be uint8
    if array.dtype == np.uint8:
        # For mask image, we don't need contrast enhancement
        if padding:  # do padding if specified
            return augmentations_resize(imPadding(array), target_pixels)
        else:
            return augmentations_resize(array, target_pixels)
    else:
        # For real data, we perform contrast enhancement
        ccount = array.shape[2]  # num of channels
        zcount = array.shape[3]  # num of z-slices
        enhanced = np.empty(array.shape, dtype=array.dtype)
        if padding:  # If padding is used, we need to store into square shape after padding.
            #  However, the shape is not final target_pixels, but rather max of x and y dimensions
            padded_xy = max(array.shape[:2])
            padded = np.empty([padded_xy, padded_xy, ccount, zcount], dtype=array.dtype)
        for c in range(ccount):
            currentStack = array[:, :, c, :]
            enhanced[:, :, c, :] = imStackAutoAdjust(currentStack, **kwargs)
            if padding:  # perform padding before resize after enhancement for image
                padded[:, :, c, :] = imPadding(enhanced[:, :, c, :])
        if padding:
            return augmentations_resize(padded, target_pixels)
        else:  # If padding is not used, we resize the enhanced image.
            return augmentations_resize(enhanced, target_pixels)


def revaug_resize(array, original_pixel_x, original_pixel_y):
    """
        Reverse augmentation "augmentations_resize". Useful for scaling segmentation model's output masks.
    :param array: Input single-channel image, shape xyz
    :param original_pixel_x: original pixel count in x direction
    :param original_pixel_y: original pixel count in y direction
    """
    array_original = np.empty((original_pixel_x, original_pixel_y, array.shape[2]), dtype=array.dtype)
    for i in range(array.shape[2]):
        # ---------- NEED TO VERIFY THIS & GENERALIZE THIS INTO A FUNCTION ------------
        array_original[:, :, i] = cv.resize(array[:, :, i],
                                            dsize=(original_pixel_y, original_pixel_x), interpolation=cv.INTER_NEAREST)
    return array_original


def revaug_resizeThenUnpad(array, original_pixel_x, original_pixel_y):
    """
        Partially reverse augmentation "augmentations_enhanceThenResize(..., padding=True)".
        Useful for rescaling segmentation model output to originally-sized segmentation mask.
        ONLY supports single-channel array input.
    :param array: Input single-channel image, shape xyz
    :param original_pixel_x: original pixel count in x direction
    :param original_pixel_y: original pixel count in y direction
    """
    resize_xy = max(original_pixel_x, original_pixel_y)  # resize target pixel size
    result = revaug_resize(array, resize_xy, resize_xy)
    return imPaddingRestore(result, original_pixel_x, original_pixel_y)
