"""
    Performs image operations needed for segmentation model.
"""


import numpy as np
import cv2 as cv
import tifffile as tf
import math


# Morphological operations on image nparrays
def Z_unravel(img):
    """Unravels 3rd or 4th dimension of a given nparray and returns a list"""
    # Input - array, shape = xyz or xycz
    # Output - list of length z containing single Z stacks, each array shape = xy
    if len(img.shape) == 3:
        zcount = img.shape[2]
        return [img[:, :, i] for i in range(zcount)]
    if len(img.shape) == 4:
        zcount = img.shape[3]
        return [img[:, :, :, i] for i in range(zcount)]


def Z_MIP(img, z_position, z_spacing, method="onesided"):
    """Performs MIP operation on single-channel image shape=xyz at z=z_position and spacing=z_spacing"""
    # Maximum intensity projection along Z direction
    # Inputs
    #    img - single-channel array of shape xyz
    #    z_position - int specifying center z-stack, 0-based
    #    z_spacing - int specifying how many z-stack to consider for MIP
    #    method - specifies what to compute
    #      "onesided" - returns two single-channel array of shape xyz
    #            [0] will be MIP with z-slices [z_position - z_spacing, z_position] inclusive
    #            [1] will be MIP with z-slices [z_position, z_position + z_spacing] inclusive
    # Output - list of single-channel array of shape xy

    if method == "onesided":
        # determine z-stacks for MIP, considering boundaries
        zmax = img.shape[2] - 1
        zmins = [max(z_position - z_spacing, 0), z_position]
        zmaxs = [z_position, min(z_position + z_spacing, zmax)]
        projs = [img[:, :, zmins[0]:(zmaxs[0] + 1)],
                 img[:, :, zmins[1]:(zmaxs[1] + 1)]]
        projs = [np.max(proj, axis=2).astype(np.float32) for proj in projs]
        return projs

    raise ValueError("Z_MIP error: method not implemented")


def Z_MIP_Physical_Stack(img, z_spacing_physical, method="onesided", physicalSizeZ=None):
    """
        This function is based on Z_MIP. It performs slice-by-slice MIP with the given physical z-spacing.
        :param img: path to the image TIF file OR a np.array. Image must be XYZ monochromatic stack.
            The tif file MUST have PhysicalSize metadata saved by TiffFile.
            If np.array is provided, physicalSizeZ MUST be passed as parameter.
        :param z_spacing_physical: z-spacing in physical unit.
        :param method: method to run the MIP.
            onesided means that for each z location, two MIPs (z - z_spacing) & (z + z_spacing) will be created.
        :return: np.array of size XYCZ, where C = 3 containing MIP information depending on method.
    """
    if type(img) == str:
        img = tf.TiffFile(img)
        try:
            physicalSizeZ = float(img.shaped_metadata[0]["PhysicalSizeZ"])
        except:
            raise ValueError("Input file does not hold valid physical size information.")
        img = img.asarray()
    elif type(img) == np.ndarray:
        if physicalSizeZ is None:
            raise ValueError("Z_MIP must have valid physicalSizeZ value given array input.")
    else:
        raise ValueError("Z_MIP must input file path OR numpy array.")
    try:
        z_spacing = math.ceil(z_spacing_physical / physicalSizeZ)
    except:
        raise ValueError("physicalSizeZ cannot be zero.")
    dimx = img.shape[0]
    dimy = img.shape[1]
    dimz = img.shape[2]
    if method == "onesided":
        newstack = np.empty((dimx, dimy, 3, dimz), dtype=img.dtype)
        for idxZ in range(dimz):
            MIPs = Z_MIP(img, idxZ, z_spacing, "onesided")
            newstack[:, :, 1, idxZ] = img[:, :, idxZ]
            newstack[:, :, 0, idxZ] = MIPs[0]
            newstack[:, :, 2, idxZ] = MIPs[1]
        return newstack
    else:
        raise ValueError("Unsupported MIP method")


def imfill(img):
    """A simple hole-filling algorithm for uint8 single-slice mask image.
    Mask value is presumed to be 1. With apparent limitations."""
    # Input: a single-slice mask image, shape = xy, dtype = uint8
    # Output: same shape filled mask image
    # Note - this function assumes that object are marked by value = 1
    # Known limitation: if you have a nucleus part of topleft or botright corner,
    # the 0-valued interior part will not be filled by this function.

    # First, identify the boundary zero-value pixels
    zeroPixelIndices = np.nonzero(img == 0)
    zeroPixelCount = zeroPixelIndices[0].shape[0]

    # --------- NEED TO CHECK THIS !!!! ------------
    topLeft = (zeroPixelIndices[1][0], zeroPixelIndices[0][0])
    botRight = (zeroPixelIndices[1][zeroPixelCount - 1], zeroPixelIndices[0][zeroPixelCount - 1])

    # topLeft = (zeroPixelIndices[0][0], zeroPixelIndices[1][0])
    # botRight = (zeroPixelIndices[0][zeroPixelCount-1], zeroPixelIndices[1][zeroPixelCount-1])

    zeroPixelsForFill = (topLeft, botRight)

    # Next, perform cv.floodfill with the pixels
    filledMasks = np.empty((img.shape[0], img.shape[1], len(zeroPixelsForFill)), dtype=img.dtype)
    img_mask = np.zeros((img.shape[0] + 2, img.shape[1] + 2),
                        np.uint8)  # This mask means that no stopping boundaries for cv.floodFill
    for i in range(len(zeroPixelsForFill)):
        pixel = zeroPixelsForFill[i]
        currentSlice = img.copy()
        # Our fill will be 255
        # Original mask values are either 0 or 1
        cv.floodFill(image=currentSlice, mask=img_mask, seedPoint=pixel, newVal=255)
        filledMasks[:, :, i] = np.uint8(
            currentSlice != 255)  # With this, the filled background will be 0 while other places will be 1
    filledMaskFinal = np.ones(img.shape, dtype=img.dtype)
    for i in range(len(zeroPixelsForFill)):
        filledMaskFinal = np.bitwise_and(filledMaskFinal, filledMasks[:, :, i])
    return filledMaskFinal


def imAutoAdjust(img, method="simple", parameters=None):
    """Implements multiple deterministic contrast enhancement algorithms.
    Requires input as grayscale, single slice image, shape = xy
    Method = simple implements MATLAB's imadjust algorithm.
    Method = CLAHE is a thin wrapper of opencv's CLAHE algorithm.
    This function returns adjusted image and do NOT modify in-place.
    Outputs uint16 image."""
    # For method = simple, a single parameter saturation_percentage needs to be provided
    # For method = CLAHE, parameters dbl clipLimit and tuple_xy tileGridSize needs to be provided
    # Leave parameters None to adopt the default settings
    if parameters is None:
        parameters = dict()
    if method == "simple":
        # Determines parameter
        if len(parameters) == 0:
            saturation_percentage = 3  # ------ bottom 3% will be dark value while top 3% will be white value
        else:
            saturation_percentage = parameters['saturation_percentage']
        # Gets the pixel intensity positions for dark and white
        dark = np.percentile(img, q=saturation_percentage)
        white = np.percentile(img, q=100 - saturation_percentage)
        # Copy image, scale and clip it, and then generate output uint16 image
        img = img.copy()
        img = (img - dark) / (white - dark)  # Note - numpy will auto-convert dtype to float64
        img = np.clip(img, 0, 1)
        img = (img * np.iinfo(np.uint16).max).astype(np.uint16)
    elif method == "CLAHE":
        # Determines parameter
        if len(parameters) == 0:
            clipLimit = 20
            tileGridSize = (10, 10)
        else:
            clipLimit = parameters['clipLimit']
            tileGridSize = parameters['tileGridSize']
        clahe = cv.createCLAHE(clipLimit, tileGridSize)
        # Copy image and perform CLAHE
        img = img.copy()
        img = clahe.apply(img)
    else:
        raise ValueError("Method not implemented: imAutoAdjust")
    return img


def imStackAutoAdjust(img, method="simple", parameters=None):
    """A simple wrapper of imAutoAdjust.
    Performs imAutoAdjust for each Z-slice independently and returns the Z-stack.
    Requires input as grayscale, shape = xyz
    """
    ztotal = img.shape[2]
    result = np.empty(img.shape, dtype=np.uint16)
    for i in range(ztotal):
        result[:, :, i] = imAutoAdjust(img[:, :, i], method, parameters)
    return result


def imPadding(img):
    """
        accepts square/rectangular images. Do nothing if square, add padding if rectangular to make square image
    :param img: multichannel z-stack image of shape xycz
    :return: padding image of size xyc where x=y.
        If x<y, will add blank pixels in all channels [(x+1)...y ,:, :]
        If x>y, will add blank pixels in all channels [:, (y+1)...x, :]
    """
    (img_x, img_y) = img.shape[:2]
    if img_x > img_y:
        padding = ((0, 0), (0, img_x - img_y), (0, 0), (0, 0))
    else:
        padding = ((0, img_y - img_x), (0, 0), (0, 0), (0, 0))
    return np.pad(img, padding, mode='constant', constant_values=0)


def imPaddingRestore(img, original_x, original_y):
    """
        reverse operation of imPadding.
    :param img: multichannel z-stack image of shape x, x, c, z (i.e., square image)
    :param original_x: original x pixel count
    :param original_y: original y pixel count
    :return: image with padding removed, shape original_x, original_y, c
    """
    (img_x, img_y) = img.shape[:2]
    if img_x != img_y:
        raise ValueError("Paddng restore must have square input.")
    if img_x == original_x:
        #  x direction is not padded -> reverse y direction padding by slicing [:, :original_y, :]
        return img[:, :original_y, :, :]
    elif img_x == original_y:
        return img[:original_x, :, :, :]
    else:
        raise ValueError("Original must have the same ")