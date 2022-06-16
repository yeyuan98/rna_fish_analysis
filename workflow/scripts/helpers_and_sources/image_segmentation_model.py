"""
    Pytorch image segmentation models
"""

import numpy as np
import torch
from torch import nn
import segmentation_models_pytorch as smp
from segmentation_models_pytorch.losses import DiceLoss
import cv2 as cv


# SMP pretrained segmentation model, flexible architecture
#   Constraints: 3 channels, binary classification
class SegmentationModel(nn.Module):
    """
    This is the segmentation model class, based on the smp package pre-trained U-Net.
    """
    def __init__(self, encoder, weights):
        super(SegmentationModel, self).__init__()  # This is needed for any model inherited from nn.Module

        # Architecture of the model
        self.arc = smp.Unet(
            encoder_name=encoder,
            encoder_weights=weights,
            in_channels=3,  # How many channels do we have, RGB image, C = 3
            classes=1,  # Binary classification, only one class except background
            activation=None
            # We do not use any sigmoid or other activation functions for final layer of generating output
        )

    def forward(self, images, masks=None):
        logits = self.arc(images)

        # If masks are provided (i.e., annotated images) calculate loss in forward run
        if masks is not None:
            loss1 = DiceLoss(mode='binary')(logits, masks)
            loss2 = nn.BCEWithLogitsLoss()(logits, masks)  # BCE = binary cross entropy
            return logits, loss1 + loss2  # Combining these two losses by adding together to get the final loss function

        # If masks are not provided (i.e., making inferences)
        # only calculate output (without activation, that would be logits in final layer)
        return logits


def getMask(model, image, device, original_xy_pixels):
    """
        Performs segmentation prediction task for a Z-stack
        :param model: pytorch nn.Module model for a single Z-slice shape=xyc
        :param image: z-stack augmented image shape=xycz where x=y=target_pixels
        :param device: what device to perform prediction
        :param original_xy_pixels: dict {x: <int>, y: <int>} for original image size.
            original sizes are needed to perform resizing
        :return: final predicted mask, shape=xyz, dtype=uint8
    """
    shapexy = image.shape[0]
    shapez = image.shape[3]
    pred_masks = np.empty((shapexy, shapexy, shapez), dtype=np.uint8)
    for i in range(shapez):
        img = image[:, :, :, i]  # shape = xyc (hwc)
        # Pytorch requires cxy
        img = np.transpose(img, (2, 0, 1))
        # Convert ndarray to tensor
        img = torch.Tensor(img.astype(np.int32))
        logits_mask = model(img.to(device).unsqueeze(0))  # this adds the batch dimension, shape = (1, c, h, w)
        pred_mask = torch.sigmoid(logits_mask)
        pred_mask = (pred_mask > 0.5) * 1.0  # Larger than 0.5 probability will give a mask value of 1.0
        pred_masks[:, :, i] = pred_mask.detach().cpu().squeeze(0).numpy().astype(np.uint8)
    original_x = original_xy_pixels["x"]
    original_y = original_xy_pixels["y"]
    pred_masks_original = np.empty((original_x, original_y, image.shape[3]), dtype=np.uint8)
    for i in range(pred_masks.shape[2]):
        # ---------- NEED TO VERIFY THIS & GENERALIZE THIS INTO A FUNCTION ------------
        pred_masks_original[:, :, i] = cv.resize(pred_masks[:, :, i],
                                                 dsize=(original_y, original_x), interpolation=cv.INTER_NEAREST)
    return pred_masks_original
