"""Helper functions for masks."""

import os
import cv2
import numpy as np


# TODO: apply custom colormap
# https://stackoverflow.com/questions/52498777/apply-matplotlib-or-custom-colormap-to-opencv-image
def mask_overlay_image(img_path, mask, dst_path, overwrite=False):
    """Overlay mask on image and save to file.

    Args:
        img_path (str): path to image.
        mask (np.ndarray): mask array of shape (row, col), background value is 0.
        dst_path (str): destination path.
        overwrite (bool): whether to overwrite existing file.

    Returns:
        str: destination path.
    """
    if os.path.exists(dst_path) and not overwrite:
        print(f'File {dst_path} already exists.')
        return dst_path
    
    # mask image
    mask_img = cv2.cvtColor(mask.astype(np.uint8), cv2.COLOR_BGR2RGB)
    mask_img = cv2.applyColorMap(mask_img, cv2.COLORMAP_JET)
    
    # background image
    bkgd_img = cv2.imread(img_path)
    
    # place filter
    mask_place = (mask != 0).astype(np.uint8)  # mask locations
    bkgd_place = (mask == 0).astype(np.uint8)  # background locations

    # overlay
    heatmap = mask_place * np.transpose(mask_img, (2, 0, 1))  # mask
    heatmap = heatmap + bkgd_place * np.transpose(bkgd_img, (2, 0, 1))  # background
    heatmap = np.transpose(heatmap, (1, 2, 0))
    cv2.imwrite(dst_path, heatmap)
    return dst_path
