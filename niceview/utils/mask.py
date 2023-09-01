"""Helper functions for masks."""

import os
import cv2
import numpy as np
import skimage.measure
from scipy.sparse import load_npz


# TODO: apply custom colormap
# https://stackoverflow.com/questions/52498777/apply-matplotlib-or-custom-colormap-to-opencv-image
def mask_overlay_image(img_path, mask, dst_path, mask_opacity=0.3, colormap='jet', overwrite=False):
    """Overlay mask on image and save to file.

    Args:
        img_path (str): path to image.
        mask (np.ndarray): mask array of shape (row, col), background value is 0.
        dst_path (str): destination path.
        mask_opacity (float): opacity of mask.
        colormap (str): colormap of mask.
        overwrite (bool): whether to overwrite existing file.

    Returns:
        str: destination path.
    """
    if os.path.exists(dst_path) and not overwrite:
        print(f'File {dst_path} already exists.')
        return dst_path
    if colormap == 'jet':
        cmap = cv2.COLORMAP_JET
    
    # mask image
    mask_img = cv2.cvtColor(mask.astype(np.uint8), cv2.COLOR_BGR2RGB)
    mask_img = cv2.applyColorMap(mask_img, cmap)
    mask_blend = cv2.bitwise_and(mask_img, mask_img, mask=mask.astype(np.uint8))
    
    # background image
    bkgd_img = cv2.imread(img_path)
    bkgd_blend = cv2.bitwise_and(bkgd_img, bkgd_img, mask=mask.astype(np.uint8))
    inv_mask = (mask == 0).astype(np.uint8)
    bkgd = cv2.bitwise_and(bkgd_img, bkgd_img, mask=inv_mask)
    
    # overlay
    mask_ovelay = cv2.addWeighted(mask_blend, mask_opacity, bkgd_blend, 1.0 - mask_opacity, 0)
    whole_img = cv2.addWeighted(mask_ovelay, 1.0, bkgd, 1.0, 0)
    
    # save
    cv2.imwrite(dst_path, whole_img)
    return dst_path


def mask_to_bbox(mask):
    """Compute bounding boxes, centroid, and labels from mask.

    Args:
        mask (np.ndarray): mask array of shape (row, col), background value is 0.

    Returns:
        dict: dictionary of bounding boxes, centroids, and labels.
    """
    props = skimage.measure.regionprops(mask)
    bboxes = np.array([prop.bbox for prop in props])
    centroids = np.array([prop.centroid for prop in props])
    centroids = np.round(centroids).astype(int)
    labels = np.array([prop.label for prop in props])
    return {
        'bboxes': bboxes,
        'centroids': centroids,
        'labels': labels,
    }


def sparse_npz_to_array(npz_path):
    """Convert sparse matrix in npz format to array.
    
    Args:
        npz_path (str): path to npz file.

    Returns:
        np.ndarray: array of shape (row, col).
    """
    mask = load_npz(npz_path)
    mask = mask.toarray()
    return mask


def relabel_mask_bbox(mask, output_label):
    """Relabel mask to output_label (apply to bounding boxes).
    
    Args:
        mask (np.ndarray): mask array of shape (row, col), background value is 0.
        output_label (np.ndarray): output labels.
    
    Returns:
        np.ndarray: relabeled mask.
    """
    table = skimage.measure.regionprops_table(mask)
    relabeled_mask = skimage.util.map_array(mask, table['label'], output_label)
    return relabeled_mask
