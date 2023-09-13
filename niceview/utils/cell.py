"""Cell."""

import numpy as np
import skimage.measure
from sklearn.neighbors import NearestNeighbors


def get_nuclei_pixels(cm, ad_cell_pos):
    """Get nuclei pixels.
    
    Args:
        cm (np.ndarray): Cell mask.
        ad_cell_pos (np.ndarray): Adherent cell positions.
    
    Returns:
        list: Nuclei region pixels.
    """
    tol = 1e-6
    regions = skimage.measure.regionprops(cm)
    xy = np.array([r.centroid[::-1] for r in regions])
    nbrs = NearestNeighbors(n_neighbors=1).fit(xy)
    distance, indices = nbrs.kneighbors(ad_cell_pos)
    c_index = indices[np.where(distance <= tol)]
    nuclei_region_pixels = [
        (regions[n].coords[:, 0], regions[n].coords[:, 1]) for n in c_index
    ]
    return np.array(nuclei_region_pixels, dtype=object)


def paint_regions(image_shape, matched_regions, cell_colors_list):
    """Paint regions.
    
    Args:
        image_shape (tuple): Image shape.
        matched_regions (list): Matched regions.
        cell_colors_list (list): Cell colors list.
    
    Returns:
        np.ndarray: Filled image.
    """
    filled = np.ma.masked_all(image_shape)
    for i, r in enumerate(matched_regions):
        cc, rr = r
        filled[cc, rr] = cell_colors_list[i]
    return filled
