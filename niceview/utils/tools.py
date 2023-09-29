"""Tools."""

from niceview.utils.cell import paint_regions
import numpy as np
from scipy.sparse import load_npz
import cv2
from shapely.geometry import Point, Polygon
import os
import rasterio
from pyproj import Geod
from scipy import ndimage
import scipy
import scanpy as sc


def txt_to_list(txt_file):
    """Read lines of a txt file to a list.

    Args:
        txt_file (str): txt file path

    Returns:
        lines (list of str): list of string of lines in the txt file
    """
    with open(txt_file, 'r') as txt:
        lines = txt.readlines()
        lines = [line.strip() for line in lines]
    return lines


def list_to_txt(lines, txt_file):
    """Write lines of a list to a txt file.

    Args:
        lines (list of str): list of string of lines to be written
        txt_file (str): txt file path
    
    Returns:
        txt_file (str): txt file path
    """
    with open(txt_file, 'w') as txt:
        for line in lines:
            txt.write(line)
            txt.write('\n')
    return txt_file


def select_col_from_name(matrix, name_list, name):
    """Select column from matrix by name.
    
    Args:
        matrix (np.ndarray): matrix of shape (row, col).
        name_list (list): list of names.
        name (str): name to select.
    
    Returns:
        np.ndarray: column of shape (row,).
    """
    idx = name_list.index(name)
    if isinstance(matrix, np.ndarray) and matrix.ndim == 2:
        return matrix[:, idx]
    return matrix.tocsr()[:, idx].todense()


def normalize_array(arr, new_min, new_max):
    """Normalize array to [new_min, new_max].
    
    Args:
        arr (np.ndarray): array to be normalized.
        new_min (float): new minimum value.
        new_max (float): new maximum value.
    
    Returns:
        np.ndarray: normalized array.
    """
    min_val = np.min(arr)
    max_val = np.max(arr)
    normalized_arr = (arr - min_val) / (max_val - min_val) * (new_max - new_min) + new_min
    return normalized_arr


def mask_filter_relabel(mask_path, matched_regions, labels):
    """Filter mask by matched regions and relabel the mask.
    
    Args:
        mask_path (str): path to the mask file.
        matched_regions (list of int): list of matched regions.
        labels (list of int): list of labels.
        
    Returns:
        np.ndarray: filtered and relabeled mask.
    """
    mask = load_npz(mask_path)
    mask = mask.tocsr()[:, :].todense()
    # TODO: increase speed of `paint_regions`
    mask_filtered_relabeled = paint_regions(mask.shape, matched_regions, cell_colors_list=labels)
    return mask_filtered_relabeled.data


def hex_to_rgb(hex_color):
    """Hexadecimal to RGB.
    
    Args:
        hex_color (str): hexadecimal color.
    
    Returns:
        tuple: RGB values.
    """
    # Remove the '#' symbol if it's present
    if hex_color.startswith('#'):
        hex_color = hex_color[1:]

    # Convert each pair from hexadecimal to decimal
    hex_max = 16
    r = int(hex_color[0:2], hex_max)
    g = int(hex_color[2:4], hex_max)
    b = int(hex_color[4:6], hex_max)

    # Return the RGB values as a tuple
    return (r, g, b)


def discrete_cmap_from_hex(id_to_hex_dict):
    """Discrete colormap from hex.
    
    Args:
        id_to_hex_dict (dict): dictionary of id to hex.
    
    Returns:
        np.ndarray: discrete colormap.
    """
    rgb_cmap = {int(k): hex_to_rgb(v) for k, v in id_to_hex_dict.items()}
    rgb_cmap = np.array([rgb_cmap[i] for i in range(1, len(rgb_cmap) + 1)])
    return rgb_cmap


def apply_custom_cmap(img_gray, cmap):
    """Apply custom colormap to gray image.
    
    Args:
        img_gray (np.ndarray): gray image.
        cmap (np.ndarray): custom colormap.
    
    Returns:
        np.ndarray: colored image.
    """
    lut = np.zeros((256, 1, 3), dtype=np.uint8)
    # rgb
    lut[1: len(cmap) + 1, 0, 0] = cmap[:, 0]
    lut[1: len(cmap) + 1, 0, 1] = cmap[:, 1]
    lut[1: len(cmap) + 1, 0, 2] = cmap[:, 2]
    # apply
    img_rgb = cv2.LUT(img_gray, lut)
    return img_rgb


def mask_to_image(mask, cmap):
    """Convert mask to image.
    
    Args:
        mask (np.ndarray): mask.
        cmap (np.ndarray): colormap.
    
    Returns:
        np.ndarray: image.
    """
    if isinstance(cmap, int):
        # TODO: increase speed of the following three lines, as they are "overlapping"
        img_rgb = cv2.cvtColor(mask.astype(np.uint8), cv2.COLOR_BGR2RGB)
        img_rgb = cv2.applyColorMap(img_rgb, cmap)
        img_rgb = cv2.bitwise_and(img_rgb, img_rgb, mask=mask.astype(np.uint8))
    else:
        img_gray = cv2.cvtColor(mask.astype(np.uint8), cv2.COLOR_GRAY2BGR)
        img_rgb = apply_custom_cmap(img_gray, cmap)
    return img_rgb


def draw_circles(img_shape, centers, diameter, colors, cmap=cv2.COLORMAP_JET, thickness=-1):
    """Draw circles on image.
    
    Args:
        img_shape (tuple): image shape.
        centers (list of tuple): list of centers.
        diameter (list of int): list of diameters.
        colors (np.ndarray): colors.
        cmap (int or np.ndarray): colormap.
        thickness (int): thickness of the circle.
    
    Returns:
        np.ndarray: image with circles.
    """
    # black background
    canvas = np.zeros((img_shape[0], img_shape[1], 3))
    
    # color
    if isinstance(cmap, int):
        colors = cv2.cvtColor(colors.astype(np.uint8), cv2.COLOR_BGR2RGB)
        colors = cv2.applyColorMap(colors, cv2.COLORMAP_JET)
        colors = np.reshape(colors, (-1, 3))
    else:
        colors = cv2.cvtColor(colors.astype(np.uint8), cv2.COLOR_BGR2RGB)
        colors = apply_custom_cmap(colors, cmap)
        colors = np.reshape(colors, (-1, 3))

    # set diameter
    if isinstance(diameter, int):
        diameter = [diameter] * len(centers)
    
    # draw circles
    for center, d, color in zip(centers, diameter, colors):
        color = tuple(map(int, color))  # convert elements to int
        center = np.round(center).astype('int')
        radius = np.round(d / 2).astype('int')
        cv2.circle(canvas, center, radius, color, thickness)
    return canvas


# TODO: speed up `blend`
def blend(img_path, mask_path, mask_opacity):
    """Blend mask and image.
    
    Args:
        img_path (str): path to the image.
        mask_path (str): path to the mask.
        mask_opacity (float): opacity of the mask.
    
    Returns:
        np.ndarray: blended image.
    """
    mask_img = cv2.imread(mask_path)
    bkgd_img = cv2.imread(img_path)
    
    # blend part of background
    mask = cv2.cvtColor(mask_img, cv2.COLOR_BGR2GRAY)
    bkgd_blend = cv2.bitwise_and(bkgd_img, bkgd_img, mask=mask)
    
    # non-blend part of background
    inv_mask = (mask == 0).astype(np.uint8)
    bkgd_non_blend = cv2.bitwise_and(bkgd_img, bkgd_img, mask=inv_mask)
    
    mask_ovelay = cv2.addWeighted(mask_img, mask_opacity, bkgd_blend, 1.0 - mask_opacity, 0)
    whole_img = cv2.addWeighted(mask_ovelay, 1.0, bkgd_non_blend, 1.0, 0)
    return whole_img


def get_bounding_box(coords):
    """Get bounding box of coordinates.
    
    Args:
        coords (list of tuple): list of coordinates.
    
    Returns:
        tuple: bounding box.
    """
    x_coords, y_coords = zip(*coords)
    x1, y1 = min(x_coords), min(y_coords)
    x2, y2 = max(x_coords), max(y_coords)
    return x1, y1, x2, y2


def save_roi_data_img(coords, adata, img, home_dir):
    """Get roi from coordinates.
    
    Args:
        coords (list of tuple): list of coordinates.
        adata (anndata.AnnData): anndata.
        img (np.ndarray): image.
        home_dir (str): home directory.
    """
    for idx, coord in enumerate(coords):
        # save adata
        roi = Polygon(coord)
        locs = list(map(lambda x: roi.contains(Point(x)), adata.obsm['spatial']))
        to_keep = adata[locs].copy()
        h5ad_path = os.path.join(home_dir, f'roi-{idx}.h5ad')
        to_keep.write_h5ad(h5ad_path)
        
        # save image
        x1, y1, x2, y2 = get_bounding_box(coord)
        pts = np.array(coord, np.int32).reshape((-1, 1, 2))
        cv2.polylines(img, [pts], isClosed=True, color=(0, 0, 255), thickness=1)
        cropped_region = img[y1:y2, x1:x2]
        cv2.imwrite(os.path.join(home_dir, f'roi-{idx}.tiff'), cropped_region)


def get_factor(gis_img_path, actual_distance=1e-6):
    """Get factor for converting pixel distance to actual distance.
    
    Args:
        gis_img_path (str): path to the GIS image.
        actual_distance (float): actual distance in micrometer.
        
    Returns:
        float: factor.
    """
    with rasterio.open(gis_img_path) as src:        
        lat1, lon1 = src.xy(0, 0)
        lat2, lon2 = src.xy(0, 1)
    g = Geod(ellps='clrk66') 
    _, _, dist = g.inv(lon1, lat1, lon2, lat2)
    factor = actual_distance / dist * 10 ** 6  # first convert to meter then convert to micrometer
    return factor


def process_data(img_path, mask_path, adata_path):
    """Process data for big image.
    
    Args:
        img_path (str): Path to image.
        mask_path (str): Path to mask.
        adata_path (str): Path to adata.
    
    Returns:
        None
    """
    img = cv2.imread(img_path)
    height, width, _ = img.shape
    max_dim = max(height, width)
    
    ref_val = 10000
    if max_dim < ref_val:
        return None
    
    # calculate resize factor
    resize_factor = ref_val / max_dim
    resized_img = cv2.resize(img, (int(width * resize_factor), int(height * resize_factor)))
    
    # read mask
    mask = np.array(load_npz(mask_path).todense())
    original_shape = mask.shape
    target_shape = (int(height * resize_factor), int(width * resize_factor))  # (height, width) in numpy format

    # resize mask array
    resized_mask = ndimage.zoom(
        mask, (target_shape[0] / original_shape[0], target_shape[1] / original_shape[1]), order=0,
    )
    
    # resize adata
    adata = sc.read_h5ad(adata_path)
    adata.obsm['spatial'] = adata.obsm['spatial'] * resize_factor
    
    # write data
    cv2.imwrite('temp.png', resized_img)
    scipy.sparse.save_npz('temp.npz', scipy.sparse.csr_matrix(resized_mask))
    adata.write('temp.h5ad')
    return None
