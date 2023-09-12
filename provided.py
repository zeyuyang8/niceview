import os
import cv2
import numpy as np
import pandas as pd
from scipy.sparse import load_npz


# helper functions
def normalize_array(arr, new_min, new_max):
    min_val = np.min(arr)
    max_val = np.max(arr)
    normalized_arr = (arr - min_val) / (max_val - min_val) * (new_max - new_min) + new_min
    return normalized_arr

def mask_filter_relabel(mask_path, matched_regions, labels):
    mask = load_npz(mask_path)
    mask = mask.tocsr()[:, :].todense()
    # TODO: increase speed of `paint_regions`
    mask_filtered_relabeled = paint_regions(mask.shape, matched_regions, cell_colors_list=labels)
    return mask_filtered_relabeled.data

def hex_to_rgb(hex_color):
    # Remove the '#' symbol if it's present
    if hex_color.startswith("#"):
        hex_color = hex_color[1:]

    # Convert each pair from hexadecimal to decimal
    r = int(hex_color[0:2], 16)
    g = int(hex_color[2:4], 16)
    b = int(hex_color[4:6], 16)

    # Return the RGB values as a tuple
    return (r, g, b)

def discrete_cmap_from_hex(id_to_hex_dict):
    rgb_cmap = {k: hex_to_rgb(v) for k, v in id_to_hex_dict.items()}
    rgb_cmap = np.array([rgb_cmap[i] for i in range(1, len(rgb_cmap) + 1)])
    return rgb_cmap

def apply_custom_cmap(img_gray, cmap):
    lut = np.zeros((256, 1, 3), dtype=np.uint8)
    # rgb
    lut[1: len(cmap) + 1, 0, 0] = cmap[:, 0]
    lut[1: len(cmap) + 1, 0, 1] = cmap[:, 1]
    lut[1: len(cmap) + 1, 0, 2] = cmap[:, 2]
    # apply
    img_rgb = cv2.LUT(img_gray, lut)
    return img_rgb

def mask_to_image(mask, cmap):
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
