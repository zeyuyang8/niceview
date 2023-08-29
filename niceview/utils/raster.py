"""functions mainly based on rasterio."""

import os
import numpy as np
import rasterio
from rasterio.crs import CRS
from rasterio.warp import calculate_default_transform, reproject, Resampling
import warnings

MAX_PIXEL_VAL = 255


def rgba2rgb(rgba):
    """Convert RGBA image to RGB.

    Args:
        rgba (np.ndarray): RGBA image array.

    Returns:
        np.ndarray: RGB image array.
        
    Raises:
        ValueError: if input image does not have 4 channels.
    """
    ch, row, col = rgba.shape
    if ch != 4:
        raise ValueError('Input image must have 4 channels.')
    rgb = np.zeros((3, row, col), dtype='float32')
    r, g, b, a = (
        rgba[0, :, :], rgba[1, :, :], rgba[2, :, :], rgba[3, :, :],
    )

    a = np.asarray(a, dtype='float32') / MAX_PIXEL_VAL

    rgb[0, :, :] = r * a + (1.0 - a) * MAX_PIXEL_VAL
    rgb[1, :, :] = g * a + (1.0 - a) * MAX_PIXEL_VAL
    rgb[2, :, :] = b * a + (1.0 - a) * MAX_PIXEL_VAL

    return np.asarray(rgb, dtype='uint8')


def geo_ref_raster(
    img_path,
    dst_path,
    src_code=32632,
    dst_code=4326,
    affine_coefs=(0.1, 0.0, 0.0, 0.0, -0.1, 0.0),
    overwrite=True,
):
    """Georefence raster image.

    Args:
        img_path (str): path to image.
        dst_path (str): destination path.
        src_code (int): source EPSG code.
        dst_code (int): destination EPSG code.
        affine_coefs (tuple): affine transform coefficients.
        overwrite (bool): whether to overwrite existing file.
    
    Returns:
        str: path to georeferenced image.
    """
    # if already exists, return path
    if os.path.exists(dst_path) and not overwrite:
        print(f'File {dst_path} already exists.')
        return dst_path
    
    # read image
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        img = rasterio.open(img_path)
    
    # get image array, crs, and affine transform
    img_array = img.read()
    
    # convert RGBA to RGB
    if img_array.shape[0] == 4:
        img_array = rgba2rgb(img_array)

    # get source crs and affine transform
    crs = CRS.from_epsg(src_code)
    affine = rasterio.Affine(*affine_coefs)

    # georeference image and write temporary file
    temp_path = os.path.join(os.path.dirname(dst_path), 'temp.tiff')
    with rasterio.open(
        temp_path,
        'w',
        driver='GTiff',
        height=img_array.shape[1],
        width=img_array.shape[2],
        count=img_array.shape[0],
        dtype=img_array.dtype,
        crs=crs,
        transform=affine,
    ) as src:
        src.write(img_array)

    # reproject image to destination crs and write to file
    dst_crs = ':'.join(['EPSG', str(dst_code)])
    with rasterio.open(temp_path) as src:
        transform, width, height = calculate_default_transform(
            src.crs, dst_crs, src.width, src.height, *src.bounds,
            dst_height=src.height, dst_width=src.width,  # very important to keep same size
        )
        kwargs = src.meta.copy()
        kwargs.update(
            {
                'crs': dst_crs,
                'transform': transform,
                'width': width,
                'height': height,
            },
        )
        with rasterio.open(dst_path, 'w', **kwargs) as dst:
            for i in range(1, src.count + 1):
                reproject(
                    source=rasterio.band(src, i),
                    destination=rasterio.band(dst, i),
                    src_transform=src.transform,
                    src_crs=src.crs,
                    dst_transform=transform,
                    dst_crs=dst_crs,
                    dst_nodata=MAX_PIXEL_VAL,
                    resampling=Resampling.nearest,
                )

    # remove temporary file
    os.remove(temp_path)
    return dst_path
