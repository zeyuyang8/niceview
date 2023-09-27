"""Functions for calculating distance between two points."""

import rasterio
from pyproj import Geod


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
