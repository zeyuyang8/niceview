import rasterio
import numpy as np

# get lon, lat from rasterio
georef_img_path = './plots/raster.tiff'
geo_ref_img = rasterio.open(georef_img_path)
lon, lat = geo_ref_img.xy(5000, 5000)
lon_min, lat_min, lon_max, lat_max = geo_ref_img.bounds
num_digits = len(str(lon)) - len(str(int(lon))) - 1
print(f'lon min: {lon_min}, lon max: {lon_max}')
print(f'lat min: {lat_min}, lat max: {lat_max}')
print(f'number of digits after decimal point: {num_digits}')

# make meshgrid
xs = np.linspace(lon_min, lon_max, geo_ref_img.width, dtype=np.float64)
ys = np.linspace(lat_min, lat_max, geo_ref_img.height, dtype=np.float64)
xx, yy = np.meshgrid(xs, ys)
y_max, x_max = xx.shape
print(f'y max: {y_max}, x max: {x_max}')