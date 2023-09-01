"""Index page for niceview."""

import os
import dash
from dash import html
import dash_core_components as dcc
from niceview.utils.mask import sparse_npz_to_array, mask_overlay_image, mask_to_bbox
from niceview.utils.raster import geo_ref_raster, geo_raster_to_meshgrid, index_to_meshgrid_coord
from niceview.pyplot.geoplot import scatter_overlay_georef_img

# configurations
DATA_PATH = '../examples/data/'

# TODO: automate this
mask_name = 'gt-iz-p9-rep2-mask.npz'  # variable
mask_overlay_name = 'gt-iz-p9-rep2-img-masked.png'  # variable
img_name = 'gt-iz-p9-rep2-img.tiff'  # variable
georef_masked_img_name = 'gt-iz-p9-rep2-img-masked.tiff'  # variable

# load data
mask = sparse_npz_to_array(os.path.join(DATA_PATH, mask_name))
mask_overlay_path = mask_overlay_image(
    os.path.join(DATA_PATH, img_name),
    mask,
    os.path.join(DATA_PATH, mask_overlay_name),
    overwrite=False,
)
georef_img_path = geo_ref_raster(
    mask_overlay_path,
    os.path.join(DATA_PATH, georef_masked_img_name),
    overwrite=False,
)
mask_info_dict = mask_to_bbox(mask)
centroids = mask_info_dict['centroids']
meshgrid = geo_raster_to_meshgrid(georef_img_path)
lons, lats = index_to_meshgrid_coord(centroids, meshgrid)

# plot image
fig = scatter_overlay_georef_img(
    georef_img_path,
    lons,
    lats,
    {
        'marker_size': 1,
        'window_width': 700,
        'window_height': 700,
    },
)

# initialize the app
app = dash.Dash(
    __name__,
    meta_tags=[
        {
            'name': 'viewport',
            'content': 'width=device-width, initial-scale=1',
        },
    ],
)
app.title = 'Nice View'

# layout
app.layout = html.Div(
    # root
    id='root',
    children=[
        # header
        html.Div(
            id='header',
            children=[
                html.H4('Nice View'),
                html.P('A simple Dash app.'),
            ],
        ),
        # app container
        html.Div(
            id='app-container',
            children=[
                # left column
                html.Div(
                    id='left-column',
                    children=[
                        html.P('Select view:'),
                        html.Div(
                            id='graph-container',
                            children=[
                                dcc.Graph(
                                    id='graph',
                                    figure=fig,
                                ),
                            ],
                        ),
                    ],
                ),
                # right column
                html.Div(
                    id='right-column',
                    children=[
                        html.P('View goes here.'),
                    ],
                ),
            ],
        ),
    ],
)


if __name__ == '__main__':
    app.run_server(debug=True)
