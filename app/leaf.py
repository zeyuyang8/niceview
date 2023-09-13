"""Index page for the app."""

import os
import numpy as np
import pandas as pd
from dash import Dash, html
import dash_leaflet as dl
from localtileserver import TileClient, get_leaflet_tile_layer
from scipy.sparse import load_npz
from niceview.utils.tools import normalize_array, select_col_from_name, txt_to_list
from niceview.utils.raster import geo_raster_to_meshgrid, index_to_meshgrid_coord

# constants
DATA_PATH = '/home/tom/github/niceview/examples/data/'
PLOTS_PATH = '/home/tom/github/niceview/examples/plots/'
CMIN = 0
CMAX = 255
SCALE_MAX_WIDTH = 200

# args.data
sample_id = 'gt-iz-p9-rep2'
cells_selected_gene_name = 'ENSG00000065534'
spots_selected_gene_name = 'ENSG00000037280'

# args.plot
window_size = (600, 600)

# file and cache
files = {
    'cells-gene-names': '-'.join([sample_id, 'cells-gene-names.txt']),
    'cells-gene': '-'.join([sample_id, 'cells-gene.npz']),
    'cells-info': '-'.join([sample_id, 'cells-info.csv']),
    'img': '-'.join([sample_id, 'img.tiff']),
    'mask-filtered-relabeled': '-'.join([sample_id, 'mask-filtered-relabeled.npz']),
    'mask': '-'.join([sample_id, 'mask.npz']),
    'spots-gene-names': '-'.join([sample_id, 'spots-gene-names.txt']),
    'spots-gene': '-'.join([sample_id, 'spots-gene.npz']),
    'spots-info': '-'.join([sample_id, 'spots-info.csv']),
}
cache = {
    'blend-cells-gene': '-'.join([sample_id, 'blend-cells-gene.png']),
    'blend-cells-type': '-'.join([sample_id, 'blend-cells-type.png']),
    'blend-spots-gene': '-'.join([sample_id, 'blend-spots-gene.png']),
    'mask-cells-type': '-'.join([sample_id, 'mask-cells-type.png']),
    'mask-cells-gene': '-'.join([sample_id, 'mask-cells-gene.png']),
    'spots-gene': '-'.join([sample_id, 'spots-gene.png']),
    'gis-blend-cells': '-'.join([sample_id, 'gis-blend-cells.tiff']),
    'gis-blend-spots': '-'.join([sample_id, 'gis-blend-spots.tiff']),
}

# cells client
cells_client = TileClient(
    os.path.join(PLOTS_PATH, cache['gis-blend-cells']),
    cors_all=True,
)
cells_tile_layer = get_leaflet_tile_layer(cells_client)

# # cells data processing
# cells_info = pd.read_csv(os.path.join(DATA_PATH, files['cells-info']))
# cells_centroids = np.round(cells_info[['y', 'x']].values).astype(int)
# cells_meshgrid = geo_raster_to_meshgrid(
#     os.path.join(os.path.join(PLOTS_PATH, cache['gis-blend-cells'])),
# )
# cells_lons, cells_lats = index_to_meshgrid_coord(cells_centroids, cells_meshgrid)
# cells_gene = load_npz(os.path.join(DATA_PATH, files['cells-gene']))  # scipy.sparse.csr.csr_matrix
# cells_gene_names = txt_to_list(os.path.join(DATA_PATH, files['cells-gene-names']))
# cells_selected_gene = select_col_from_name(cells_gene, cells_gene_names, cells_selected_gene_name)
# cells_selected_gene_normalized = normalize_array(cells_selected_gene, 1, CMAX)

# app starts here
app = Dash(
    __name__,
    meta_tags=[
        {
            'name': 'viewport',
            'content': 'width=device-width, initial-scale=1',
        },
    ],
)
app.title = 'Nice View'
app.layout = html.Div(
    id='root',
    children=[
        html.Div(
            id='header',
            children=[
                html.H4('Nice View'),
                html.P('A simple Dash app.'),
            ],
        ),
        html.Div(
            id='app-container',
            children=[
                dl.Map(
                    children=[
                        dl.TileLayer(
                            url=cells_tile_layer.url,
                            maxZoom=cells_client.max_zoom,
                            minZoom=cells_client.default_zoom,
                        ),
                        dl.FullScreenControl(),
                        dl.ScaleControl(imperial=False, metric=True, maxWidth=SCALE_MAX_WIDTH),
                        dl.FeatureGroup(
                            [
                                dl.EditControl(),
                            ],
                        ),
                    ],
                    center=[cells_client.center()[0], cells_client.center()[1]],
                    zoom=cells_client.default_zoom + 0.5,
                    style={'height': '70vh', 'margin': 'auto', 'display': 'block'},
                    attributionControl=False,
                ),
            ],
            style={'width': '70vh', 'height': '70vh'},
        ),
    ],
)


# run app
if __name__ == '__main__':
    app.run_server(debug=True)
