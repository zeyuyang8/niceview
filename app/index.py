"""Index page for the app."""

import os
import dash
from dash import html
from dash import dcc
from localtileserver import TileClient, get_leaflet_tile_layer
import numpy as np
import pandas as pd
from niceview.pyplot.geoplot import scatter_overlay_gis
from scipy.sparse import load_npz
from niceview.utils.tools import normalize_array, select_col_from_name, txt_to_list
from niceview.utils.raster import geo_raster_to_meshgrid, index_to_meshgrid_coord

# constants
DATA_PATH = '/home/tom/github/niceview/examples/data/'
PLOTS_PATH = '/home/tom/github/niceview/examples/plots/'
CMIN = 0
CMAX = 255

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

# cells data processing
cells_info = pd.read_csv(os.path.join(DATA_PATH, files['cells-info']))
cells_centroids = np.round(cells_info[['y', 'x']].values).astype(int)
cells_meshgrid = geo_raster_to_meshgrid(
    os.path.join(os.path.join(PLOTS_PATH, cache['gis-blend-cells'])),
)
cells_lons, cells_lats = index_to_meshgrid_coord(cells_centroids, cells_meshgrid)
cells_gene = load_npz(os.path.join(DATA_PATH, files['cells-gene']))  # scipy.sparse.csr.csr_matrix
cells_gene_names = txt_to_list(os.path.join(DATA_PATH, files['cells-gene-names']))
cells_selected_gene = select_col_from_name(cells_gene, cells_gene_names, cells_selected_gene_name)
cells_selected_gene_normalized = normalize_array(cells_selected_gene, 1, CMAX)

# cells figure
cells_fig = scatter_overlay_gis(
    cells_client, cells_tile_layer,
    cells_lons, cells_lats,
    marker_size=5,
    marker_color=np.squeeze(np.array(cells_selected_gene)),  # ravel col vector
    marker_color_scale='jet',
    marker_opacity=0.0,
    window_height=window_size[0],
    window_width=window_size[1],
    title='Cells',
)

# spots client
spots_client = TileClient(
    os.path.join(PLOTS_PATH, cache['gis-blend-spots']),
    cors_all=True,
)
spots_tile_layer = get_leaflet_tile_layer(spots_client)

# spots data processing
spots_info = pd.read_csv(os.path.join(DATA_PATH, files['spots-info']))
spots_pos = spots_info[['x', 'y']].values
spots_diameter = spots_info['diameter'].values
spots_gene = load_npz(os.path.join(DATA_PATH, files['spots-gene']))  # scipy.sparse.csr.csr_matrix
spots_gene_names = txt_to_list(os.path.join(DATA_PATH, files['spots-gene-names']))
spots_selected_gene = select_col_from_name(spots_gene, spots_gene_names, spots_selected_gene_name)
spots_selected_gene_normalized = normalize_array(spots_selected_gene, 1, CMAX)

# spots figure
spots_fig = scatter_overlay_gis(
    spots_client, spots_tile_layer,
    spots_pos[:, 0], spots_pos[:, 1],
    marker_size=spots_diameter,
    marker_color=np.squeeze(np.array(spots_selected_gene)),  # ravel col vector
    marker_color_scale='jet',
    marker_opacity=0.0,
    window_height=window_size[0],
    window_width=window_size[1],
    title='Spots',
)

# app starts here
app = dash.Dash(
    __name__,
    meta_tags=[
        {
            'name': 'viewport',
            'content': 'width=device-width, initial-scale=1',
        },
    ],
    external_stylesheets=[
        'https://unpkg.com/leaflet.markercluster@1.4.1/dist/MarkerCluster.Default.css',
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
                dcc.Graph(
                    id='cells-fig',
                    figure=cells_fig,
                ),
                dcc.Graph(
                    id='spots-fig',
                    figure=spots_fig,
                ),
            ],
        ),
    ],
)

if __name__ == '__main__':
    app.run_server(debug=True)
