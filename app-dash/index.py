"""Index page for the app."""

import os
import dash
from dash import html
import dash_leaflet as dl
from localtileserver import TileClient, get_leaflet_tile_layer

# constants
DATA_PATH = '/home/tom/github/niceview/examples/data/'
PLOTS_PATH = '/home/tom/github/niceview/examples/plots/'

# args
sample_id = 'gt-iz-p9-rep2'
cells_selected_gene_name = 'ENSG00000065534'
spots_selected_gene_name = 'ENSG00000037280'

# cache
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

# client
cells_client = TileClient(
    os.path.join(PLOTS_PATH, cache['gis-blend-cells']),
    cors_all=True,
)
cells_tile_layer = get_leaflet_tile_layer(cells_client)

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
                dl.Map(
                    [
                        dl.TileLayer(url=cells_tile_layer.url),
                        dl.FullScreenControl(),
                        dl.ScaleControl(position='bottomleft'),
                    ],
                    center=[cells_client.center()[0], cells_client.center()[1]],
                    zoom=cells_client.default_zoom,
                    style={'width': '60vh', 'height': '60vh', 'margin': 'auto', 'display': 'block'},
                    attributionControl=False,
                ),
            ],
        ),
    ],
)

if __name__ == '__main__':
    app.run_server(debug=True)
