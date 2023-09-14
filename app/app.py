"""App."""

import toml
import json
import dash_leaflet as dl
import dash_uploader as du
from niceview.utils.dataset import ThorQuery
from dash import Dash, html, dcc

# config
config = toml.load('user/config.toml')
data_path = config['path']['data']
cache_path = config['path']['cache']
max_file_size = config['constant']['max_file_size']

# database information
with open('./db/db-info.json', 'r') as json_file:
    db_info = json.load(json_file)
data_extension = db_info['data_extension']
cache_extension = db_info['cache_extension']
cell_label_encoder = db_info['cell_label_encoder']
cell_label_cmap = db_info['cell_label_cmap']
primary_key_list = db_info['primary_key_list']

# arguments
with open('./user/args.json') as f:
    args = json.load(f)
sample_id = args['sampleId']
selected_cell_gene_name = args['selectedCellGeneName']
selected_spot_gene_name = args['selectedSpotGeneName']

# create instance and analyze data
thor = ThorQuery(
    data_path,
    cache_path,
    data_extension,
    cache_extension,
    cell_label_encoder,
    cell_label_cmap,
    primary_key_list,
)
thor.cell_gis(
    sample_id,
    selected_cell_gene_name,
    label_analysis=True,
)
thor.spot_gis(
    sample_id,
    selected_spot_gene_name,
)
thor.wsi_gis(sample_id)

# gis
cell_gene_client, cell_gene_layer = thor.gis_client_and_layer(sample_id, 'gis-blend-cell-gene-img')
cell_type_client, cell_type_layer = thor.gis_client_and_layer(sample_id, 'gis-blend-cell-type-img')
spot_gene_client, spot_gene_layer = thor.gis_client_and_layer(sample_id, 'gis-blend-spot-gene-img')
wsi_client, wsi_layer = thor.gis_client_and_layer(sample_id, 'gis-wsi-img')

# app
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
            id='user-input',
            children=[
                html.Div(
                    id='data-upload',
                    children=[
                        du.Upload(id='upload-image', max_file_size=max_file_size),
                    ],
                ),
                html.Div(
                    id='data-select',
                    children=[
                        dcc.Dropdown(placeholder='Select sample'),
                        dcc.Dropdown(placeholder='Select type of visualization'),
                        dcc.Dropdown(placeholder='Select gene'),
                    ],
                ),
            ],
        ),
        html.Div(
            id='trigger',
            children=[
                html.Button('Submit', id='submit-button', n_clicks=0),
            ],
        ),
        html.Div(
            id='app-container',
            children=[
                dl.Map(
                    id='cell-gene',
                    children=[
                        dl.TileLayer(
                            url=wsi_layer.url,
                            maxZoom=wsi_client.max_zoom,
                            minZoom=wsi_client.default_zoom,
                        ),
                        dl.LayersControl(
                            [
                                dl.Overlay(
                                    dl.TileLayer(
                                        url=cell_gene_layer.url,
                                        maxZoom=cell_gene_client.max_zoom,
                                        minZoom=cell_gene_client.default_zoom,
                                    ),
                                    name=selected_cell_gene_name,
                                ),
                            ],
                        ),
                        dl.FullScreenControl(),
                        dl.FeatureGroup(
                            [
                                dl.EditControl(),
                            ],
                        ),
                    ],
                    center=[cell_gene_client.center()[0], cell_gene_client.center()[1]],
                    zoom=cell_gene_client.default_zoom + 0.5,
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
