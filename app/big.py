"""App."""

import os
import json
import shutil
import toml
import dash_uploader as du
from niceview.utils.dataset import ThorQuery
from niceview.pyplot.leaflet import create_leaflet_map
from niceview.utils.tools import save_roi_data_img
from dash import Dash, html, dcc, Input, Output, State

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
selected_pathway = args['selectedPathway']

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

# compute
thor.cell_gis(
    sample_id,
    selected_cell_gene_name,
    label_analysis=False,
    heatmap_analysis=False,
    selected_pathway=None,
)
# thor.spot_gis(
#     sample_id,
#     selected_spot_gene_name,
# )
thor.wsi_gis(sample_id)
cell_adata, wsi_img = thor.get_cell_adata_and_img(sample_id)

# gis
_, cell_gene_layer = thor.gis_client_and_layer(sample_id, 'gis-blend-cell-gene-img')
# _, cell_type_layer = thor.gis_client_and_layer(sample_id, 'gis-blend-cell-type-img')
# _, cell_gene_heatmap_layer = thor.gis_client_and_layer(sample_id, 'gis-blend-cell-gene-heatmap-img')
# _, cell_pathway_heatmap_layer = thor.gis_client_and_layer(sample_id, 'gis-blend-cell-pathway-heatmap-img')
# _, spot_gene_layer = thor.gis_client_and_layer(sample_id, 'gis-blend-spot-gene-img')
wsi_client, wsi_layer = thor.gis_client_and_layer(sample_id, 'gis-wsi-img')
mapper = thor.get_coord_mapping(sample_id)
height, width = thor.get_sample_img_shape(sample_id)

# crate map
fig = create_leaflet_map(
    'map',
    wsi_client,
    wsi_layer,
    [
        # (spot_gene_layer, 'spot gene'), 
        (cell_gene_layer, 'cell gene'), 
        # (cell_type_layer, 'cell type'),
        # (cell_gene_heatmap_layer, 'cell gene heatmap'),
        # (cell_pathway_heatmap_layer, 'cell pathway heatmap'),
    ],
    cmax=thor.get_gene_max(sample_id, selected_cell_gene_name),
)

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
            id='userInput',
            children=[
                html.Div(
                    id='dataUpload',
                    children=[
                        du.Upload(id='uploadImage', max_file_size=max_file_size),
                    ],
                ),
                html.Div(
                    id='dataSelect',
                    children=[
                        dcc.Dropdown(placeholder='Select sample'),
                        dcc.Dropdown(placeholder='Select type of visualization'),
                        dcc.Dropdown(placeholder='Select gene'),
                    ],
                ),
            ],
        ),
        html.Div(
            id='mainView',
            children=fig,
            style={'width': '70vh', 'height': '70vh'},
        ),
        html.Br(),
        html.Div(
            [
                dcc.Input(
                    id='newFilename',
                    type='text',
                    placeholder='Please enter a name',
                ),
                html.Button('Save', id='saveButton', n_clicks=0),
                html.Div(id='outputMessage'),
            ],
        ),
        html.Div(
            id='roi',
            children=[
                html.Div(
                    id='geoCoords',
                    style={'display': 'none'},
                ),
                html.Div(
                    id='coordinates',
                ),
            ],
        ),
        html.Div(id='currentInfo'),
    ],
)


@app.callback(
    Output('outputMessage', 'children'),
    Input('saveButton', 'n_clicks'),
    State('newFilename', 'value'),
    prevent_initial_call=True,
)
def copy_and_rename_file(n_clicks, new_name):
    """Copy and rename file.
    
    Args:
        n_clicks (int): Number of clicks.
        new_name (str): New name.
    
    Returns:
        str: Output message.
    """
    dir_path = f'./user/{new_name}/'
    os.mkdir(dir_path)
    coord_path = os.path.join(dir_path, 'coords.json')
    shutil.copyfile('./user/coords.json', coord_path)
    with open(coord_path, 'r') as f:
        coords = json.load(f)
    
    save_roi_data_img(coords, cell_adata, wsi_img, dir_path)
    return f'data of {new_name} is saved.'


@app.callback(
    Output('geoCoords', 'children'),
    Output('coordinates', 'children'),
    Input('editControl', 'geojson'),
    prevent_initial_call=True,
)
def save_roi(drawn_geojson):
    """Save coordinates.
    
    Args:
        drawn_geojson (dict): GeoJSON.

    Returns:
        json file.
    """
    coords = []
    for region in drawn_geojson['features']:
        temp = region['geometry']['coordinates'][0]
        temp = [mapper(*point) for point in temp]
        temp = [[point[1], point[0]] for point in temp]  # y, x -> x, y
        coords.append(temp)
    
    geojson_name = './user/roi.json'
    with open(geojson_name, 'w') as f:
        json.dump(drawn_geojson, f)
    
    coordjson_name = './user/coords.json'
    with open(coordjson_name, 'w') as f:
        json.dump(coords, f)
    
    return f'{drawn_geojson}', f'{coords}'


@app.callback(
    Output('currentInfo', 'children'),
    Input('map', 'bounds'),
)
def display_current_info(viewport):
    """Display current zoom level and center.
    
    Args:
        viewport (dict): Viewport.
    
    Returns:
        str: Current zoom level and center.
    """
    viewport = [[point[1], point[0]] for point in viewport]
    viewport = [mapper(*point) for point in viewport]
    return f'{viewport}'


# run app
if __name__ == '__main__':
    app.run_server(debug=True)
