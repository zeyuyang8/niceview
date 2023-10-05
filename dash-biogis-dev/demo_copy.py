import sys
sys.path.append("../")
import toml
import json
import dash
from dash import html
from dash import dcc,ctx
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate
import pandas as pd
import os
import io
import dash_uploader as du
import shutil
import dash_leaflet as dl
import uuid
from niceview.utils.dataset import ThorQuery
from niceview.pyplot.leaflet import create_leaflet_map,create_leaflet_map_output
from niceview.utils.convert import h5ad_converter
import dash_daq as daq



# session_id=uuid.uuid1()
# session_id=str(session_id)
session_id="hello-rep"
config = toml.load('../user/config.toml')
data_path = config['path']['data']
cache_path = config['path']['cache']
max_file_size = config['constant']['max_file_size']
# args.plot

with open('../user/args-default.json') as default_f:
        args_default = json.load(default_f)
with open('../user/args.json', 'w') as f:
        json.dump(args_default, f)


with open('../db/db-info-default.json', 'r') as json_file_d:
        db_info_default = json.load(json_file_d)

with open('../db/db-info.json', 'w') as json_file:
        json.dump(db_info_default, json_file)





# Initialize app
app = dash.Dash(
    __name__,
    meta_tags=[
        {"name": "viewport"}
    ],
    
)
app.title = "BioGIS"
server = app.server

du.configure_upload(app, r"../data_input_temp/tmp/",use_upload_id=True)

# Map tile client#
# Create TileClients and tile layers for both maps
def files_generate(sample_id):
    files = {
        'cells-gene-names': '-'.join([sample_id, 'cell-gene-name.txt']),
        'img': '-'.join([sample_id, 'wsi-img.tiff']),
        'spots-gene-names': '-'.join([sample_id, 'spot-gene-name.txt'])
    }
    return files

def cache_generate(sample_id,sample_id_file='',sample_id_gene_cell='',sample_id_gene_spot=''):
    cache = {
        'gis-img': '-'.join([sample_id, 'gis-wsi-img.tiff']),
        'gis-blend-cells': '-'.join([sample_id, 'gis-blend-cell-gene-img.tiff']),
        'gis-blend-spots': '-'.join([sample_id, 'gis-blend-spot-gene-img.tiff']),
        'gis-blend-cell-type': '-'.join([sample_id, 'gis-blend-cell-type-img.tiff']),

        'gis-img-file': '-'.join([sample_id_file, 'gis-wsi-img.tiff']),
        'gis-blend-cells-gene': '-'.join([sample_id_gene_cell, 'gis-blend-cell-gene-img.tiff']),
        'gis-blend-spots-gene': '-'.join([sample_id_gene_spot, 'gis-blend-spot-gene-img.tiff']),
        'gis-blend-cell-type-file': '-'.join([sample_id_file, 'gis-blend-cell-type-img.tiff'])
    }
    return cache

def get_parameter():
    id_list=[]
    with open('../db/db-info.json', 'r') as json_file:
        db_info = json.load(json_file)
    data_extension = db_info['data_extension']
    cache_extension = db_info['cache_extension']
    cell_label_encoder = db_info['cell_label_encoder']
    cell_label_cmap = db_info['cell_label_cmap']
    primary_key_list = db_info['primary_key_list']
    with open('../user/args.json') as f:
        args = json.load(f)
    sample_id = args['sampleId']
    sample_id_file=args['sampleIdFile']
    sample_id_gene_spot=args['sampleIdSpotGene']
    sample_id_gene_cell=args['sampleIdSpotGene']

    # arguments

    thor = ThorQuery(
        data_path,
        cache_path,
        data_extension,
        cache_extension,
        cell_label_encoder,
        cell_label_cmap,
        primary_key_list,
    )
    
    return thor,args

def get_wsi():
    thor,args=get_parameter()
    sample_id = args['sampleId']
    sample_id_file=args['sampleId']+'-'+args['fileName']
    sample_id_gene_spot=args['sampleIdSpotGene']
    sample_id_gene_cell=args['sampleIdSpotGene']

    args["sampleIdFile"]=sample_id_file

    with open('../user/args.json', 'w') as f:
        json.dump(args, f)

    thor.wsi_gis(sample_id)
    cache=cache_generate(sample_id,sample_id_file=sample_id_file)
    shutil.copy(os.path.join(cache_path,cache["gis-img"]),os.path.join(cache_path,cache["gis-img-file"]))

def calculation_cell():
    thor,args=get_parameter()
    

    sample_id = args['sampleId']
    sample_id_file=args['sampleIdFile']
    selected_cell_gene_name = args['selectedCellGeneName']
    
    sample_id_gene_cell=sample_id+"-"+selected_cell_gene_name
    
    args['sampleIdCellGene']=sample_id_gene_cell

    with open('../user/args.json', 'w') as f:
        json.dump(args, f)

    thor.cell_gis(
        sample_id,
        selected_cell_gene_name,
        label_analysis=True,
    )
    
    cache=cache_generate(sample_id,sample_id_gene_cell=sample_id_gene_cell,sample_id_file=sample_id_file)
    shutil.copy(os.path.join(cache_path,cache["gis-blend-cells"]),os.path.join(cache_path,cache["gis-blend-cells-gene"]))
    shutil.copy(os.path.join(cache_path,cache["gis-blend-cell-type"]),os.path.join(cache_path,cache["gis-blend-cell-type-file"]))

def calculation_spot():
    thor,args=get_parameter()

    sample_id = args['sampleId']
    selected_spot_gene_name = args['selectedSpotGeneName']
    sample_id_gene_spot=sample_id+"-"+selected_spot_gene_name
    print(sample_id_gene_spot,"calculation spot")

    args['sampleIdSpotGene']=sample_id_gene_spot

    with open('../user/args.json', 'w') as f:
        json.dump(args, f)

    thor.spot_gis(
        sample_id,
        selected_spot_gene_name,
    )
    
    cache=cache_generate(sample_id,sample_id_gene_spot=sample_id_gene_spot)
    shutil.copy(os.path.join(cache_path,cache["gis-blend-spots"]),os.path.join(cache_path,cache["gis-blend-spots-gene"]))


def visualization_img_input(data_path, cache_path):
    thor, args = get_parameter()
    sample_id_file = args['sampleIdFile']
    wsi_client, wsi_layer = thor.gis_client_and_layer(sample_id_file, 'gis-wsi-img')
    input_map = create_leaflet_map(
        'map-input',
        wsi_client,
        wsi_layer,
        [],
        'editControl'
    )
    return input_map

def visualization_img_all(data_path,cache_path):
    thor,args=get_parameter()
    sample_id_file=args['sampleIdFile']
    sample_id_gene_spot=args['sampleIdSpotGene']
    sample_id_gene_cell=args['sampleIdSpotGene']

    cell_gene_client, cell_gene_layer = thor.gis_client_and_layer(sample_id_gene_cell, 'gis-blend-cell-gene-img')
    wsi_client, wsi_layer = thor.gis_client_and_layer(sample_id_file, 'gis-wsi-img')
    cell_type_client, cell_type_layer = thor.gis_client_and_layer(sample_id_file, 'gis-blend-cell-type-img')
    spot_gene_client, spot_gene_layer = thor.gis_client_and_layer(sample_id_gene_spot, 'gis-blend-spot-gene-img')
    output_map = create_leaflet_map_output(
        'map-output',
        wsi_client,
        wsi_layer,
        [(spot_gene_layer, 'spot gene'), (cell_gene_layer, 'cell gene'), (cell_type_layer, 'cell type')],
        'geojson'
    )
    return output_map

def visualization_img_cell(data_path,cache_path,cell_type=True):
    thor,args=get_parameter()
    sample_id_file=args['sampleIdFile']
    sample_id_gene_cell=args['sampleIdSpotGene']

    cell_gene_client, cell_gene_layer = thor.gis_client_and_layer(sample_id_gene_cell, 'gis-blend-cell-gene-img')
    wsi_client, wsi_layer = thor.gis_client_and_layer(sample_id_file, 'gis-wsi-img')
    if cell_type == True:
        cell_type_client, cell_type_layer = thor.gis_client_and_layer(sample_id_file, 'gis-blend-cell-type-img')
        output_map = create_leaflet_map_output(
            'map-output',
            wsi_client,
            wsi_layer,
            [(cell_gene_layer, 'cell gene'), (cell_type_layer, 'cell type')],
            'geojson'
        )
    elif cell_type == False:
        cell_type_client, cell_type_layer = thor.gis_client_and_layer(sample_id_file, 'gis-blend-cell-type-img')
        output_map = create_leaflet_map_output(
            'map-output',
            wsi_client,
            wsi_layer,
            [(cell_gene_layer, 'cell gene')],
            'geojson'
        )
    return output_map


def visualization_img_spot(data_path,cache_path):
    thor,args=get_parameter()
    sample_id_file=args['sampleIdFile']
    sample_id_gene_spot=args['sampleIdSpotGene']

    # print(sample_id_gene_spot)

    wsi_client, wsi_layer = thor.gis_client_and_layer(sample_id_file, 'gis-wsi-img')
    spot_gene_client, spot_gene_layer = thor.gis_client_and_layer(sample_id_gene_spot, 'gis-blend-spot-gene-img')
    output_map = create_leaflet_map_output(
        'map-output',
        wsi_client,
        wsi_layer,
        [(spot_gene_layer, 'spot gene')],
        "geojson"
    )
    return output_map
    

# Create map components
get_wsi()
map_input = visualization_img_input(data_path,cache_path)

calculation_cell()
calculation_spot()
map_output = visualization_img_all(data_path,cache_path)

# App layout
app.layout = html.Div(
    
    children=[
        html.Link(rel='stylesheet', href="https://use.fontawesome.com/releases/v5.5.0/css/all.css", integrity="sha384-B4dIYHKNBt8Bc12p+WXckhzcICo0wtJAoU8YZTY5qE0Id1GSseTk6S+L3BlXeVIU", crossOrigin="anonymous" ),
        html.Link(rel='stylesheet', href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/4.7.0/css/font-awesome.min.css"),
        html.Div(id='header',
        children=[
        
        html.H4(id="logo",children='BioGIS'),
        html.Br(),
            html.A('A web application for visualizing and analyzing biological data',
                    id='description'
                ),

    ]),
        html.Div(
            className='container clearfix',id="graph",
            children=[
                html.Div(id='submit-container',children=[  
                        html.H5("Upload H&E image:"),
                        html.Br(),
                        html.Div(className="upload-data",children=[
                            du.Upload(id='upload-data-image',max_file_size=5000),
                        ]),
                        html.Div(id="status1"),
                        html.Br(),
                        html.Div(children=[
                            dcc.Dropdown(['Spot data','Cell data'],id="spot-cell-option", className='dropdown-input',placeholder="Choose Cell or Spot data type")
                            ]),
                        html.Br(),
                        html.Div(id="additional-data-box"),
                        html.Div(id="status2"),
                        html.Div(id="status3"),
                        html.Br(),
                        html.Div(id='output-data-upload'),
                        html.Br(),
                        html.H5("Choose type of visualization"),
                        html.Div(children=[dcc.Dropdown(['GO pathway','Gene Expression'],id="visual-type-container", className='dropdown-input',placeholder="Select Type")]),
                        html.Br(),html.Br(),html.Br(),
                        html.Div(id="gene-dropdown"),
                        html.Div(id='gene-chosen-output'),
                        html.Button('Start Visualize', className="button", id="analyze-button",n_clicks=0),
                        html.Div(id="hello")
                        ]),
                html.Div(
                    id='left-column',
                    children=[
                        html.Div(
                            id="input-image",
                            children=[
                                # TODO: Add something here
                                map_input,
                            ],
                        ),
                    ],
                ),
                html.Div(
                    id='right-column',
                    children=[
                        html.Div(
                            id="output-image",
                            children=[
                                # TODO: Add something here
                                map_output,
                            ],
                        ),
                    ],
                ),
            ],
        ),
    html.Div(
        id="footer",children=[
            html.A(html.I(className="fa fa-github", style={"font-size": "24px"}),
           href="https://github.com/GuangyuWangLab2021", target="_blank"),
    
            html.A(html.I(className="fa fa-linkedin-square", style={"font-size": "24px"}),
                href="https://www.linkedin.com/in/guangyu-wang-27696819b/", target="_blank"),

            html.A(html.I(className="fa fa-twitter", style={"font-size": "24px"}),
                href="https://twitter.com/Guangyu_Wang01", target="_blank"),
            html.P("©2023 by Wang lab.", className="text")
        ]
    )
    
    
    ],
)

            

#upload HE image
@du.callback(
    output=Output('input-image', 'children'),
    id='upload-data-image',
)
def upload_image(filenames):
    sample_id=session_id
    basename=os.path.splitext(os.path.basename(filenames[0]))[0]
    thor,args=get_parameter()

    # # Update args and write back to args.json
    args['sampleId']= sample_id
    args['fileName']= basename

    with open('../user/args.json', 'w') as f:
        json.dump(args, f)

    # with open('../db/db-info.json', 'r') as json_file:
    #     db_info = json.load(json_file)

    # # Update db_info and write back to db-info.json
    # db_info['primary_key_list'] = [sample_id]
    # with open('../db/db-info.json', 'w') as json_file:
    #     json.dump(db_info, json_file)
    
    files=files_generate(sample_id)
    shutil.copy(filenames[0],os.path.join(data_path,files["img"]))
    shutil.rmtree("../data_input_temp/tmp/")
    get_wsi()
    map_input=visualization_img_input(data_path,cache_path)
    return html.Div(
                    id="input-image",
                    children=[
                        # TODO: Add something here
                        map_input,
                    ]
                )

#choose cell or spot data
@app.callback(
    Output('additional-data-box', 'children'),
    Input('spot-cell-option', 'value')
)
def show_cell_spot_upload(value):
    if value == "Spot data":
        return html.Div(className="upload-data",children=[
            html.Br(),html.Br(),
            html.H5("Upload addition data for spot:"),
            html.Br(),
            du.Upload(id='upload-data-addition-spot',max_file_size=5000)
        ])
    elif value == "Cell data" :
        return html.Div(className="upload-data",children=[
            html.Br(),html.Br(),
            # daq.ToggleSwitch(color="lightgreen",id='cell-type-switch',label='Showing cell type?',labelPosition='top'),
            # html.Div(id='my-toggle-switch-output'),
            # html.Br(),
            html.H5("Upload addition data for cell:"),
            html.Br(),
            du.Upload(id='upload-data-addition-cell',max_file_size=5000)
        ])

@app.callback(
    Output('my-toggle-switch-output', 'children'),
    Input('cell-type-switch', 'value')
)
# def update_output(value):
#     if value == False:
#         return f'No'
#     elif value == True:
#         return f'Yes'


#upload aditional data        
@du.callback(
    output=Output('status2', 'children'),
    id='upload-data-addition-spot',
)
def upload_spot_data(filenames):
    thor,args=get_parameter()
    sample_id = args['sampleId']
    h5ad_converter(data_path,sample_id,h5ad_spot=filenames[0])
    shutil.rmtree("../data_input_temp/tmp/")
    return None 

@du.callback(
    output=Output('status3', 'children'),
    id='upload-data-addition-cell',
)
def upload_cell_data(filenames):
    thor,args=get_parameter()
    sample_id = args['sampleId']
    h5ad_converter(data_path,sample_id,h5ad_cell=filenames[0])
    shutil.rmtree("../data_input_temp/tmp/")
    return None 



    

# choose type of visualization
@app.callback(
    Output('gene-dropdown', 'children'),
    Input('spot-cell-option', 'value'),
    Input('visual-type-container', 'value')
)
def update_output_visual(value1,value2):
    if value2 == "Gene Expression":
        sample_id=session_id
        if value1 == "Spot data":
            files=files_generate(sample_id)
            gene_name=pd.read_csv(os.path.join(data_path,files["spots-gene-names"]),header=None,index_col=0)
            gene_list=list(gene_name.index)
            
        elif value1 == "Cell data":
            files=files_generate(sample_id)
            gene_name=pd.read_csv(os.path.join(data_path,files["cells-gene-names"]),header=None,index_col=0)
            gene_list=list(gene_name.index)
        return html.Div(children=[
            html.H5("Choose gene"),
            dcc.Dropdown(gene_list, id="gene-input-container",className='dropdown-input',placeholder="Select Gene")
            ])
    else:
        return None


# choose gene
@app.callback(
    Output('output-image', 'children'),
    Input('spot-cell-option', 'value'),
    Input('gene-input-container', 'value')
)
def get_gene(value1,value2):

    if value1 == "Spot data":
        if value2 is not None:
            with open('../user/args.json') as f:
                args = json.load(f)
            args['selectedSpotGeneName']=value2
            with open('../user/args.json', 'w') as f:
                json.dump(args, f)
            print("here")
            calculation_spot()
            output_image = visualization_img_spot(data_path,cache_path)
            output_content = html.Div(
            id="output-image",
            children=[
                # TODO: Add something here
                output_image
            ]
             )
            return output_content
    elif value1 == "Cell data":
        if value2 is not None:
            with open('../user/args.json') as f:
                args = json.load(f)
            args['selectedCellGeneName']=value2
            with open('../user/args.json', 'w') as f:
                json.dump(args, f)
            calculation_cell()
            output_image = visualization_img_cell(data_path,cache_path,cell_type=True)
            output_content = html.Div(
            id="output-image",
            children=[
                # TODO: Add something here
                output_image
            ]
             )
            return output_content
        

# @app.callback(
#     Output('output-image', 'children'),
#     Input("analyze-button", "n_clicks"),
#     Input('spot-cell-option', 'value'),

#     prevent_initial_call=True
# )
# def start_visualize(n_clicks,value):
#     if n_clicks:
#         print(value)
#         print("start visualize")
#         if value == "Spot data":
#             output_image = visualization_img_spot(data_path,cache_path)
#         if value == "Cell data":
#             output_image = visualization_img_cell(data_path,cache_path,cell_type=True)

#         output_content = html.Div(
#             id="output-image",
#             children=[
#                 # TODO: Add something here
#                 output_image
#             ]
#         )

#         return output_content
#     return dash.no_update, dash.no_update




@app.callback(Output("geojson", "data"), Input("editControl", "geojson"))
def mirror(x):
    return x




# Run app
if __name__ == "__main__":
    app.run_server(host="localhost", port=8081,debug=True,dev_tools_hot_reload=False)