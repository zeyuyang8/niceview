import sys
sys.path.append("../")
from niceview.dash.callback import *
from niceview.dash.interface import *
import dash
from dash import html
from dash import dcc
from dash.dependencies import Input, Output
import dash_uploader as du
import toml
import dash_loading_spinners as dls
import shutil

# get info
config = toml.load('../user/config.toml')
data_path = config['path']['data']
cache_path = config['path']['cache']
max_file_size = config['constant']['max_file_size']


# Initialize app
app = dash.Dash(
    __name__,
    meta_tags=[
        {"name": "viewport"}
    ],
    
)

app.title = "BioGIS"
server = app.server
du.configure_upload(app, r"../data_input_temp/tmp/", use_upload_id=True)


# App layout
def app_layout():
    try:
        shutil.rmtree("../data_input_temp/tmp/")    
    except FileNotFoundError:
        pass
    try:
        shutil.rmtree("../user/selected_area/")
    except FileNotFoundError:
        pass
    clear_cache()
    clear_data()
    get_default_para()
    return html.Div( id="body",
        
        children=[
            html.Link(rel='stylesheet', href="https://use.fontawesome.com/releases/v5.5.0/css/all.css", integrity="sha384-B4dIYHKNBt8Bc12p+WXckhzcICo0wtJAoU8YZTY5qE0Id1GSseTk6S+L3BlXeVIU", crossOrigin="anonymous"),
            html.Link(rel='stylesheet', href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/4.7.0/css/font-awesome.min.css"),
            html.Div(id='header',
            children=[
                html.H4(id="logo", children='BioGIS', className="text")
            ]),
            html.Div(
                className='container clearfix', id="graph",
                children=[
                    html.Div(id='submit-container', children=[  
                            html.Div(className="upload-data", children=[
                                html.H5("Upload H&E image:", className="text"),
                                du.Upload(id='upload-data-image', max_file_size=10000),
                            ]),
                            html.Br(),
                            html.Div(children=[
                                html.H5("Choose type of data", className="text"),
                                dcc.Dropdown(['Spot data', 'Cell data'], id="spot-cell-option", className='dropdown-input', placeholder="Choose Cell or Spot data type")
                                ]),
                            html.Br(),
                            html.Div(id="additional-data-box"),
                            html.Br(), html.Br(),
                            html.H5("Choose type of visualization", className="text"),
                            html.Div(children=[dcc.Dropdown(['Pathway Enrichment Analysis', 'Gene Expression', 'CNV'], id="visual-type-container", className='dropdown-input', placeholder="Select Type")]),
                            html.Br(), html.Br(), html.Br(),
                            html.Div(id="gene-dropdown"),
                             html.Br(),
                            html.Div(id="pathway-dropdown"),
                            # html.Button('Done(clear cache)', className="button", id="clear-cache", n_clicks=0),
                            dls.Hash(
                                html.Div(id="status1"),
                                color="#ffffff",
                                speed_multiplier=2,
                                size=100,
                            ),
                            dls.Hash(
                                html.Div(id="status2"),
                                color="#ffffff",
                                speed_multiplier=2,
                                size=100,
                            ),
                            dls.Hash(
                                html.Div(id="status3"),
                                color="#ffffff",
                                speed_multiplier=2,
                                size=100,
                            ),
                            dls.Hash(
                                dcc.Download(id="download"),
                                color="#ffffff",
                                speed_multiplier=2,
                                size=100,
                            ),
                            dls.Hash(
                                html.Div(id="status5"),
                                color="#ffffff",
                                speed_multiplier=2,
                                size=100,
                            ),
                            ]),
                    html.Div(
                        id='left-column',
                        children=[
                            dls.Hash(
                                html.Div(
                                    id="input-image",
                                    children=[
                                        visualization_img_input(data_path, cache_path),
                                    ],
                                ),
                                color="#ffffff",
                                speed_multiplier=2,
                                size=100,
                            ),
                        ],
                    ),
                    html.Div(
                        id='right-column',
                        children=[
                            html.Div(
                                id="output-image",
                                children=[
                                    # map_output,
                                ],
                            ),
                        ],
                    ),
                ],
            ),
        html.Div(
            id="footer", children=[
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


app.layout = app_layout


# upload HE image
@du.callback(
    output=Output('status1', 'children'),
    id='upload-data-image',
)
def callback_upload_image(filenames_upload_image):
    return upload_image(filenames_upload_image)


# choose cell or spot data
@app.callback(
    Output('additional-data-box', 'children'),
    Input('spot-cell-option', 'value')
)
def callback_show_cell_spot_upload(spot_cell_option):
    return show_cell_spot_upload(spot_cell_option)


# upload aditional data 
@du.callback(
    output=Output('status2', 'children'),
    id='upload-data-addition-spot',
)
def callback_upload_spot_data(filenames_upload_spot_data):
    return upload_spot_data(filenames_upload_spot_data)


@du.callback(
    output=Output('status3', 'children'),
    id='upload-data-addition-cell',
)
def callback_upload_cell_data(filenames_upload_cell_data):
    return upload_cell_data(filenames_upload_cell_data)


# choose type of visualization
@app.callback(
    Output('gene-dropdown', 'children'),
    Input('spot-cell-option', 'value'),
    Input('visual-type-container', 'value')
)
def callback_update_output_visual(spot_cell_option, visualize_option):
    return update_output_visual(spot_cell_option, visualize_option)


# upload cnv file
@du.callback(
    output=Output('pathway-dropdown', 'children'),
    id="upload-data-pathway"
)
def callback_upload_pathway(filenames_upload_pathway):
    return upload_pathway(filenames_upload_pathway)


@app.callback(
    Output('input-image', 'children', allow_duplicate=True),
    Input('spot-cell-option', 'value'),
    Input('pathway-input-container', 'value'),
    prevent_initial_call='initial_duplicate'
)
def callback_get_pathway_output(spot_cell_option, pathway_value):
    return get_pathway_output(spot_cell_option, pathway_value)


# upload cnv file
@du.callback(
    output=Output('input-image', 'children'),
    id="upload-data-cnv"
)
def callback_upload_cnv(filenames_upload_cnv):
    return upload_cnv(filenames_upload_cnv)


# choose gene
@app.callback(
    Output('input-image', 'children', allow_duplicate=True),
    Input('spot-cell-option', 'value'),
    Input('gene-input-container', 'value'),
    prevent_initial_call='initial_duplicate'
)
def callback_get_gene(spot_cell_option, gene_chosen):
    return get_gene(spot_cell_option, gene_chosen)


@app.callback(
    Output('input-image', 'children', allow_duplicate=True),
    Input("btn_home", "n_clicks"),
    Input('spot-cell-option', 'value'),
    Input('visual-type-container', 'value'),    
    prevent_initial_call='initial_duplicate'
)
def callback_reset(n_clicks, spot_cell_option, visual_type):
    return reset(n_clicks, spot_cell_option, visual_type)


@app.callback(
    Output('download', 'data'),
    Input('btn_save', 'n_clicks'),
    prevent_initial_call=True
)
def callback_copy_and_rename_file(n_clicks):
    return copy_and_rename_file(n_clicks)


@app.callback(
    Output('status5', 'children'),
    Input('editControl', 'geojson'),
    prevent_initial_call=True,
)
def callback_save_roi(drawn_geojson):
    return save_roi(drawn_geojson)

# Run app
if __name__ == "__main__":
    app.run_server(host="0.0.0.0", port=8081, debug=False, dev_tools_hot_reload=False)
