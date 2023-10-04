#python app.py --port=8082
import sys
sys.path.append("./")
from interface.biogis_inter.callback import *
import dash
from dash import html
from dash import dcc
from dash.dependencies import Input, Output,State
import dash_uploader as du
import toml
import dash_loading_spinners as dls
import shutil
import argparse
import uuid



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

folder_id = ""


# App layout
def app_layout():
    # new_factor = "e*10.094100024003849"
    # update_javascript(new_factor)
    try:
        os.remove("../user_/selected_area.zip")    
    except FileNotFoundError:
        pass
    try:
        shutil.rmtree("../user_/selected_area/")
    except FileNotFoundError:
        pass
    clear_cache(folder_id)
    clear_data(folder_id)
    dump_default_para_arg()
    return html.Div( id="body",
        
        children=[
            html.Link(rel='stylesheet', href="https://use.fontawesome.com/releases/v5.5.0/css/all.css", integrity="sha384-B4dIYHKNBt8Bc12p+WXckhzcICo0wtJAoU8YZTY5qE0Id1GSseTk6S+L3BlXeVIU", crossOrigin="anonymous"),
            html.Link(rel='stylesheet', href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/4.7.0/css/font-awesome.min.css"),
            html.Div(id='header',
            children=[
                html.H4(id="logo", children='BioGIS', className="text")
            ]),
            html.Div(
                className='container clearfix', id="display",
                children=[
                    html.Div(id='submit-container', children=[  
                            html.Div(className="upload-data", children=[
                                html.H5("Upload H&E image:", className="text"),
                                du.Upload(id='upload-data-image', max_file_size=10000),
                            ]),
                            html.Br(),
                            dls.Hash(
                                html.Div(id="status1"),
                                color="#ffffff",
                                speed_multiplier=2,
                                size=100,
                            ),
                            html.Br(),
                            html.Div(children=[
                                html.H5("Choose type of data", className="text"),
                                dcc.Dropdown(['Spot data', 'Cell data'], id="spot-cell-option", className='dropdown-input', placeholder="Choose Cell or Spot data type")
                                ]),
                            html.Br(),
                            html.Br(),
                            dls.Hash(
                                html.Div(id="additional-data-box"),
                                color="#ffffff",
                                speed_multiplier=2,
                                size=100,
                            ),
                            html.Br(), html.Br(),
                            html.H5("Choose type of visualization", className="text"),
                            html.Div(children=[dcc.Dropdown(['Pathway Enrichment Analysis', 'Gene Expression', 'CNV', 'Cell Detection Check'], id="visual-type-container", className='dropdown-input', placeholder="Select Type")]),
                            html.Br(), html.Br(), html.Br(),
                            dls.Hash(
                                html.Div(id="gene-dropdown"),
                                color="#ffffff",
                                speed_multiplier=2,
                                size=100,
                            ),
                            html.Br(),
                            html.Div(id="pathway-dropdown"),
                            # html.Button('Done(clear cache)', className="button", id="clear-cache", n_clicks=0),
                            html.Br(),html.Br(),html.Br(),html.Br(),
                            html.H5("Change number of columns of graph", className="text"), 
                            dcc.Input(
                                id='colIdx',
                                className="input-column-num",
                                type='number',
                                placeholder=1,
                            ),
                            html.Button('Plot', id='plotButton',className="button", n_clicks=0),
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
                                size=100
                            ),
                            ]),
                    html.Div(
                        id='left-column-temp',
                        children=[
                            html.H2("Visualize HE image:", className="text"),
                            dls.Hash(
                                html.Div(
                                    id="input-image",
                                    children=[
                                        visualization_img_input(folder_id, data_path, cache_path),
                                    ],
                                ),
                                color="#ffffff",
                                speed_multiplier=2,
                                size=100,
                            ),
                        ],
                    ),
                    html.Div(
                        id='right-column-temp',
                        children=[
                            html.Div(
                                id="output-image",
                                children=[
                                    # map_output,
                                ],
                            ),
                        ],
                    ),
                    html.Div(id ="graph", children=[
                        html.Div(
                            id='left-column',
                            children=[
                                dls.Hash(
                                    html.Div(
                                        id="hist"
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
                                dls.Hash(
                                    html.Div(
                                        id="stats"
                                    ),
                                    color="#ffffff",
                                    speed_multiplier=2,
                                    size=100,
                                ),
                            ],
                        ),
                        ]
                    )
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
    return upload_image(filenames_upload_image, folder_id)


# choose cell or spot data
@app.callback(
    Output('additional-data-box', 'children'),
    Input('spot-cell-option', 'value'),
)
def callback_show_cell_spot_upload(spot_cell_option):
    return show_cell_spot_upload(spot_cell_option, folder_id)


# upload aditional data 
@du.callback(
    output=Output('status2', 'children'),
    id='upload-data-addition-spot',
)
def callback_upload_spot_data(filenames_upload_spot_data):
    return upload_spot_data(filenames_upload_spot_data, folder_id)


@du.callback(
    output=Output('status3', 'children'),
    id='upload-data-addition-cell',
)
def callback_upload_cell_data(filenames_upload_cell_data):
    return upload_cell_data(filenames_upload_cell_data, folder_id)


# choose type of visualization
@app.callback(
    Output('gene-dropdown', 'children'),
    Input('spot-cell-option', 'value'),
    Input('visual-type-container', 'value')
)
def callback_update_output_visual(spot_cell_option, visualize_option):
    return update_output_visual(spot_cell_option, visualize_option, folder_id)

@app.callback(
    Output('input-image', 'children', allow_duplicate=True),
    Input("cell-detection", "n_clicks"),
    prevent_initial_call='initial_duplicate'
)
def callback_show_cell_detection(n_clicks):
    return show_cell_detection(n_clicks, folder_id)

# upload cnv file
@du.callback(
    output=Output('pathway-dropdown', 'children'),
    id="upload-data-pathway"
)
def callback_upload_pathway(filenames_upload_pathway):
    return upload_pathway(filenames_upload_pathway, folder_id)


@app.callback(
    Output('input-image', 'children', allow_duplicate=True),
    Input('spot-cell-option', 'value'),
    Input('pathway-input-container', 'value'),
    prevent_initial_call='initial_duplicate'
)
def callback_get_pathway_output(spot_cell_option, pathway_value):
    return get_pathway_output(spot_cell_option, pathway_value, folder_id)


# upload cnv file
@du.callback(
    output=Output('input-image', 'children'),
    id="upload-data-cnv"
)
def callback_upload_cnv(filenames_upload_cnv):
    return upload_cnv(filenames_upload_cnv, folder_id)


# choose gene
@app.callback(
    Output('input-image', 'children', allow_duplicate=True),
    Input('spot-cell-option', 'value'),
    Input('gene-input-container', 'value'),
    prevent_initial_call='initial_duplicate'
)
def callback_get_gene(spot_cell_option, gene_chosen):
    return get_gene(spot_cell_option, gene_chosen, folder_id)


@app.callback(
    Output('input-image', 'children', allow_duplicate=True),
    Input("btn_home", "n_clicks"),
    Input('spot-cell-option', 'value'),
    Input('visual-type-container', 'value'),    
    prevent_initial_call='initial_duplicate'
)
def callback_reset(n_clicks, spot_cell_option, visual_type):
    return reset(n_clicks, spot_cell_option, visual_type, folder_id)


@app.callback(
    Output('download', 'data'),
    Input('btn_save', 'n_clicks'),
    prevent_initial_call=True
)
def callback_copy_and_rename_file(n_clicks):
    return copy_and_rename_file(n_clicks, folder_id)


@app.callback(
    Output('status5', 'children'),
    Input('editControl', 'geojson'),
    prevent_initial_call=True,
)
def callback_save_roi(drawn_geojson):
    return save_roi(drawn_geojson, folder_id)

@app.callback(
    Output('hist', 'children'),
    Output('stats', 'children'),
    Input('plotButton', 'n_clicks'),
    Input('editControl', 'geojson'),
    State('colIdx', 'value'),
    prevent_initial_call=True,
)
def callback_plot_stats(n_clicks, drawn_geojson, idx):
    return plot_stats(n_clicks, drawn_geojson, idx, folder_id)


parser = argparse.ArgumentParser(description='Run Dash app.')
parser.add_argument('--port', type=int, default=8081, help='Port to run the app on')
args = parser.parse_args()
# Run app
if __name__ == "__main__":
    app.run_server(host="localhost", port=args.port, debug=True, dev_tools_hot_reload=False)
