
import sys
sys.path.append("./")
from niceview.utils.convert import h5ad_converter
from niceview.utils.dataset import ThorQuery
from interface.biogis_inter.interface import *
import toml
import shutil
import json
from dash import html
from dash import dcc
import dash_uploader as du
import pandas as pd
import os
import time
import dash
from niceview.utils.tools import save_roi_data_img
import plotly.graph_objects as go
import numpy as np
from shapely.geometry import Point, Polygon

# get info
config = toml.load('../user/config.toml')
data_path = config['path']['data']
cache_path = config['path']['cache']
max_file_size = config['constant']['max_file_size']


# upload HE image
def upload_image(filenames_upload_image, folder_id):
    """
    Uploads the HE image and copy it to data path, then create client.

    Parameters:
        filenames_upload_image (list): List of uploaded filenames.

    Returns:
        None
    """

    # get thor parameter and all user input info
    thor, args, p_input_json = get_parameter(folder_id)

    # get sample id base on uuid from dash-uploader
    sample_id = os.path.split(os.path.split(filenames_upload_image[0])[0])[1]

    # get basename of upload image
    basename = os.path.splitext(os.path.basename(filenames_upload_image[0]))[0]

    # wait and chech if file exist in folder
    while not os.path.exists(filenames_upload_image[0]):
        time.sleep(5)

    # get the height and width of image to calculate max dimension
    height, width = thor.process_data(sample_id, img_path=filenames_upload_image[0])
    max_dim = max(height, width)

    # update it in user argument file
    args["heightWidth"] = [height, width]

    # check if max dimension > 10000, if yes put "temp" in their sample id, then update in arg file
    if max_dim > 10000:
        args['sampleId'] = "temp" + "-" + sample_id
        args['fileName'] = basename
        dumpjson_parameter_from_user_input(folder_id, args=args)

    # else keep the originale sample id, then update in arg file
    else:
        args['sampleId'] = sample_id
        args['fileName'] = basename
        dumpjson_parameter_from_user_input(folder_id, args=args)

    # get the name rules for further calculation, then copy it to data path
        files = files_generate(sample_id)
        shutil.copy(filenames_upload_image[0], os.path.join(data_path, files["img"]))
    # remove temp folder that dash-uploader created
    shutil.rmtree(f"../data_input_temp/tmp/{folder_id}")

    # calculate wsi of image
    get_wsi(folder_id)

    # update the scale factor of data
    thor, args, p_input_json = get_parameter(folder_id)
    sample_id_file = args["sampleIdFile"]
    cache = cache_generate(sample_id, sample_id_file=sample_id_file)
    gis_img_path = os.path.join(cache_path, cache["gis-img-file"])
    gis_img_path = os.path.abspath(gis_img_path)
    factor = thor.get_factor(gis_img_path)
    new_factor = f"e*{factor}"
    update_javascript(new_factor)
    return html.H5("Click the home button on the map to see input image", className="text") 


# choose cell or spot data
def show_cell_spot_upload(spot_cell_option, folder_id):
    """
    Displays the upload options based on the selected data type.

    Parameters:
        spot_cell_option (str): Selected data type ('Spot data' or 'Cell data').

    Returns:
        html.Div: Div containing upload options.
    """
    if spot_cell_option == "Spot data":
        # show spot upload box
        return html.Div(className="upload-data", children=[
            html.Br(), html.Br(),
            html.H5("Upload mask(.npz file, the file have to contain 'mask' in name) and gene expression data(.h5ad file) for spot:", className="text"),
            du.Upload(id='upload-data-addition-spot', max_file_size=10000)
        ])
    elif spot_cell_option == "Cell data":
        # show cell upload box
        return html.Div(className="upload-data", children=[
            html.Br(), html.Br(),
            html.H5("Upload mask(.npz file, the file have to contain 'mask' in name) and gene expression data(.h5ad file) for cell:", className="text"),
            du.Upload(id='upload-data-addition-cell', max_file_size=10000)
        ])
    else:
        # show nothing
        return html.Div(children=[])


# upload aditional data        
def upload_spot_data(filenames_upload_spot_data, folder_id):
    """
    Uploads additional spot data and performs necessary operations.

    Parameters:
        filenames_upload_spot_data (list): List of uploaded filenames.

    Returns:
        None
    """
    # get parameter in argument file
    thor, args, p_input_json = get_parameter(folder_id)
    sample_id = args['sampleId']
    files = files_generate(sample_id)

    # check if 'mask' in the file name, if yes copy to data path and rename it
    if "mask" in filenames_upload_spot_data[0]:
        while not os.path.exists(filenames_upload_spot_data[0]):
            time.sleep(5)
        shutil.copy(filenames_upload_spot_data[0], os.path.join(data_path, files["mask"]))
    # if it not mask file h5ad convert
    else:
        while not os.path.exists(filenames_upload_spot_data[0]):
            time.sleep(5)
        h5ad_converter(data_path, sample_id, h5ad_spot=filenames_upload_spot_data[0])

    # remove temp folder that dash-uploader created
    shutil.rmtree(f"../data_input_temp/tmp/{folder_id}")
    return None 


def upload_cell_data(filenames_upload_cell_data, folder_id):
    """
    Uploads additional cell data and performs necessary operations.

    Parameters:
        filenames_upload_cell_data (list): List of uploaded filenames.

    Returns:
        None
    """
    # get parameter in argument file, and get max dimension
    thor, args, p_input_json = get_parameter(folder_id)
    sample_id = args['sampleId']
    height = args["heightWidth"][0]
    width = args["heightWidth"][1]
    height_width = args["heightWidth"]
    max_dim = max(height, width)
    files = files_generate(sample_id)

    # check if 'mask' in the file name, if yes copy to data path and rename it
    if "mask" in filenames_upload_cell_data[0]:
        while not os.path.exists(filenames_upload_cell_data[0]):
            time.sleep(5)
    # check if max dim > 10000 reduce size
        if max_dim > 10000:
            max_dim = thor.process_data(sample_id, height_width=height_width, mask_path=filenames_upload_cell_data[0])
        else:
            shutil.copy(filenames_upload_cell_data[0], os.path.join(data_path, files["mask"]))
    
    # same thing with h5ad, then convert
    else:
        while not os.path.exists(filenames_upload_cell_data[0]):
            time.sleep(5)
        if max_dim > 10000:
            max_dim = thor.process_data(sample_id, height_width=height_width, adata_path=filenames_upload_cell_data[0])
            h5ad_converter(data_path, sample_id, h5ad_cell=os.path.join(data_path, files["cell-temp-h5ad"]))
        else:
            h5ad_converter(data_path, sample_id, h5ad_cell=filenames_upload_cell_data[0])

    # remove temp folder that dash-uploader created         
    shutil.rmtree(f"../data_input_temp/tmp/{folder_id}")

    # check and remind user to upload mask file
    while not os.path.exists(os.path.join(data_path, files["mask"])):
        time.sleep(5)
        return html.H5("Remember to upload mask data", className="text")
    return None 
    

# choose type of visualization
def update_output_visual(spot_cell_option, visualize_option, folder_id):
    """
    Updates the visualization options based on selected data type and visualization type.

    Parameters:
        spot_cell_option (str): Selected data type ('Spot data' or 'Cell data').
        visualize_option (str): Selected visualization type.

    Returns:
        html.Div: Div containing updated visualization options.
    """
    # if choose gene expression
    if visualize_option == "Gene Expression":

        # get parameter
        thor, args, p_input_json = get_parameter(folder_id)
        sample_id = args['sampleId']

        # if spot data
        if spot_cell_option == "Spot data":

            # get name follow the rule
            files = files_generate(sample_id)

            # check if file exist
            while not os.path.exists(os.path.join(data_path, files["spots-gene-names"])):
                time.sleep(5)
                return html.H5("Data has not successfuly uploaded yet", className="text")
            
            # read gene name file and convert to list
            gene_name = pd.read_csv(os.path.join(data_path, files["spots-gene-names"]), header=None, index_col=0)
            gene_list = list(gene_name.index)
        
        # if cell data
        elif spot_cell_option == "Cell data":

            # get name follow the rule
            files = files_generate(sample_id)

            # check file exist
            while not os.path.exists(os.path.join(data_path, files["cells-gene-names"])):
                time.sleep(5)
                return html.H5("Data has not successfuly uploaded yet", className="text")
            
            # read gene name file and convert to list 
            gene_name = pd.read_csv(os.path.join(data_path, files["cells-gene-names"]), header=None, index_col=0)
            gene_list = list(gene_name.index)

        # return the dropdown contain gene list
        return html.Div(children=[
            html.H5("Choose Gene", className="text"),
            dcc.Dropdown(gene_list, id="gene-input-container", className='dropdown-input', placeholder="Input or Select Gene")])
    
    # if chose Pathway
    elif visualize_option == "Pathway Enrichment Analysis":

        # if cell data
        if spot_cell_option == "Cell data":

            # return upload container
            return html.Div(children=[
                html.Div(className="upload-data", children=[
                    html.H5("Upload pathway", className="text"),
                    du.Upload(id='upload-data-pathway', max_file_size=10000)
                ])
            ])
        # or show not support
        else:
            return html.H5("Only support Cell data", className="text")
        
    # if chose CNV    
    elif visualize_option == "CNV":

        # if cell data
        if spot_cell_option == "Cell data":

            # return upload container
            return html.Div(children=[
                html.Div(className="upload-data", children=[
                    html.H5("Upload CNV", className="text"),
                    du.Upload(id='upload-data-cnv', max_file_size=10000)
                ])
            ])
        
        # or show not support
        else:
            return html.H5("Only support Cell data", className="text")
        
    # if chose cell detection   
    elif visualize_option == 'Cell Detection Check':

        # if cell data
        if spot_cell_option == "Cell data":

            # get parameter
            thor, args, p_input_json = get_parameter(folder_id)
            sample_id = args['sampleId']
            sample_id_file = args['sampleIdFile']

            # calculation 
            thor.cell_gis(sample_id)

            # get name by rule and copy to cache folder
            cache = cache_generate(sample_id, sample_id_file=sample_id_file)
            shutil.copy(os.path.join(cache_path, cache["gis-blend-cell-random-img"]), os.path.join(cache_path, cache["gis-blend-cell-random-img-file"]))

            # return button
            return html.Button('Show Cell Detection', className="button", id="cell-detection", n_clicks=0)
        
        # or show not support
        else:
            return html.H5("Only support Cell data", className="text") 
    

# cell detection button
def show_cell_detection(n_clicks, folder_id):

    # if click
    if n_clicks:

        # generate map
        map_input = visualization_img_cell(folder_id, data_path, cache_path, cell_detect=True)

        # visualize result
        return html.Div(id="input-image", children=[map_input])
    
    # if no click keep previous map
    else:
        raise dash.exceptions.PreventUpdate


# upload pathway
def upload_pathway(filenames_upload_pathway, folder_id):
    """
    Uploads pathway data and performs necessary operations.

    Parameters:
        filenames_upload_pathway (list): List of uploaded filenames.

    Returns:
        None
    """
    # get parameter
    thor, args, p_input_json = get_parameter(folder_id)
    sample_id = args['sampleId']

    # remove previous selected area folder, if not found pass
    try:
        os.remove(f"../user{folder_id}/selected_area.zip")
    except FileNotFoundError:
        pass

    # check if file exist
    while not os.path.exists(filenames_upload_pathway[0]):
        time.sleep(5)

    # convert h5ad to processed file    
    h5ad_converter(data_path, sample_id, h5ad_cell_pathway=filenames_upload_pathway[0])

    # remove temp folder create by dash uploader
    shutil.rmtree(f"../data_input_temp/tmp/{folder_id}")

    # get file name by rules
    files = files_generate(sample_id)

    # read pathway file convert to list
    pathway_name = pd.read_csv(os.path.join(data_path, files["cell-pathway-name"]), header=None, index_col=0)
    pathway_list = list(pathway_name.index)

    # return dropdown with all pathway
    return html.Div(children=[
            html.H5("Choose Pathway", className="text"),
            dcc.Dropdown(pathway_list, id="pathway-input-container", className='dropdown-input', placeholder="Select Pathway")
    ])


# choose pathway
def get_pathway_output(spot_cell_option, pathway_value, folder_id):
    """
    Handles the selection of a pathway and performs necessary calculation.

    Parameters:
        spot_cell_option (str): Selected data type ('Spot data' or 'Cell data').
        pathway_value (str): Selected pathway name.

    Returns:
        html.Div: Div containing the updated visualization.
    """

    # if cell data
    if spot_cell_option == "Cell data":

        # if pathway been chosen
        if pathway_value is not None:

            # get parameter
            thor, args, p_input_json = get_parameter(folder_id)
            sample_id = args['sampleId']

            # remove previous selected area folder, if not found pass
            try:
                os.remove(f"../user{folder_id}/selected_area.zip")
            except FileNotFoundError:
                pass

            # get pathway name from user input and dump it to json file
            args["selectedPathway"] = pathway_value
            dumpjson_parameter_from_user_input(folder_id, args=args)

            # check if user already chose that pathway, if not calculate and dump it to json file, if yes show previous result
            if pathway_value not in p_input_json['Pathway']:

                calculation_pathway(folder_id)

                p_input_json['Pathway'].append(pathway_value)
                dumpjson_parameter_from_user_input(folder_id, p_input_json=p_input_json)

                map_input = visualization_img_cell(folder_id, data_path, cache_path, pathway=True)

            else:
                sample_id_pathway = sample_id + "-" + pathway_value
                args['sampleIdPathway'] = sample_id_pathway

                dumpjson_parameter_from_user_input(folder_id, args=args)
                map_input = visualization_img_cell(folder_id, data_path, cache_path, pathway=True)

            # visualize result    
            return html.Div(id="input-image", children=[map_input])
        
        # or keep previous map
        else:
            raise dash.exceptions.PreventUpdate


# upload cnv
def upload_cnv(filenames_upload_cnv, folder_id):
    """
    Uploads cnv data and performs necessary operations.

    Parameters:
        filenames_upload_cnv (list): List of uploaded filenames.

    Returns:
        None
    """

    # get parameter
    thor, args, p_input_json = get_parameter(folder_id)
    sample_id = args['sampleId']

    # remove previous selected area folder, if not found pass
    try:
        os.remove(f"../user{folder_id}/selected_area.zip")
    except FileNotFoundError:
        pass

    # get file name by rules
    files = files_generate(sample_id)

    # check if all file needed exist
    while not os.path.exists(os.path.join(data_path, files["cell-barcode"])):
        time.sleep(5)
    while not os.path.exists(os.path.join(data_path, files["cell-info"])):
        time.sleep(5)
    while not os.path.exists(os.path.join(data_path, files["mask"])):
        time.sleep(5)
    cell_order = pd.read_csv(os.path.join(data_path, files["cell-barcode"]), index_col=0, sep="\t", header=None)

    # check extension of cnv file, to get appropriate read
    if ".csv" in filenames_upload_cnv[0]:
        input_cnv = pd.read_csv(filenames_upload_cnv[0], index_col=0)
    else:
        input_cnv = pd.read_csv(filenames_upload_cnv[0], index_col=0, sep="\t")

    # reorder cell order
    input_cnv = input_cnv.reindex(list(cell_order.index))

    # open cell info and put cnv info under label
    cell_info = pd.read_csv(os.path.join(data_path, files["cell-info"]), index_col=0)
    cell_info["label"] = list(input_cnv.iloc[:, 0])
    cell_info.to_csv(os.path.join(data_path, files["cell-info"]))

    # calculate cnv
    calculation_CNV(folder_id, label_analysis=True)
    map_input = visualization_img_cell(folder_id, data_path, cache_path, cell_type=True)
    shutil.rmtree(f"../data_input_temp/tmp/{folder_id}")
    return html.Div(id="input-image", children=[map_input])


# choose gene
def get_gene(spot_cell_option, gene_chosen, folder_id):
    """
    Handles the selection of a gene and performs necessary calculation.

    Parameters:
        spot_cell_option (str): Selected data type ('Spot data' or 'Cell data').
        gene_chosen (str): Selected gene name.

    Returns:
        html.Div: Div containing the updated visualization.
    """
    if spot_cell_option == "Spot data":
        if gene_chosen is not None:

            # Get Thor, arguments, and previous input JSON
            thor, args, p_input_json = get_parameter(folder_id)
            sample_id = args['sampleId']

            # remove previous selected area folder, if not found pass
            try:
                os.remove(f"../user{folder_id}/selected_area.zip")
            except FileNotFoundError:
                pass
            
            # Set the selected spot gene name and dump to args json
            args['selectedSpotGeneName'] = gene_chosen
            dumpjson_parameter_from_user_input(folder_id, args=args)

            # If gene not present in previous input JSON
            if gene_chosen not in p_input_json['SpotGene']:

                # Perform calculation for Spot data
                calculation_spot(folder_id)

                # Update previous input JSON with new gene
                p_input_json['SpotGene'].append(gene_chosen)
                dumpjson_parameter_from_user_input(folder_id, p_input_json=p_input_json)

                # Generate visualization for Spot data
                map_input = visualization_img_spot(folder_id, data_path, cache_path)

            # If gene already present in previous input JSON
            else:

                # get file name with gene name in it to prevent dash leaflet cache
                sample_id_gene_spot = sample_id + "-" + gene_chosen
                args['sampleIdSpotGene'] = sample_id_gene_spot
                dumpjson_parameter_from_user_input(folder_id, args=args)

                # Generate visualization for Spot data
                map_input = visualization_img_spot(folder_id, data_path, cache_path)

            # Return updated visualization
            return html.Div(id="input-image", children=[map_input])
        
        # If gene is None, prevent update
        else:
            raise dash.exceptions.PreventUpdate
         
    elif spot_cell_option == "Cell data":
        if gene_chosen is not None:

            # Get Thor, arguments, and previous input JSON
            thor, args, p_input_json = get_parameter(folder_id)
            sample_id = args['sampleId']

            # Attempt to remove a file, if it exists
            try:
                os.remove(f"../user{folder_id}/selected_area.zip")
            except FileNotFoundError:
                pass
            
            # Set the selected cell gene name and dump to args json
            args['selectedCellGeneName'] = gene_chosen
            dumpjson_parameter_from_user_input(folder_id, args=args)

            # If gene not present in previous input JSON
            if gene_chosen not in p_input_json['CellGene']:

                # Perform calculation for cell data
                calculation_cell(folder_id, label_analysis=False)

                # dump gene name in previous input json
                p_input_json['CellGene'].append(gene_chosen)
                dumpjson_parameter_from_user_input(folder_id, p_input_json=p_input_json)

                # Generate visualization for cell data
                map_input = visualization_img_cell(folder_id, data_path, cache_path)
            else:

                # get file name with gene name in it to prevent dash leaflet cache
                sample_id_gene_cell = sample_id + "-" + gene_chosen
                args['sampleIdCellGene'] = sample_id_gene_cell
                dumpjson_parameter_from_user_input(folder_id, args=args)

                # Generate visualization for cell data
                map_input = visualization_img_cell(folder_id, data_path, cache_path)

            # Return updated visualization
            return html.Div(id="input-image", children=[map_input])
        
        # If gene is None, prevent update
        else:
            raise dash.exceptions.PreventUpdate
        
    # If spot_cell_option is None, prevent update
    else:
        raise dash.exceptions.PreventUpdate


def reset(n_clicks, spot_cell_option, visual_type, folder_id):
    """
    Reset image to original center

    Parameters:
        spot_cell_option (str): Selected data type ('Spot data' or 'Cell data').
        gene_chosen (str): Selected gene name.
        n_clicks (int): number of clicks

    Returns:
        html.Div: Div containing the updated visualization.
    """
    if n_clicks:
        if spot_cell_option == "Spot data":
            if visual_type == "Gene Expression":
                map_input = visualization_img_spot(folder_id, data_path, cache_path)
        else:
            if visual_type == "Gene Expression":
                map_input = visualization_img_cell(folder_id, data_path, cache_path)
            elif visual_type == "CNV":
                map_input = visualization_img_cell(folder_id, data_path, cache_path, cell_type=True)
            elif visual_type == "Pathway Enrichment Analysis":
                map_input = visualization_img_cell(folder_id, data_path, cache_path, pathway=True)
            elif visual_type == "Cell Detection Check":
                map_input = visualization_img_cell(folder_id, data_path, cache_path, cell_detect=True)
            else:
                map_input = visualization_img_input(folder_id, data_path, cache_path)
    return html.Div(id="input-image", children=[map_input])


def copy_and_rename_file(n_clicks, folder_id):
    """Copy and rename file.
    
    Args:
        n_clicks (int): Number of clicks.
    
    Returns:
        Download zip folder
    """
    if n_clicks:
        thor, args, p_input_json = get_parameter(folder_id)

        try:
            os.remove(f"../user{folder_id}/selected_area.zip")
        except FileNotFoundError:
            pass
        sample_id = args['sampleId']
        dir_path = f'../user{folder_id}/selected_area/'
        os.makedirs(dir_path, exist_ok=True)
        coord_path = os.path.join(dir_path, 'coords.json')
        shutil.copyfile(f'../user{folder_id}/coords.json', coord_path)
        with open(coord_path, 'r') as f:
            coords = json.load(f)
        cell_adata, wsi_img = thor.get_cell_adata_and_img(sample_id)
        
        save_roi_data_img(coords, cell_adata, wsi_img, dir_path)

        drawn_geojson = []
        coords = []
        geojson_name = f'../user{folder_id}/roi.json'
        with open(geojson_name, 'w') as f:
            json.dump(drawn_geojson, f)
        coordjson_name = f'../user{folder_id}/coords.json'
        with open(coordjson_name, 'w') as f:
            json.dump(coords, f)
        zip_folder(dir_path, f"../user{folder_id}/selected_area.zip")
    try:
        return dcc.send_file(f"../user{folder_id}/selected_area.zip")
    except FileNotFoundError:
        return None


def save_roi(drawn_geojson, folder_id):
    """Save coordinates.
    
    Args:
        drawn_geojson (dict): GeoJSON.

    Returns:
        json file.
    """
    thor, args, p_input_json = get_parameter(folder_id)

    sample_id_file = args['sampleIdFile']
    mapper = thor.get_coord_mapping(sample_id_file)
    coords = []
    for region in drawn_geojson['features']:
        temp = region['geometry']['coordinates'][0]
        temp = [mapper(*point) for point in temp]
        temp = [[point[1], point[0]] for point in temp]  # y, x -> x, y
        coords.append(temp)
    
    geojson_name = f'../user{folder_id}/roi.json'
    with open(geojson_name, 'w') as f:
        json.dump(drawn_geojson, f)
    
    coordjson_name = f'../user{folder_id}/coords.json'
    with open(coordjson_name, 'w') as f:
        json.dump(coords, f)
    
    return None


def plot_stats(n_clicks, drawn_geojson, idx, folder_id):
    """Plot stats.
    
    Args:
        n_clicks (int): number of clicks.
        drawn_geojson (dict): GeoJSON.
        idx (int): index.
    
    Returns:
        fig: Figure.
    """
    if idx is None:
        idx = 0
    coords = []
    thor, args, p_input_json = get_parameter(folder_id)
    sample_id_file = args['sampleIdFile']
    mapper = thor.get_coord_mapping(sample_id_file)
    sample_id = args['sampleId']
    cell_adata, wsi_img = thor.get_cell_adata_and_img(sample_id)
    for region in drawn_geojson['features']:
        temp = region['geometry']['coordinates'][0]
        temp = [mapper(*point) for point in temp]
        temp = [[point[1], point[0]] for point in temp]  # y, x -> x, y
        coords.append(temp)
    
    target = np.array([])
    for coord in coords:
        roi = Polygon(coord)
        locs = list(map(lambda x: roi.contains(Point(x)), cell_adata.obsm['spatial']))
        to_keep = cell_adata[locs].copy()
        try:
            array = to_keep.X[:, idx].ravel()
        except AttributeError:
            array = to_keep.X[:, idx].toarray().ravel()
        target = np.concatenate((target, array))
    hist = go.Figure(
        data=go.Histogram(
            x=target,
        ),
        layout=go.Layout(
            title='Histogram',
            xaxis={'title': 'Gene expression'},
            yaxis={'title': 'Number of cells'},
        ),
    )
    table = go.Figure(
        data=go.Table(
            header={
                'values': ['Mean', 'Median', 'Std'],
                'align': 'center',
            },
            cells={
                'values': [
                    np.round(np.mean(target), 2),
                    np.round(np.median(target), 2),
                    np.round(np.std(target), 2),
                ],
                'align': 'center',
            },
        ),
        layout=go.Layout(
            title='Statistics',
        ),
    )
    
    return dcc.Graph(figure=hist), dcc.Graph(figure=table)
