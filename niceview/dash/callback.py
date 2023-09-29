import sys
sys.path.append("../")
from niceview.utils.convert import h5ad_converter
from niceview.utils.dataset import ThorQuery
from niceview.dash.interface import *
import toml
import shutil
import json
from dash import html
from dash import dcc
import dash_uploader as du
import pandas as pd
import os
import time
from niceview.utils.tools import save_roi_data_img


# session_id = "hello-rep"

# get info
config = toml.load('../user/config.toml')
data_path = config['path']['data']
cache_path = config['path']['cache']
max_file_size = config['constant']['max_file_size']


# upload HE image
def upload_image(filenames_upload_image):
    """
    Uploads the HE image and copy it to data path, then create client.

    Parameters:
        filenames_upload_image (list): List of uploaded filenames.

    Returns:
        None
    """
    # sample_id="68c2d72e-58d1-11ee-8a1b-e450ebb96603"
    sample_id = os.path.split(os.path.split(filenames_upload_image[0])[0])[1]
    basename = os.path.splitext(os.path.basename(filenames_upload_image[0]))[0]
    thor, args = get_parameter()
    height, width = thor.process_data(sample_id, img_path=filenames_upload_image[0])
    max_dim = max(height, width)
    args["heightWidth"] = [height, width]

    if max_dim > 10000:
        args['sampleId'] = "temp" + "-" + sample_id
        args['fileName'] = basename
        with open('../user/args.json', 'w') as f:
            json.dump(args, f)
    else:
        args['sampleId'] = sample_id
        args['fileName'] = basename

        with open('../user/args.json', 'w') as f:
            json.dump(args, f)
        thor, args = get_parameter()
        sample_id = args['sampleId']
        files = files_generate(sample_id)
        while not os.path.exists(filenames_upload_image[0]):
            time.sleep(5)
        shutil.copy(filenames_upload_image[0], os.path.join(data_path, files["img"]))
    shutil.rmtree("../data_input_temp/tmp/")
    get_wsi()
    return None 


# choose cell or spot data
def show_cell_spot_upload(spot_cell_option):
    """
    Displays the upload options based on the selected data type.

    Parameters:
        spot_cell_option (str): Selected data type ('Spot data' or 'Cell data').

    Returns:
        html.Div: Div containing upload options.
    """
    if spot_cell_option == "Spot data":
        return html.Div(className="upload-data", children=[
            html.Br(), html.Br(),
            html.H5("Upload gene expression data(.h5ad file) for spot:", className="text"),
            du.Upload(id='upload-data-addition-spot', max_file_size=10000)
        ])
    elif spot_cell_option == "Cell data":
        return html.Div(className="upload-data", children=[
            html.Br(), html.Br(),
            html.H5("Upload mask(.npz file, the file have to contain 'mask' in name) and gene expression data(.h5ad file) for cell:", className="text"),
            du.Upload(id='upload-data-addition-cell', max_file_size=10000)
        ])
    else:
        return html.Div(children=[])


# upload aditional data        
def upload_spot_data(filenames_upload_spot_data):
    """
    Uploads additional spot data and performs necessary operations.

    Parameters:
        filenames_upload_spot_data (list): List of uploaded filenames.

    Returns:
        None
    """
    thor, args = get_parameter()
    sample_id = args['sampleId']
    files = files_generate(sample_id)
    if "mask" in filenames_upload_spot_data[0]:
        while not os.path.exists(filenames_upload_spot_data[0]):
            time.sleep(5)
        shutil.copy(filenames_upload_spot_data[0], os.path.join(data_path, files["mask"]))
    else:
        while not os.path.exists(filenames_upload_spot_data[0]):
            time.sleep(5)
        h5ad_converter(data_path, sample_id, h5ad_spot=filenames_upload_spot_data[0])
    shutil.rmtree("../data_input_temp/tmp/")
    return None 


def upload_cell_data(filenames_upload_cell_data):
    """
    Uploads additional cell data and performs necessary operations.

    Parameters:
        filenames_upload_cell_data (list): List of uploaded filenames.

    Returns:
        None
    """
    thor, args = get_parameter()
    sample_id = args['sampleId']
    height = args["heightWidth"][0]
    width = args["heightWidth"][1]
    height_width = args["heightWidth"]
    max_dim = max(height, width)
    files = files_generate(sample_id)
    if "mask" in filenames_upload_cell_data[0]:
        while not os.path.exists(filenames_upload_cell_data[0]):
            time.sleep(5)
        if max_dim > 10000:
            max_dim = thor.process_data(sample_id, height_width=height_width, mask_path=filenames_upload_cell_data[0])
        else:
            shutil.copy(filenames_upload_cell_data[0], os.path.join(data_path, files["mask"]))
    else:
        while not os.path.exists(filenames_upload_cell_data[0]):
            time.sleep(5)
        if max_dim > 10000:
            max_dim = thor.process_data(sample_id, height_width=height_width, adata_path=filenames_upload_cell_data[0])
            h5ad_converter(data_path, sample_id, h5ad_cell=os.path.join(data_path, files["cell-temp-h5ad"]))
        else:
            h5ad_converter(data_path, sample_id, h5ad_cell=filenames_upload_cell_data[0])
    shutil.rmtree("../data_input_temp/tmp/")
    return None 
    

# choose type of visualization
def update_output_visual(spot_cell_option, visualize_option):
    """
    Updates the visualization options based on selected data type and visualization type.

    Parameters:
        spot_cell_option (str): Selected data type ('Spot data' or 'Cell data').
        visualize_option (str): Selected visualization type.

    Returns:
        html.Div: Div containing updated visualization options.
    """
    if visualize_option == "Gene Expression":
        thor, args = get_parameter()
        sample_id = args['sampleId']
        if spot_cell_option == "Spot data":
            files = files_generate(sample_id)
            gene_name = pd.read_csv(os.path.join(data_path, files["spots-gene-names"]), header=None, index_col=0)
            gene_list = list(gene_name.index)
            
        elif spot_cell_option == "Cell data":
            files = files_generate(sample_id)
            gene_name = pd.read_csv(os.path.join(data_path, files["cells-gene-names"]), header=None, index_col=0)
            gene_list = list(gene_name.index)
        return html.Div(children=[
            html.H5("Choose Gene"),
            dcc.Dropdown(gene_list, id="gene-input-container", className='dropdown-input', placeholder="Select Gene")])
    elif visualize_option == "Pathway Enrichment Analysis":
        if spot_cell_option == "Cell data":
            return html.Div(children=[
                html.Div(className="upload-data", children=[
                    html.H5("Upload pathway", className="text"),
                    du.Upload(id='upload-data-pathway', max_file_size=10000)
                ])
            ])
        else:
            return html.H5("Only support Cell data", className="text")
        
    elif visualize_option == "CNV":
        if spot_cell_option == "Cell data":
            return html.Div(children=[
                html.Div(className="upload-data", children=[
                    html.H5("Upload CNV", className="text"),
                    du.Upload(id='upload-data-cnv', max_file_size=10000)
                ])
            ])
        else:
            return html.H5("Only support Cell data", className="text")


# upload pathway
def upload_pathway(filenames_upload_pathway):
    """
    Uploads pathway data and performs necessary operations.

    Parameters:
        filenames_upload_pathway (list): List of uploaded filenames.

    Returns:
        None
    """
    try:
        os.remove("../user/selected_area.zip")
    except FileNotFoundError:
        pass
    thor, args = get_parameter()
    sample_id = args['sampleId']
    while not os.path.exists(filenames_upload_pathway[0]):
        time.sleep(5)
    h5ad_converter(data_path, sample_id, h5ad_cell_pathway=filenames_upload_pathway[0])
    shutil.rmtree("../data_input_temp/tmp/")
    files = files_generate(sample_id)
    pathway_name = pd.read_csv(os.path.join(data_path, files["cell-pathway-name"]), header=None, index_col=0)
    pathway_list = list(pathway_name.index)
    return html.Div(children=[
            html.H5("Choose Pathway", className="text"),
            dcc.Dropdown(pathway_list, id="pathway-input-container", className='dropdown-input', placeholder="Select Pathway")
    ])


# choose pathway
def get_pathway_output(spot_cell_option, pathway_value):
    """
    Handles the selection of a pathway and performs necessary calculation.

    Parameters:
        spot_cell_option (str): Selected data type ('Spot data' or 'Cell data').
        pathway_value (str): Selected pathway name.

    Returns:
        html.Div: Div containing the updated visualization.
    """
    if spot_cell_option == "Cell data":
        if pathway_value is not None:
            try:
                os.remove("../user/selected_area.zip")
            except FileNotFoundError:
                pass
            thor, args = get_parameter()
            sample_id = args['sampleId']
            args["selectedPathway"] = pathway_value
            with open('../user/args.json', 'w') as f:
                json.dump(args, f)

            with open('../user/previous-input.json', 'r') as p_input:
                p_input_json = json.load(p_input)
            if pathway_value not in p_input_json['Pathway']:

                calculation_pathway()

                p_input_json['Pathway'].append(pathway_value)
                with open('../user/previous-input.json', 'w') as p_input:
                    json.dump(p_input_json, p_input)

                map_input = visualization_img_cell(data_path, cache_path, pathway=True)

            else:
                sample_id_pathway = sample_id + "-" + pathway_value
                args['sampleIdPathway'] = sample_id_pathway

                with open('../user/args.json', 'w') as f:
                    json.dump(args, f)
                map_input = visualization_img_cell(data_path, cache_path, pathway=True)
            
    return html.Div(id="input-image", children=[map_input])


# upload cnv
def upload_cnv(filenames_upload_cnv):
    """
    Uploads cnv data and performs necessary operations.

    Parameters:
        filenames_upload_cnv (list): List of uploaded filenames.

    Returns:
        None
    """
    try:
        os.remove("../user/selected_area.zip")
    except FileNotFoundError:
        pass
    thor, args = get_parameter()
    sample_id = args['sampleId']
    files = files_generate(sample_id)
    cell_order = pd.read_csv(os.path.join(data_path, files["cell-barcode"]), index_col=0, sep="\t", header=None)
    if ".csv" in filenames_upload_cnv[0]:
        input_cnv = pd.read_csv(filenames_upload_cnv[0], index_col=0)
    else:
        input_cnv = pd.read_csv(filenames_upload_cnv[0], index_col=0, sep="\t")
    input_cnv = input_cnv.reindex(list(cell_order.index))
    cell_info = pd.read_csv(os.path.join(data_path, files["cell-info"]), index_col=0)
    cell_info["label"] = list(input_cnv["copykat.pred"])
    cell_info.to_csv(os.path.join(data_path, files["cell-info"]))
    calculation_CNV(label_analysis=True)
    map_input = visualization_img_cell(data_path, cache_path, cell_type=True)
    shutil.rmtree("../data_input_temp/tmp/")
    return html.Div(id="input-image", children=[map_input])


# choose gene
def get_gene(spot_cell_option, gene_chosen):
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
            try:
                os.remove("../user/selected_area.zip")
            except FileNotFoundError:
                pass
            thor, args = get_parameter()
            sample_id = args['sampleId']
            args['selectedSpotGeneName'] = gene_chosen
            with open('../user/args.json', 'w') as f:
                json.dump(args, f)
            
            with open('../user/previous-input.json', 'r') as p_input:
                p_input_json = json.load(p_input)
            if gene_chosen not in p_input_json['SpotGene']:

                calculation_spot()

                p_input_json['SpotGene'].append(gene_chosen)
                with open('../user/previous-input.json', 'w') as p_input:
                    json.dump(p_input_json, p_input)

                map_input = visualization_img_spot(data_path, cache_path)
                # map_input = visualization_img_all(data_path,cache_path)

            else:
                sample_id_gene_spot = sample_id + "-" + gene_chosen
                args['sampleIdSpotGene'] = sample_id_gene_spot

                with open('../user/args.json', 'w') as f:
                    json.dump(args, f)
                map_input = visualization_img_spot(data_path, cache_path)
            return html.Div(id="input-image", children=[map_input])
    elif spot_cell_option == "Cell data":
        if gene_chosen is not None:
            try:
                os.remove("../user/selected_area.zip")
            except FileNotFoundError:
                pass
            thor, args = get_parameter()
            sample_id = args['sampleId']
            args['selectedCellGeneName'] = gene_chosen
            file_name = args["sampleIdFile"]
            with open('../user/args.json', 'w') as f:
                json.dump(args, f)

            with open('../user/previous-input.json', 'r') as p_input:
                p_input_json = json.load(p_input)

            if gene_chosen not in p_input_json['CellGene']:
                if file_name not in p_input_json['CellGene']:
                    calculation_cell(label_analysis=False)
                else:
                    calculation_cell(label_analysis=False)

                p_input_json['CellGene'].append(gene_chosen)
                with open('../user/previous-input.json', 'w') as p_input:
                    json.dump(p_input_json, p_input)

                map_input = visualization_img_cell(data_path, cache_path)
                # map_input = visualization_img_all(data_path,cache_path)
            else:
                sample_id_gene_cell = sample_id + "-" + gene_chosen
                args['sampleIdCellGene'] = sample_id_gene_cell

                with open('../user/args.json', 'w') as f:
                    json.dump(args, f)
                map_input = visualization_img_cell(data_path, cache_path)
            return html.Div(id="input-image", children=[map_input])


def reset(n_clicks, spot_cell_option, visual_type):
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
                map_input = visualization_img_spot(data_path, cache_path)
        elif spot_cell_option == "Cell data":
            if visual_type == "Gene Expression":
                map_input = visualization_img_cell(data_path, cache_path)
            elif visual_type == "CNV":
                map_input = visualization_img_cell(data_path, cache_path, cell_type=True)
            elif visual_type == "Pathway Enrichment Analysis":
                map_input = visualization_img_cell(data_path, cache_path, pathway=True)
        else:
            map_input = visualization_img_input(data_path, cache_path)
    return html.Div(id="input-image", children=[map_input])


def copy_and_rename_file(n_clicks):
    """Copy and rename file.
    
    Args:
        n_clicks (int): Number of clicks.
    
    Returns:
        Download zip folder
    """
    if n_clicks:
        try:
            os.remove("../user/selected_area.zip")
        except FileNotFoundError:
            pass
        thor, args = get_parameter()
        sample_id = args['sampleId']
        dir_path = '../user/selected_area/'
        os.makedirs(dir_path, exist_ok=True)
        coord_path = os.path.join(dir_path, 'coords.json')
        shutil.copyfile('../user/coords.json', coord_path)
        with open(coord_path, 'r') as f:
            coords = json.load(f)
        cell_adata, wsi_img = thor.get_cell_adata_and_img(sample_id)
        
        save_roi_data_img(coords, cell_adata, wsi_img, dir_path)
        drawn_geojson = []
        coords = []
        geojson_name = '../user/roi.json'
        with open(geojson_name, 'w') as f:
            json.dump(drawn_geojson, f)
        coordjson_name = '../user/coords.json'
        with open(coordjson_name, 'w') as f:
            json.dump(coords, f)
        zip_folder(dir_path, "../user/selected_area.zip")

    return dcc.send_file("../user/selected_area.zip")


def save_roi(drawn_geojson):
    """Save coordinates.
    
    Args:
        drawn_geojson (dict): GeoJSON.

    Returns:
        json file.
    """
    thor, args = get_parameter()
    sample_id_file = args['sampleIdFile']
    mapper = thor.get_coord_mapping(sample_id_file)
    coords = []
    for region in drawn_geojson['features']:
        temp = region['geometry']['coordinates'][0]
        temp = [mapper(*point) for point in temp]
        temp = [[point[1], point[0]] for point in temp]  # y, x -> x, y
        coords.append(temp)
    
    geojson_name = '../user/roi.json'
    with open(geojson_name, 'w') as f:
        json.dump(drawn_geojson, f)
    
    coordjson_name = '../user/coords.json'
    with open(coordjson_name, 'w') as f:
        json.dump(coords, f)
    
    return None
