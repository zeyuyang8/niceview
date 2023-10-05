
import sys
sys.path.append("./")
from niceview.utils.dataset import ThorQuery
from niceview.pyplot.leaflet import create_leaflet_map
#from interface.biogis_inter.leaflet import *
import shutil
import json
import toml
import os
import zipfile
import re


def dump_default_para_config():
    """
    Sets the config parameters to their default values.

    Parameters:
        None

    Returns:
        None
    """
    with open('../user/config-default.toml', 'r') as conf:
        configs = toml.load(conf)
    with open('../user/config.toml', 'w') as conf:
        toml.dump(configs, conf)


def update_real_path_data_cache():
    """
    Get the realpath of data and cache folder

    Parameters:
        None

    Returns:
        None
    """
    with open('../user/config.toml', 'r') as conf:
        configs = toml.load(conf)
    configs['path']['data'] = os.path.abspath(configs['path']['data']) + "/"
    configs['path']['cache'] = os.path.abspath(configs['path']['cache']) + "/"
    with open('../user/config.toml', 'w') as conf:
        toml.dump(configs, conf)


dump_default_para_config()
update_real_path_data_cache()
configs = toml.load('../user/config.toml')
data_path = configs['path']['data']
cache_path = configs['path']['cache']
max_file_size = configs['constant']['max_file_size']


def dump_default_para_arg():
    """
    Sets the application parameters to their default values.

    Parameters:
        None

    Returns:
        None
    """
    with open('../user/args-default.json') as f:
        args_default = json.load(f)
    with open('../user/args.json', 'w') as f:
        json.dump(args_default, f)
    with open('../user/previous-input-default.json', 'r') as p_input:
        p_input_default = json.load(p_input)
    with open('../user/previous-input.json', 'w') as p_input:
        json.dump(p_input_default, p_input)


def dumpjson_parameter_from_user_input(folder_id, args=None, p_input_json=None):
    """
    Dump user-provided parameters and input JSON to corresponding files.

    Parameters:
        folder_id (str): Unique folder identifier for each sample.

        args (dict, optional): A dictionary containing user-provided parameters. Default is None.

        p_input_json (dict, optional): A dictionary containing previous input JSON data.
            Default is None.

    Returns:
        None
    """
    if args is not None:
        with open(f'../user{folder_id}/args.json', 'w') as f:
            json.dump(args, f)
    if p_input_json is not None:
        with open(f'../user{folder_id}/previous-input.json', 'w') as p_input:
            json.dump(p_input_json, p_input)


def files_generate(sample_id):
    """
    Generate file names for various data and cache files based on the sample_id.

    Parameters:
        sample_id (str): The sample ID.

    Returns:
        dict: A dictionary containing file names for cells, images, and spots.
    """
    files = {
        'cell-h5ad': '-'.join([sample_id, 'cell.h5ad']),
        'cell-temp-h5ad': '-'.join([sample_id, 'cell-temp.h5ad']),
        'cells-gene-names': '-'.join([sample_id, 'cell-gene-name.txt']),
        'cell-info': '-'.join([sample_id, 'cell-info.csv']),
        'cell-barcode': '-'.join([sample_id, 'cell-barcode.txt']),
        'cell-pathway-name': '-'.join([sample_id, 'cell-pathway-name.txt']),
        'img': '-'.join([sample_id, 'wsi-img.tiff']),
        'spots-gene-names': '-'.join([sample_id, 'spot-gene-name.txt']),
        'mask': '-'.join([sample_id, 'cell-mask.npz']),
        

    }
    return files


def cache_generate(sample_id, sample_id_file='', sample_id_gene_cell='', sample_id_gene_spot='', sample_id_pathway=''):
    """
    Generate cache file names based on sample IDs and file names.

    Parameters:
        sample_id (str): The sample ID.
        sample_id_file (str): The sample ID with the file name.
        sample_id_gene_cell (str): The sample ID with the cell gene name.
        sample_id_gene_spot (str): The sample ID with the spot gene name.

    Returns:
        dict: A dictionary containing cache file names.
    """
    cache = {
        'gis-img': '-'.join([sample_id, 'gis-wsi-img.tiff']),
        'gis-blend-cells': '-'.join([sample_id, 'gis-blend-cell-gene-img.tiff']),
        'gis-blend-spots': '-'.join([sample_id, 'gis-blend-spot-gene-img.tiff']),
        'gis-blend-cell-type': '-'.join([sample_id, 'gis-blend-cell-type-img.tiff']),
        'gis-blend-cell-heatmap': '-'.join([sample_id, 'gis-blend-cell-pathway-heatmap-img.tiff']),
        'gis-blend-cell-random-img': '-'.join([sample_id, 'gis-blend-cell-random-img.tiff']),

        'gis-img-file': '-'.join([sample_id_file, 'gis-wsi-img.tiff']),
        'gis-blend-cells-gene': '-'.join([sample_id_gene_cell, 'gis-blend-cell-gene-img.tiff']),
        'gis-blend-spots-gene': '-'.join([sample_id_gene_spot, 'gis-blend-spot-gene-img.tiff']),
        'gis-blend-cell-type-file': '-'.join([sample_id_file, 'gis-blend-cell-type-img.tiff']),
        'gis-blend-cell-pathway-heatmap': '-'.join([sample_id_pathway, 'gis-blend-cell-pathway-heatmap-img.tiff']),
        'gis-blend-cell-random-img-file': '-'.join([sample_id_file, 'gis-blend-cell-random-img.tiff']),
    }
    return cache


def get_parameter(folder_id):
    """
    Get parameters from configuration files and create a ThorQuery object.

    Parameters:
        folder_id: unique folder id for each sample

    Returns:
        tuple: A tuple containing ThorQuery object and parameters.
    """
    with open('../db/db-info.json', 'r') as json_file:
        db_info = json.load(json_file)
    data_extension = db_info['data_extension']
    cache_extension = db_info['cache_extension']
    cell_label_encoder = db_info['cell_label_encoder']
    cell_label_cmap = db_info['cell_label_cmap']
    primary_key_list = db_info['primary_key_list']
    with open(f'../user{folder_id}/args.json') as f:
        args = json.load(f)
    with open(f'../user{folder_id}/previous-input.json', 'r') as p_input:
        p_input_json = json.load(p_input)
    thor = ThorQuery(
        data_path,
        cache_path,
        data_extension,
        cache_extension,
        cell_label_encoder,
        cell_label_cmap,
        primary_key_list,
    )
    
    return thor, args, p_input_json


def update_javascript(new_factor):
    """
    Update a specific line in a JavaScript file with a new mathematical expression.

    Parameters:
        new_factor (str): The new mathematical expression to replace the existing one.
            This should be a string representing a JavaScript expression. For example,
            'e*2.5' will replace the existing expression with 'Math.round(e*2.5)'.

    Returns:
        None
    """
    file_path = "./assets/dash_leaflet.js"
    # Read the file content
    with open(file_path, 'r') as file:
        content = file.read()
    # Define the pattern to match the line you want to replace
    pattern = r'_updateMetric:function\(t\){var e=this\._getRoundNum\(t\),n=Math\.round\(e\*[-+]?\d*\.\d+\)'
    # Define the replacement string
    replacement = r'_updateMetric:function(t){var e=this._getRoundNum(t),n=Math.round(' + new_factor + ')'

    # Perform the replacement
    new_content = re.sub(pattern, replacement, content)

    # Save the modified content back to the file
    with open(file_path, 'w') as file:
        file.write(new_content)


def get_wsi(folder_id):
    """
    Get client for wsi image and perform caching.

    Parameters:
        folder_id: unique folder id for each sample

    Returns:
        None
    """
    thor, args, p_input_json = get_parameter(folder_id)
    sample_id = args['sampleId']
    sample_id_file = args['sampleId'] + '-' + args['fileName']

    args["sampleIdFile"] = sample_id_file

    with open(f'../user{folder_id}/args.json', 'w') as f:
        json.dump(args, f)
    thor.wsi_gis(sample_id)
    cache = cache_generate(sample_id, sample_id_file=sample_id_file)
    shutil.copy(os.path.join(cache_path, cache["gis-img"]), os.path.join(cache_path, cache["gis-img-file"]))
    

def calculation_cell(folder_id, label_analysis=True):
    """
    Perform cell gene analysis and caching.

    Parameters:
        folder_id: unique folder id for each sample
        label_analysis (bool): Whether to perform cell type analysis.

    Returns:
        None
    """
    thor, args, p_input_json = get_parameter(folder_id)
    sample_id = args['sampleId']
    sample_id_file = args['sampleIdFile']
    selected_cell_gene_name = args['selectedCellGeneName']
    sample_id_gene_cell = sample_id + "-" + selected_cell_gene_name

    thor.empty_cache_cell(sample_id, gene=True, label=True)
    
    args['sampleIdCellGene'] = sample_id_gene_cell

    with open(f'../user{folder_id}/args.json', 'w') as f:
        json.dump(args, f)

    thor.cell_gis(
        sample_id,
        selected_cell_gene_name,
        label_analysis=label_analysis,
    )
    
    cache = cache_generate(sample_id, sample_id_gene_cell=sample_id_gene_cell, sample_id_file=sample_id_file)
    shutil.copy(os.path.join(cache_path, cache["gis-blend-cells"]), os.path.join(cache_path, cache["gis-blend-cells-gene"]))
    # shutil.copy(os.path.join(cache_path, cache["gis-blend-cell-type"]), os.path.join(cache_path, cache["gis-blend-cell-type-file"]))


def calculation_pathway(folder_id, label_analysis=False):
    """
    Perform pathway analysis and caching.

    Parameters:
        folder_id: unique folder id for each sample
        label_analysis (bool): Whether to perform cell type analysis.

    Returns:
        None
    """
    thor, args, p_input_json = get_parameter(folder_id)
    sample_id = args['sampleId']
    selected_pathway = args["selectedPathway"]
    
    sample_id_pathway = sample_id + "-" + selected_pathway

    thor.empty_cache_cell(sample_id, pathway=True)
    
    args['sampleIdPathway'] = sample_id_pathway

    with open(f'../user{folder_id}/args.json', 'w') as f:
        json.dump(args, f)

    thor.cell_gis(
        sample_id=sample_id,
        selected_pathway=selected_pathway,
        label_analysis=label_analysis,
    )
    
    cache = cache_generate(sample_id, sample_id_pathway=sample_id_pathway)
    shutil.copy(os.path.join(cache_path, cache["gis-blend-cell-heatmap"]), os.path.join(cache_path, cache["gis-blend-cell-pathway-heatmap"]))


def calculation_CNV(folder_id, label_analysis=True):
    """
    Perform CNV analysis and caching.

    Parameters:
        folder_id: unique folder id for each sample
        label_analysis (bool): Whether to perform cell type analysis.

    Returns:
        None
    """
    thor, args, p_input_json = get_parameter(folder_id)
    sample_id = args['sampleId']
    sample_id_file = args['sampleIdFile']

    with open(f'../user{folder_id}/previous-input.json', 'r') as p_input:
        p_input_json = json.load(p_input)

    if sample_id_file not in p_input_json["CellType"]:
        p_input_json["CellType"].append(sample_id_file)

    with open(f'../user{folder_id}/previous-input.json', 'w') as p_input:
        json.dump(p_input_json, p_input)

    thor.empty_cache_cell(sample_id, label=True)

    thor.cell_gis(
        sample_id=sample_id,
        label_analysis=label_analysis,
    )
    
    cache = cache_generate(sample_id, sample_id_file=sample_id_file)
    shutil.copy(os.path.join(cache_path, cache["gis-blend-cell-type"]), os.path.join(cache_path, cache["gis-blend-cell-type-file"]))


def calculation_spot(folder_id):
    """
    Perform spot gene analysis and caching.

    Parameters:
        folder_id: unique folder id for each sample

    Returns:
        None
    """
    thor, args, p_input_json = get_parameter(folder_id)

    sample_id = args['sampleId']
    selected_spot_gene_name = args['selectedSpotGeneName']
    sample_id_gene_spot = sample_id + "-" + selected_spot_gene_name

    thor.empty_cache_spot(sample_id, gene=True)

    args['sampleIdSpotGene'] = sample_id_gene_spot

    with open(f'../user{folder_id}/args.json', 'w') as f:
        json.dump(args, f)

    thor.spot_gis(
        sample_id=sample_id,
        selected_spot_gene_name=selected_spot_gene_name,
    )
    
    cache = cache_generate(sample_id, sample_id_gene_spot=sample_id_gene_spot)
    shutil.copy(os.path.join(cache_path, cache["gis-blend-spots"]), os.path.join(cache_path, cache["gis-blend-spots-gene"]))


def visualization_img_input(folder_id, data_path=data_path, cache_path=cache_path):
    """
    Visualize input image.

    Parameters:
        data_path (str): The path to the data directory.
        cache_path (str): The path to the cache directory.

    Returns:
        obj: The input map object.
    """
    thor, args, p_input_json = get_parameter(folder_id)
    sample_id_file = args['sampleIdFile']
    wsi_client, wsi_layer = thor.gis_client_and_layer(sample_id_file, 'gis-wsi-img')
    input_map = create_leaflet_map(
        'map-input',
        wsi_client,
        wsi_layer,
        [],
        cmax=1
    )
    return input_map


def visualization_img_all(folder_id, data_path=data_path, cache_path=cache_path):
    """
    Visualize all images including cell gene, cell type, and spot gene.

    Parameters:
        folder_id: unique folder id for each sample
        data_path (str): The path to the data directory.
        cache_path (str): The path to the cache directory.

    Returns:
        obj: The output map object.
    """
    thor, args, p_input_json = get_parameter(folder_id)
    sample_id_file = args['sampleIdFile']
    sample_id_gene_spot = args['sampleIdSpotGene']
    sample_id_gene_cell = args['sampleIdCellGene']

    cell_gene_client, cell_gene_layer = thor.gis_client_and_layer(sample_id_gene_cell, 'gis-blend-cell-gene-img')
    wsi_client, wsi_layer = thor.gis_client_and_layer(sample_id_file, 'gis-wsi-img')
    # cell_type_client, cell_type_layer = thor.gis_client_and_layer(sample_id_file, 'gis-blend-cell-type-img')
    _, spot_gene_layer = thor.gis_client_and_layer(sample_id_gene_spot, 'gis-blend-spot-gene-img')
    output_map = create_leaflet_map(
        'map-output',
        wsi_client,
        wsi_layer,
        [(spot_gene_layer, 'spot gene'), (cell_gene_layer, 'cell gene')],
        cmax=0 
        # (cell_type_layer, 'cell type')],
    )
    return output_map


all_gene_cell_layer = []
gene_cell = []


def visualization_img_cell(folder_id, data_path=data_path, cache_path=cache_path, cell_type=False, pathway=False, cell_detect=False):
    """
    Visualize cell images including cell gene and optionally cell type.

    Parameters:
        folder_id: unique folder id for each sample
        data_path (str): The path to the data directory.
        cache_path (str): The path to the cache directory.
        cell_type (bool): Whether to include cell type visualization.

    Returns:
        obj: The output map object.
    """
    thor, args, p_input_json = get_parameter(folder_id)
    selected_cell_gene_name = args['selectedCellGeneName']
    sample_id_file = args['sampleIdFile']
    sample_id_gene_cell = args['sampleIdCellGene']
    sample_id_pathway = args['sampleIdPathway']
    if len(all_gene_cell_layer) > 10:
        all_gene_cell_layer.clear()
        gene_cell.clear()
    if cell_type is True:
        wsi_client, wsi_layer = thor.gis_client_and_layer(sample_id_file, 'gis-wsi-img')
        _, cell_type_layer = thor.gis_client_and_layer(sample_id_file, 'gis-blend-cell-type-img')
        output_map = create_leaflet_map(
            'map-output',
            wsi_client,
            wsi_layer,
            [(cell_type_layer, 'cell type')],
            cmax=0
            
        )
    
    elif pathway is True:
        wsi_client, wsi_layer = thor.gis_client_and_layer(sample_id_file, 'gis-wsi-img')
        _, cell_pathway_heatmap_layer = thor.gis_client_and_layer(sample_id_pathway, 'gis-blend-cell-pathway-heatmap-img')
        output_map = create_leaflet_map(
            'map-output',
            wsi_client,
            wsi_layer,
            [(cell_pathway_heatmap_layer, f'{sample_id_pathway}')],
            cmax=0
        )

    elif cell_detect is True:
        _, cell_random_layer = thor.gis_client_and_layer(sample_id_file, 'gis-blend-cell-random-img')
        wsi_client, wsi_layer = thor.gis_client_and_layer(sample_id_file, 'gis-wsi-img')
        output_map = create_leaflet_map(
            'map-output',
            wsi_client,
            wsi_layer,
            [(cell_random_layer, 'cell detection')],
            cmax=0
        )
        
    elif cell_type is False and pathway is False and cell_detect is False:
        _, globals()[f"cell_gene_layer_{selected_cell_gene_name}"] = thor.gis_client_and_layer(sample_id_gene_cell, 'gis-blend-cell-gene-img')
        wsi_client, wsi_layer = thor.gis_client_and_layer(sample_id_file, 'gis-wsi-img')
        if selected_cell_gene_name in gene_cell:
            pass
        else:
            all_gene_cell_layer.append((globals()[f"cell_gene_layer_{selected_cell_gene_name}"], f'cell data: {selected_cell_gene_name}'))
            gene_cell.append(selected_cell_gene_name)
        output_map = create_leaflet_map(
            'map-output',
            wsi_client,
            wsi_layer,
            all_gene_cell_layer,
            cmax=0
        )
        
    return output_map


all_gene_spot_layer = []
gene_spot = []


def visualization_img_spot(folder_id, data_path=data_path, cache_path=cache_path):
    """
    Visualize spot images including spot gene.

    Parameters:
        folder_id: unique folder id for each sample
        data_path (str): The path to the data directory.
        cache_path (str): The path to the cache directory.

    Returns:
        obj: The output map object.
    """
    thor, args, p_input_json = get_parameter(folder_id)
    selected_spot_gene_name = args['selectedSpotGeneName']
    sample_id_file = args['sampleIdFile']
    sample_id_gene_spot = args['sampleIdSpotGene']
    if len(all_gene_spot_layer) > 10:
        all_gene_spot_layer.clear()
        gene_spot.clear()
    wsi_client, wsi_layer = thor.gis_client_and_layer(sample_id_file, 'gis-wsi-img')
    _, globals()[f"spot_gene_layer_{selected_spot_gene_name}"] = thor.gis_client_and_layer(sample_id_gene_spot, 'gis-blend-spot-gene-img')
    if selected_spot_gene_name in gene_spot:
        pass
    else:
        all_gene_spot_layer.append((globals()[f"spot_gene_layer_{selected_spot_gene_name}"], f'spot data: {selected_spot_gene_name}'))
        gene_spot.append(selected_spot_gene_name)
    output_map = create_leaflet_map(
        'map-output',
        wsi_client,
        wsi_layer,
        all_gene_spot_layer,
    )
    return output_map


# zip folder
def zip_folder(folder_path, output_path):
    """
    Compresses the contents of a specified folder into a ZIP file.

    Parameters:
        folder_path (str): The path to the folder containing the files to be compressed.
        output_path (str): The path to the output ZIP file.
        
    Returns:
        None
    """
    with zipfile.ZipFile(output_path, 'w', zipfile.ZIP_DEFLATED) as zipf:
        for root, dirs, files in os.walk(folder_path):
            for file in files:
                file_path = os.path.join(root, file)
                arcname = os.path.relpath(file_path, folder_path)
                zipf.write(file_path, arcname)
                

# clear cache    
def clear_cache(folder_id):
    """
    Clears cache files.

    Parameters:
        folder_id: unique folder id for each sample

    Returns:
        None or dash.no_update: If no cache files meet the conditions, returns dash.no_update.
            Otherwise, deletes cache files and returns None.
    """
    thor, args, p_input_json = get_parameter(folder_id)
    sample_id = args['sampleId']
    files = os.listdir(cache_path)
    for del_file in files:
        if not del_file.startswith("gt-iz-p9-rep2-file-name-gis-wsi-img.tiff") and del_file.startswith(sample_id):
            file_path = os.path.join(cache_path, del_file)
            os.remove(file_path)


# clear data
def clear_data(folder_id):
    """
    Clears previous data folder.

    Parameters:
        folder_id: unique folder id for each sample

    Returns:
        None or dash.no_update: If no cache files meet the conditions, returns dash.no_update.
            Otherwise, deletes cache files and returns None.
    """
    thor, args, p_input_json = get_parameter(folder_id)
    sample_id = args['sampleId']
    files = os.listdir(data_path)
    for del_file in files:
        if not del_file.startswith("gt-iz-p9-rep2") and del_file.startswith(sample_id):
            file_path = os.path.join(data_path, del_file)
            os.remove(file_path)

    try:
        if folder_id == "":
            pass
        else:
            shutil.rmtree(f"../user{folder_id}")
    except FileNotFoundError:
        pass
