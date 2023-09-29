import sys
sys.path.append("../")
from niceview.utils.dataset import ThorQuery
from niceview.pyplot.leaflet import create_leaflet_map
import shutil
import json
import toml
import os
import zipfile

config = toml.load('../user/config.toml')
data_path = config['path']['data']
cache_path = config['path']['cache']
max_file_size = config['constant']['max_file_size']


def get_default_para():
    """
    Sets the application parameters to their default values.

    This function loads default parameter values from respective JSON files and overwrites the current
    application parameters with these default values.

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
    # with open('../db/db-info-default.json', 'r') as info:
    #     db_info_default = json.load(info)
    # with open('../db/db-info.json', 'w') as info:
    #     json.dump(db_info_default, info)


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

        'gis-img-file': '-'.join([sample_id_file, 'gis-wsi-img.tiff']),
        'gis-blend-cells-gene': '-'.join([sample_id_gene_cell, 'gis-blend-cell-gene-img.tiff']),
        'gis-blend-spots-gene': '-'.join([sample_id_gene_spot, 'gis-blend-spot-gene-img.tiff']),
        'gis-blend-cell-type-file': '-'.join([sample_id_file, 'gis-blend-cell-type-img.tiff']),
        'gis-blend-cell-pathway-heatmap': '-'.join([sample_id_pathway, 'gis-blend-cell-pathway-heatmap-img.tiff']),
    }
    return cache


def get_parameter():
    """
    Get parameters from configuration files and create a ThorQuery object.

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
    with open('../user/args.json') as f:
        args = json.load(f)

    thor = ThorQuery(
        data_path,
        cache_path,
        data_extension,
        cache_extension,
        cell_label_encoder,
        cell_label_cmap,
        primary_key_list,
    )
    
    return thor, args


def get_wsi():
    """
    Get client for wsi image and perform caching.

    Returns:
        None
    """
    thor, args = get_parameter()
    sample_id = args['sampleId']
    sample_id_file = args['sampleId'] + '-' + args['fileName']

    args["sampleIdFile"] = sample_id_file

    with open('../user/args.json', 'w') as f:
        json.dump(args, f)
    thor.wsi_gis(sample_id)
    cache = cache_generate(sample_id, sample_id_file=sample_id_file)
    shutil.copy(os.path.join(cache_path, cache["gis-img"]), os.path.join(cache_path, cache["gis-img-file"]))
    

def calculation_cell(label_analysis=True):
    """
    Perform cell gene analysis and caching.

    Parameters:
        label_analysis (bool): Whether to perform cell type analysis.

    Returns:
        None
    """
    thor, args = get_parameter()
    sample_id = args['sampleId']
    sample_id_file = args['sampleIdFile']
    selected_cell_gene_name = args['selectedCellGeneName']
    sample_id_gene_cell = sample_id + "-" + selected_cell_gene_name

    thor.empty_cache_cell(sample_id, gene=True, label=True)
    
    args['sampleIdCellGene'] = sample_id_gene_cell

    with open('../user/args.json', 'w') as f:
        json.dump(args, f)

    thor.cell_gis(
        sample_id,
        selected_cell_gene_name,
        label_analysis=label_analysis,
    )
    
    cache = cache_generate(sample_id, sample_id_gene_cell=sample_id_gene_cell, sample_id_file=sample_id_file)
    shutil.copy(os.path.join(cache_path, cache["gis-blend-cells"]), os.path.join(cache_path, cache["gis-blend-cells-gene"]))
    # shutil.copy(os.path.join(cache_path, cache["gis-blend-cell-type"]), os.path.join(cache_path, cache["gis-blend-cell-type-file"]))


def calculation_pathway(label_analysis=False):
    """
    Perform cell gene analysis and caching.

    Parameters:
        label_analysis (bool): Whether to perform cell type analysis.

    Returns:
        None
    """
    thor, args = get_parameter()
    sample_id = args['sampleId']
    selected_pathway = args["selectedPathway"]
    
    sample_id_pathway = sample_id + "-" + selected_pathway

    thor.empty_cache_cell(sample_id, pathway=True)
    
    args['sampleIdPathway'] = sample_id_pathway

    with open('../user/args.json', 'w') as f:
        json.dump(args, f)

    thor.cell_gis(
        sample_id=sample_id,
        selected_pathway=selected_pathway,
        label_analysis=label_analysis,
    )
    
    cache = cache_generate(sample_id, sample_id_pathway=sample_id_pathway)
    shutil.copy(os.path.join(cache_path, cache["gis-blend-cell-heatmap"]), os.path.join(cache_path, cache["gis-blend-cell-pathway-heatmap"]))
    # shutil.copy(os.path.join(cache_path, cache["gis-blend-cell-type"]), os.path.join(cache_path, cache["gis-blend-cell-type-file"]))


def calculation_CNV(label_analysis=True):
    """
    Perform cell gene analysis and caching.

    Parameters:
        label_analysis (bool): Whether to perform cell type analysis.

    Returns:
        None
    """
    thor, args = get_parameter()
    sample_id = args['sampleId']
    sample_id_file = args['sampleIdFile']

    with open('../user/previous-input.json', 'r') as p_input:
        p_input_json = json.load(p_input)

    if sample_id_file not in p_input_json["CellType"]:
        p_input_json["CellType"].append(sample_id_file)

    with open('../user/previous-input.json', 'w') as p_input:
        json.dump(p_input_json, p_input)

    thor.empty_cache_cell(sample_id, label=True)

    thor.cell_gis(
        sample_id=sample_id,
        label_analysis=label_analysis,
    )
    
    cache = cache_generate(sample_id, sample_id_file=sample_id_file)
    shutil.copy(os.path.join(cache_path, cache["gis-blend-cell-type"]), os.path.join(cache_path, cache["gis-blend-cell-type-file"]))


def calculation_spot():
    """
    Perform spot gene analysis and caching.

    Returns:
        None
    """
    thor, args = get_parameter()

    sample_id = args['sampleId']
    selected_spot_gene_name = args['selectedSpotGeneName']
    sample_id_gene_spot = sample_id + "-" + selected_spot_gene_name

    thor.empty_cache_spot(sample_id, gene=True)

    args['sampleIdSpotGene'] = sample_id_gene_spot

    with open('../user/args.json', 'w') as f:
        json.dump(args, f)

    thor.spot_gis(
        sample_id=sample_id,
        selected_spot_gene_name=selected_spot_gene_name,
    )
    
    cache = cache_generate(sample_id, sample_id_gene_spot=sample_id_gene_spot)
    shutil.copy(os.path.join(cache_path, cache["gis-blend-spots"]), os.path.join(cache_path, cache["gis-blend-spots-gene"]))


def visualization_img_input(data_path, cache_path):
    """
    Visualize input image.

    Parameters:
        data_path (str): The path to the data directory.
        cache_path (str): The path to the cache directory.

    Returns:
        obj: The input map object.
    """
    thor, args = get_parameter()
    sample_id_file = args['sampleIdFile']
    wsi_client, wsi_layer = thor.gis_client_and_layer(sample_id_file, 'gis-wsi-img')
    input_map = create_leaflet_map(
        'map-input',
        wsi_client,
        wsi_layer,
        [],
        cmax=None
    )
    return input_map


def visualization_img_all(data_path, cache_path):
    """
    Visualize all images including cell gene, cell type, and spot gene.

    Parameters:
        data_path (str): The path to the data directory.
        cache_path (str): The path to the cache directory.

    Returns:
        obj: The output map object.
    """
    thor, args = get_parameter()
    sample_id_file = args['sampleIdFile']
    sample_id_gene_spot = args['sampleIdSpotGene']
    sample_id_gene_cell = args['sampleIdCellGene']

    cell_gene_client, cell_gene_layer = thor.gis_client_and_layer(sample_id_gene_cell, 'gis-blend-cell-gene-img')
    wsi_client, wsi_layer = thor.gis_client_and_layer(sample_id_file, 'gis-wsi-img')
    # cell_type_client, cell_type_layer = thor.gis_client_and_layer(sample_id_file, 'gis-blend-cell-type-img')
    spot_gene_client, spot_gene_layer = thor.gis_client_and_layer(sample_id_gene_spot, 'gis-blend-spot-gene-img')
    output_map = create_leaflet_map(
        'map-output',
        wsi_client,
        wsi_layer,
        [(spot_gene_layer, 'spot gene'), (cell_gene_layer, 'cell gene')] 
        # (cell_type_layer, 'cell type')],
    )
    return output_map


def visualization_img_cell(data_path, cache_path, cell_type=False, pathway=False):
    """
    Visualize cell images including cell gene and optionally cell type.

    Parameters:
        data_path (str): The path to the data directory.
        cache_path (str): The path to the cache directory.
        cell_type (bool): Whether to include cell type visualization.

    Returns:
        obj: The output map object.
    """
    thor, args = get_parameter()
    sample_id = args['sampleId']
    selected_cell_gene_name = args['selectedCellGeneName']
    sample_id_file = args['sampleIdFile']
    sample_id_gene_cell = args['sampleIdCellGene']
    sample_id_pathway = args['sampleIdPathway']

    if cell_type is True:
        wsi_client, wsi_layer = thor.gis_client_and_layer(sample_id_file, 'gis-wsi-img')
        _, cell_type_layer = thor.gis_client_and_layer(sample_id_file, 'gis-blend-cell-type-img')
        output_map = create_leaflet_map(
            'map-output',
            wsi_client,
            wsi_layer,
            [(cell_type_layer, 'cell type')]
            
        )
    
    elif pathway is True:
        wsi_client, wsi_layer = thor.gis_client_and_layer(sample_id_file, 'gis-wsi-img')
        _, cell_pathway_heatmap_layer = thor.gis_client_and_layer(sample_id_pathway, 'gis-blend-cell-pathway-heatmap-img')
        output_map = create_leaflet_map(
            'map-output',
            wsi_client,
            wsi_layer,
            [(cell_pathway_heatmap_layer, f'{sample_id_pathway}')],
            cmax=thor.get_gene_max(sample_id, selected_cell_gene_name)
        )
    elif cell_type is False and pathway is False:
        _, cell_gene_layer = thor.gis_client_and_layer(sample_id_gene_cell, 'gis-blend-cell-gene-img')
        wsi_client, wsi_layer = thor.gis_client_and_layer(sample_id_file, 'gis-wsi-img')
        output_map = create_leaflet_map(
            'map-output',
            wsi_client,
            wsi_layer,
            [(cell_gene_layer, f'{sample_id_gene_cell}')]
        )
    return output_map


def visualization_img_spot(data_path, cache_path):
    """
    Visualize spot images including spot gene.

    Parameters:
        data_path (str): The path to the data directory.
        cache_path (str): The path to the cache directory.

    Returns:
        obj: The output map object.
    """
    thor, args = get_parameter()
    sample_id = args['sampleId']
    selected_spot_gene_name = args['selectedSpotGeneName']
    sample_id_file = args['sampleIdFile']
    sample_id_gene_spot = args['sampleIdSpotGene']

    wsi_client, wsi_layer = thor.gis_client_and_layer(sample_id_file, 'gis-wsi-img')
    _, spot_gene_layer = thor.gis_client_and_layer(sample_id_gene_spot, 'gis-blend-spot-gene-img')
    output_map = create_leaflet_map(
        'map-output',
        wsi_client,
        wsi_layer,
        [(spot_gene_layer, f'{sample_id_gene_spot}')],
        cmax=thor.get_gene_max(sample_id, selected_spot_gene_name)
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
def clear_cache():
    """
    Clears cache files.

    Returns:
        None or dash.no_update: If no cache files meet the conditions, returns dash.no_update.
            Otherwise, deletes cache files and returns None.
    """
    files = os.listdir(cache_path)
    for del_file in files:
        if not del_file.startswith("gt-iz-p9-rep2-file-name-gis-wsi-img.tiff") and os.path.isfile(os.path.join(cache_path, del_file)):
            file_path = os.path.join(cache_path, del_file)
            print(file_path)
            os.remove(file_path)


# clear data
def clear_data():
    """
    Clears previous data folder.

    Returns:
        None or dash.no_update: If no cache files meet the conditions, returns dash.no_update.
            Otherwise, deletes cache files and returns None.
    """
    with open('../user/args-default.json') as f:
        args_default = json.load(f)
    with open('../user/args.json', 'w') as f:
        json.dump(args_default, f)
    files = os.listdir(data_path)
    for del_file in files:
        if not del_file.startswith("gt-iz-p9-rep2") and os.path.isfile(os.path.join(data_path, del_file)):
            file_path = os.path.join(data_path, del_file)
            os.remove(file_path)
