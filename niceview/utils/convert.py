"""Convert."""

import json
import os
import shutil
import scanpy as sc
import scipy
import pandas as pd
from niceview.utils.tools import list_to_txt


def h5ad_converter(
    data_path, db_info_path, sample_id,
    h5ad_cell, h5ad_spot, cell_mask,
    delete_original=False,
):
    """Convert h5ad file to database format.
    
    Args:
        data_path (str): data path.
        db_info_path (str): database information path.
        sample_id (str): sample id.
        h5ad_cell (str): cell-wise h5ad file path.
        h5ad_spot (str): spot-wise h5ad file path.
        cell_mask (str): cell mask file path.
        delete_original (bool, optional): whether to delete original files. Defaults to False.

    Raises:
        ValueError: sample id already exists in database.
    """
    with open(db_info_path, 'r') as json_file:
        db_info = json.load(json_file)
    
    # updae primary key list in database information
    primary_key_list = db_info['primary_key_list']
    if sample_id not in primary_key_list:
        primary_key_list.append(sample_id)
        db_info['primary_key_list'] = primary_key_list
        with open(db_info_path, 'w') as json_file:
            json.dump(db_info, json_file)
    else:
        raise ValueError('sample id already exists in database.')
    
    data_extension = db_info['data_extension']
    data_file_names = {}
    for key, ext in data_extension.items():
        data_file_names[key] = f'{data_path}{sample_id}-{key}.{ext}'
    
    # rename h5ad file for cell-wise data
    shutil.copy2(h5ad_cell, data_file_names['cell'])
    
    # cell-wise data
    cell = sc.read_h5ad(data_file_names['cell'])
    scipy.sparse.save_npz(data_file_names['cell-gene'], cell.X)
    cell_gene_name = cell.var_names.to_list()
    list_to_txt(cell_gene_name, data_file_names['cell-gene-name'])
    cell_centroid = cell.obsm['spatial']
    cell_type = cell.obs['cell_type'].to_list()
    cell_info = pd.DataFrame(
        {
            'x': cell_centroid[:, 0],
            'y': cell_centroid[:, 1],
            'label': cell_type,
        },
    )
    cell_info.to_csv(data_file_names['cell-info'], index=False)
    shutil.copy2(cell_mask, data_file_names['cell-mask'])
    
    # rename h5ad file for spot-wise data
    shutil.copy2(h5ad_spot, data_file_names['spot'])
    
    # spot-wise data
    spot = sc.read_h5ad(data_file_names['spot'])
    scipy.sparse.save_npz(data_file_names['spot-gene'], spot.X)
    spot_gene_name = spot.var_names.to_list()
    list_to_txt(spot_gene_name, data_file_names['spot-gene-name'])
    spot_centroids = spot.obsm['spatial']
    spot_diameter = spot.uns['spatial']['Visium_19_CK297']['scalefactors']['spot_diameter_fullres']
    spot_info = pd.DataFrame(
        {
            'x': spot_centroids[:, 0],
            'y': spot_centroids[:, 1],
            'diameter': spot_diameter,
        },
    )
    spot_info.to_csv(data_file_names['spot-info'], index=False)
    
    if delete_original:
        os.remove(h5ad_cell)
        os.remove(h5ad_spot)
        os.remove(cell_mask)
