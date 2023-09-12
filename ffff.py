"""File."""

import toml

config = toml.load('config.toml')
primary_key_list = ['gt-iz-p9-rep2', 'gt-iz-p7-rep2']
data_extension = {
    'cell': 'h5ad',
    'cell-gene': 'npz',
    'cell-gene-name': 'txt',
    'cell-info': 'csv',
    'cell-mask': 'npz',
    'wsi-img': 'tiff',
    'spot': 'h5ad',
    'spot-gene': 'npz',
    'spot-gene-name': 'txt',
    'spot-info': 'csv',
}
cache_extension = {
    'blend-cell-gene-img': 'png',
    'blend-cell-type-img': 'png',
    'blend-spot-gene-img': 'png',
    'mask-cell-gene-img': 'png',
    'mask-cell-type-img': 'png',
    'circle-spot-gene-img': 'png',
    'gis-blend-cell-gene-img': 'tiff',
    'gis-blend-cell-type-img': 'tiff',
    'gis-blend-spot-gene-img': 'tiff',
}
