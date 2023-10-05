"""Dataset utilities."""

import os
import cv2
import rasterio
import scipy
import pandas as pd
import numpy as np
import scanpy as sc
import PIL
from PIL import Image
from scipy.sparse import load_npz
from localtileserver import TileClient, get_leaflet_tile_layer
from niceview.utils.tools import txt_to_list, select_col_from_name, normalize_array
from niceview.utils.tools import mask_filter_relabel, mask_to_image, discrete_cmap_from_hex
from niceview.utils.tools import blend, draw_circles
from niceview.utils.raster import geo_ref_raster
from niceview.utils.cell import get_nuclei_pixels
from niceview.pyplot.heatmap import heatmap_from_scatter
from scipy import ndimage
from pyproj import Geod

Image.MAX_IMAGE_PIXELS = None
CMAX = 255
CMIN = 1  # avoid zero to distinguish from background


class AristotleDataset:
    """Aristotle dataset."""
    
    def __init__(self, data_dir, data_extension, cache_dir, cache_extension, primary_key_list):
        """Initialize Aristotle dataset.
        
        Args:
            data_dir (str): data directory.
            data_extension (dict): data extension.
            cache_dir (str): cache directory.
            cache_extension (dict): cache extension.
            primary_key_list (list of str): list of primary keys.
        """
        self.data_dir = data_dir
        self.data_extension = data_extension
        self.cache_dir = cache_dir
        self.cache_extension = cache_extension
        self.primary_key_list = primary_key_list
    
    def get_data_field(self, primary_key, data_field):
        """Get data field.
        
        Raise:
            ValueError: bad input primary key.
        
        Args:
            primary_key (str): primary key.
            data_field (str): data field.
        
        Returns:
            str: data field path.
        """
        # if primary_key not in self.primary_key_list:
        #     raise ValueError('Bad input primary key')
        
        filename = self._unparse_filename(primary_key, data_field, self.data_extension[data_field])
        filepath = os.path.join(self.data_dir, filename)
        return filepath
    
    def get_cache_field(self, primary_key, cache_field):
        """Get cache field.

        Raise:
            ValueError: bad input primary key.
        
        Args:
            primary_key (str): primary key.
            cache_field (str): cache field.
            
        Returns:
            str: cache field path.
        """
        # if primary_key not in self.primary_key_list:
        #     raise ValueError('Bad input primary key')
        
        filename = self._unparse_filename(
            primary_key, cache_field, self.cache_extension[cache_field],
        )
        filepath = os.path.join(self.cache_dir, filename)
        return filepath
    
    def _unparse_filename(self, primary_key, field_name, extension):
        """Unparse filename.
        
        Args:
            primary_key (str): primary key.
            field_name (str): field name.
            extension (str): extension.
        
        Returns:
            str: filename.
        """
        filename = '-'.join([primary_key, field_name])
        filename = '.'.join([filename, extension])
        return filename


class ThorQuery:
    """Container for query."""
    
    def __init__(
        self,
        data_path,
        cache_path, 
        data_extension, 
        cache_extension, 
        cell_label_encoder, 
        cell_label_cmap, 
        primary_key_list,
    ):
        """Initialize query.
        
        Args:
            data_path (str): data path.
            cache_path (str): cache path.
            data_extension (dict): data extension.
            cache_extension (dict): cache extension.
            cell_label_encoder (dict): cell label encoder.
            cell_label_cmap (dict): cell label colormap.
            primary_key_list (list of str): list of primary keys.
        """
        self._data_path = data_path
        self._cache_path = cache_path
        self._data_extension = data_extension
        self._cache_extension = cache_extension
        self._cell_label_encoder = cell_label_encoder
        self._cell_label_cmap = cell_label_cmap
        self._primary_key_list = primary_key_list
    
        self.dataset = AristotleDataset(
            data_path,
            data_extension,
            cache_path,
            cache_extension,
            primary_key_list,
        )
    
    def cell_analysis(self, sample_id, selected_cell_gene_name=None, label_analysis=False, heatmap=False, selected_pathway=None):
        """Cell gene analysis.
        
        Args:
            sample_id (str): sample id.
            selected_cell_gene_name (str): list of selected cell gene name.
            label_analysis (bool): whether to label the cell.
            heatmap (bool): whether to generate heatmap.
            selected_pathway (str): pathway name.
        """
        cell_info = pd.read_csv(
            self.dataset.get_data_field(sample_id, 'cell-info'),
        )
        cell_pos = cell_info[['x', 'y']].values
        if os.path.exists(self.dataset.get_cache_field(sample_id, 'mask-cell-match-region')):
            cell_matched_region = np.load(
                self.dataset.get_cache_field(sample_id, 'mask-cell-match-region'),
                allow_pickle=True,
            )
        else:
            cell_matched_region = get_nuclei_pixels(
                load_npz(
                    self.dataset.get_data_field(sample_id, 'cell-mask'),
                ).tocsr()[:, :].todense(),
                cell_pos,
            )
            np.save(
                self.dataset.get_cache_field(sample_id, 'mask-cell-match-region'),
                cell_matched_region,
                allow_pickle=True,
            )
        
        # random color
        if not os.path.exists(self.dataset.get_cache_field(sample_id, 'mask-cell-random-img')):
            cv2.imwrite(
                self.dataset.get_cache_field(sample_id, 'mask-cell-random-img'),
                mask_to_image(
                    mask_filter_relabel(
                        self.dataset.get_data_field(sample_id, 'cell-mask'),
                        cell_matched_region,
                        np.random.randint(CMIN, CMAX, len(cell_matched_region)),
                    ),
                    cv2.COLORMAP_JET,
                ),
            )

        # gene
        if selected_cell_gene_name:
            cell_gene = load_npz(
                self.dataset.get_data_field(sample_id, 'cell-gene'),
            )
            cell_gene_name = txt_to_list(
                self.dataset.get_data_field(sample_id, 'cell-gene-name'),
            )
            cell_selected_gene = select_col_from_name(
                cell_gene, cell_gene_name, selected_cell_gene_name,
            )
            cell_selected_gene_norm = normalize_array(cell_selected_gene, CMIN, CMAX)
            if not os.path.exists(self.dataset.get_cache_field(sample_id, 'mask-cell-gene-img')):
                cv2.imwrite(
                    self.dataset.get_cache_field(sample_id, 'mask-cell-gene-img'),
                    mask_to_image(
                        mask_filter_relabel(
                            self.dataset.get_data_field(sample_id, 'cell-mask'),
                            cell_matched_region,
                            cell_selected_gene_norm,
                        ),
                        cv2.COLORMAP_JET,
                    ),
                )
        
        # label
        if label_analysis:
            cell_label = [
                self._cell_label_encoder[x] for x in cell_info['label'].values
            ]
            if not os.path.exists(self.dataset.get_cache_field(sample_id, 'mask-cell-type-img')):
                cv2.imwrite(
                    self.dataset.get_cache_field(sample_id, 'mask-cell-type-img'),
                    mask_to_image(
                        mask_filter_relabel(
                            self.dataset.get_data_field(sample_id, 'cell-mask'),
                            cell_matched_region,
                            cell_label,
                        ),
                        discrete_cmap_from_hex(self._cell_label_cmap),
                    ),
                )
        
        # heatmap
        if heatmap:
            cell_gene = load_npz(
                self.dataset.get_data_field(sample_id, 'cell-gene'),
            )
            cell_gene_name = txt_to_list(
                self.dataset.get_data_field(sample_id, 'cell-gene-name'),
            )
            cell_selected_gene = select_col_from_name(
                cell_gene, cell_gene_name, selected_cell_gene_name,
            )
            cell_selected_gene_norm = normalize_array(cell_selected_gene, CMIN, CMAX)
            cell_gene_color = np.array(cell_selected_gene_norm).ravel()
            if not os.path.exists(self.dataset.get_cache_field(sample_id, 'cell-gene-heatmap-img')):
                image = PIL.Image.open(self.dataset.get_data_field(sample_id, 'wsi-img'))
                xmax, ymax = image.size
                heatmap_from_scatter(
                    (xmax, ymax), np.round(cell_pos).astype(int), cell_gene_color, 
                    dst_path=self.dataset.get_cache_field(sample_id, 'cell-gene-heatmap-img'),
                )

        if selected_pathway:
            cell_pathway_matrix = np.load(self.dataset.get_data_field(sample_id, 'cell-pathway-matrix'))
            cell_pathway_name = txt_to_list(self.dataset.get_data_field(sample_id, 'cell-pathway-name'))
            cell_selected_pathway = select_col_from_name(cell_pathway_matrix, cell_pathway_name, selected_pathway)
            cell_selected_pathway_norm = normalize_array(cell_selected_pathway, CMIN, CMAX)
            cell_pathway_color = np.array(cell_selected_pathway_norm).ravel()
            if not os.path.exists(self.dataset.get_cache_field(sample_id, 'cell-pathway-heatmap-img')):
                image = PIL.Image.open(self.dataset.get_data_field(sample_id, 'wsi-img'))
                xmax, ymax = image.size
                heatmap_from_scatter(
                    (xmax, ymax), np.round(cell_pos).astype(int), cell_pathway_color, 
                    dst_path=self.dataset.get_cache_field(sample_id, 'cell-pathway-heatmap-img'),
                )
    
    def cell_blend(
        self, sample_id, selected_cell_gene_name=None, label_analysis=False, 
        heatmap_analysis=False, selected_pathway=None, mask_opacity=1,
    ):
        """Cell blend.

        Args:
            sample_id (str): sample id.
            selected_cell_gene_name (str): list of selected cell gene name.
            label_analysis (bool): whether to label the cell.
            heatmap_analysis (bool): whether to generate heatmap.
            selected_pathway (str): pathway name.
            mask_opacity (float): mask opacity.
        """
        # analysis
        self.cell_analysis(sample_id, selected_cell_gene_name, label_analysis, heatmap_analysis, selected_pathway)

        # random color
        if not os.path.exists(self.dataset.get_cache_field(sample_id, 'blend-cell-random-img')):
            cv2.imwrite(
                self.dataset.get_cache_field(sample_id, 'blend-cell-random-img'),
                blend(
                    self.dataset.get_data_field(sample_id, 'wsi-img'),
                    self.dataset.get_cache_field(sample_id, 'mask-cell-random-img'),
                    mask_opacity,
                ),
            )

        if selected_cell_gene_name:
            if not os.path.exists(self.dataset.get_cache_field(sample_id, 'blend-cell-gene-img')):
                cv2.imwrite(
                    self.dataset.get_cache_field(sample_id, 'blend-cell-gene-img'),
                    blend(
                        self.dataset.get_data_field(sample_id, 'wsi-img'),
                        self.dataset.get_cache_field(sample_id, 'mask-cell-gene-img'),
                        mask_opacity,
                    ),
                )
        
        if label_analysis:
            if not os.path.exists(self.dataset.get_cache_field(sample_id, 'blend-cell-type-img')):
                cv2.imwrite(
                    self.dataset.get_cache_field(sample_id, 'blend-cell-type-img'),
                    blend(
                        self.dataset.get_data_field(sample_id, 'wsi-img'),
                        self.dataset.get_cache_field(sample_id, 'mask-cell-type-img'),
                        mask_opacity,
                    ),
                )

        if heatmap_analysis:
            if not os.path.exists(self.dataset.get_cache_field(sample_id, 'blend-cell-gene-heatmap-img')):
                cv2.imwrite(
                    self.dataset.get_cache_field(sample_id, 'blend-cell-gene-heatmap-img'),
                    blend(
                        self.dataset.get_data_field(sample_id, 'wsi-img'),
                        self.dataset.get_cache_field(sample_id, 'cell-gene-heatmap-img'),
                        0.5,
                    ),
                )
        
        if selected_pathway:
            if not os.path.exists(self.dataset.get_cache_field(sample_id, 'blend-cell-pathway-heatmap-img')):
                cv2.imwrite(
                    self.dataset.get_cache_field(sample_id, 'blend-cell-pathway-heatmap-img'),
                    blend(
                        self.dataset.get_data_field(sample_id, 'wsi-img'),
                        self.dataset.get_cache_field(sample_id, 'cell-pathway-heatmap-img'),
                        0.5,
                    ),
                )

    def cell_gis(self, sample_id, selected_cell_gene_name=None, label_analysis=False, heatmap_analysis=False, selected_pathway=None, mask_opacity=1):
        """Cell blend.

        Args:
            sample_id (str): sample id.
            selected_cell_gene_name (str): list of selected cell gene name.
            label_analysis (bool): whether to label the cell.
            heatmap_analysis (bool): whether to generate heatmap.
            selected_pathway (str): pathway name.
            mask_opacity (float): mask opacity.
        """
        # blend
        self.cell_blend(sample_id, selected_cell_gene_name, label_analysis, heatmap_analysis, selected_pathway, mask_opacity)
        
        # random color
        if not os.path.exists(self.dataset.get_cache_field(sample_id, 'gis-blend-cell-random-img')):
            geo_ref_raster(
                self.dataset.get_cache_field(sample_id, 'blend-cell-random-img'),
                self.dataset.get_cache_field(sample_id, 'gis-blend-cell-random-img'),
            )
            
        # georeference images for blended cell selected gene and cell type
        if selected_cell_gene_name:
            if not os.path.exists(self.dataset.get_cache_field(sample_id, 'gis-blend-cell-gene-img')):
                geo_ref_raster(
                    self.dataset.get_cache_field(sample_id, 'blend-cell-gene-img'),
                    self.dataset.get_cache_field(sample_id, 'gis-blend-cell-gene-img'),
                )
        
        if label_analysis:
            if not os.path.exists(self.dataset.get_cache_field(sample_id, 'gis-blend-cell-type-img')):
                geo_ref_raster(
                    self.dataset.get_cache_field(sample_id, 'blend-cell-type-img'),
                    self.dataset.get_cache_field(sample_id, 'gis-blend-cell-type-img'),
                )
        if heatmap_analysis:
            if not os.path.exists(self.dataset.get_cache_field(sample_id, 'gis-blend-cell-gene-heatmap-img')):
                geo_ref_raster(
                    self.dataset.get_cache_field(sample_id, 'blend-cell-gene-heatmap-img'),
                    self.dataset.get_cache_field(sample_id, 'gis-blend-cell-gene-heatmap-img'),
                )
        if selected_pathway:
            if not os.path.exists(self.dataset.get_cache_field(sample_id, 'gis-blend-cell-pathway-heatmap-img')):
                geo_ref_raster(
                    self.dataset.get_cache_field(sample_id, 'blend-cell-pathway-heatmap-img'),
                    self.dataset.get_cache_field(sample_id, 'gis-blend-cell-pathway-heatmap-img'),
                )

    def spot_analysis(self, sample_id, selected_spot_gene_name, thickness=-1):
        """Spot analysis.
        
        Args:
            sample_id (str): sample id
            selected_spot_gene_name (str): list of selected spot gene name.
            thickness (int): thickness of the circle.
        """
        # read image shape
        img_shape = load_npz(self.dataset.get_data_field(sample_id, 'cell-mask')).shape
        
        if selected_spot_gene_name:
            # spot info
            spot_info = pd.read_csv(self.dataset.get_data_field(sample_id, 'spot-info'))
            spot_pos = spot_info[['x', 'y']].values
            spot_diameter = spot_info['diameter'].values
            spot_gene = load_npz(self.dataset.get_data_field(sample_id, 'spot-gene'))
            spot_gene_name = txt_to_list(self.dataset.get_data_field(sample_id, 'spot-gene-name'))
            spot_selected_gene = select_col_from_name(
                spot_gene, spot_gene_name, selected_spot_gene_name,
            )
            spot_selected_gene_norm = normalize_array(spot_selected_gene, CMIN, CMAX)
            
            # draw circles
            if not os.path.exists(self.dataset.get_cache_field(sample_id, 'circle-spot-gene-img')):
                cv2.imwrite(
                    self.dataset.get_cache_field(sample_id, 'circle-spot-gene-img'),
                    draw_circles(
                        img_shape,
                        spot_pos,
                        spot_diameter,
                        spot_selected_gene_norm,
                        cmap=cv2.COLORMAP_JET,
                        thickness=-1,
                    ),
                )
    
    def spot_blend(self, sample_id, selected_spot_gene_name, thickness=-1, mask_opacity=1):
        """Spot blend.
        
        Args:
            sample_id (str): sample id.
            selected_spot_gene_name (str): list of selected spot gene name.
            thickness (int): thickness of the circle.
            mask_opacity (float): mask opacity.
        """
        # analysis
        self.spot_analysis(sample_id, selected_spot_gene_name, thickness)
        
        if selected_spot_gene_name:
            if not os.path.exists(self.dataset.get_cache_field(sample_id, 'blend-spot-gene-img')):
                cv2.imwrite(
                    self.dataset.get_cache_field(sample_id, 'blend-spot-gene-img'),
                    blend(
                        self.dataset.get_data_field(sample_id, 'wsi-img'),
                        self.dataset.get_cache_field(sample_id, 'circle-spot-gene-img'),
                        mask_opacity,
                    ),
                )
    
    def spot_gis(self, sample_id, selected_spot_gene_name, thickness=-1, mask_opacity=1):
        """Spot GIS.
        
        Args:
            sample_id (str): sample id.
            selected_spot_gene_name (str): list of selected spot gene name.
            thickness (int): thickness of the circle.
            mask_opacity (float): mask opacity.
        """
        # blend
        self.spot_blend(sample_id, selected_spot_gene_name, thickness, mask_opacity)
        
        # georeference images for blended spot selected gene
        if selected_spot_gene_name:
            if not os.path.exists(self.dataset.get_cache_field(sample_id, 'gis-blend-spot-gene-img')):
                geo_ref_raster(
                    self.dataset.get_cache_field(sample_id, 'blend-spot-gene-img'),
                    self.dataset.get_cache_field(sample_id, 'gis-blend-spot-gene-img'),
                )
    
    def wsi_gis(self, sample_id):
        """WSI GIS.
        
        Args:
            sample_id (str): sample id.
        """
        if not os.path.exists(self.dataset.get_cache_field(sample_id, 'gis-wsi-img')):
            geo_ref_raster(
                self.dataset.get_data_field(sample_id, 'wsi-img'),
                self.dataset.get_cache_field(sample_id, 'gis-wsi-img'),
            )

    def empty_cache(self, sample_id, cache_field):
        """Empty cache.
        
        Args:
            sample_id (str): sample id.
            cache_field (str): cache field.
        """
        os.remove(self.dataset.get_cache_field(sample_id, cache_field))

    def empty_cache_cell(self, sample_id, gene=False, label=False, heatmap=False, pathway=False):
        """Empty cell gene.
        
        Args:
            sample_id (str): sample id.
            gene (bool): whether to empty cell gene.
            label (bool): whether to empty cell label.
            heatmap (bool): whether to empty cell heatmap.
            pathway (bool): whether to empty cell pathway.
        """
        try:
            if gene:
                self.empty_cache(sample_id, 'mask-cell-gene-img')
                self.empty_cache(sample_id, 'blend-cell-gene-img')
                self.empty_cache(sample_id, 'gis-blend-cell-gene-img')

            if label:
                self.empty_cache(sample_id, 'mask-cell-type-img')
                self.empty_cache(sample_id, 'blend-cell-type-img')
                self.empty_cache(sample_id, 'gis-blend-cell-type-img')

            if heatmap:
                self.empty_cache(sample_id, 'cell-gene-heatmap-img')
                self.empty_cache(sample_id, 'blend-cell-gene-heatmap-img')
                self.empty_cache(sample_id, 'gis-blend-cell-gene-heatmap-img')

            if pathway:
                self.empty_cache(sample_id, 'cell-pathway-heatmap-img')
                self.empty_cache(sample_id, 'blend-cell-pathway-heatmap-img')
                self.empty_cache(sample_id, 'gis-blend-cell-pathway-heatmap-img')
                
        except FileNotFoundError:
            pass
    
    def empty_cache_spot(self, sample_id, gene=False):
        """Empty spot gene.
        
        Args:
            sample_id (str): sample id.
            gene (bool): whether to empty spot gene.
        """
        try:
            if gene:
                self.empty_cache(sample_id, 'circle-spot-gene-img')
                self.empty_cache(sample_id, 'blend-spot-gene-img')
                self.empty_cache(sample_id, 'gis-blend-spot-gene-img')

        except FileNotFoundError:
            pass

    def gis_client_and_layer(self, sample_id, cache_field):
        """GIS client.
        
        Args:
            sample_id (str): sample id.
            cache_field (str): cache field.
        
        Returns:
            client (TileClient): tile client.
            layer (TileLayer): tile layer.
        """
        client = TileClient(
            self.dataset.get_cache_field(sample_id, cache_field),
            cors_all=True,
        )
        layer = get_leaflet_tile_layer(client)
        return client, layer

    def get_coord_mapping(self, sample_id):
        """Get coordinate mapping from geographic to pixel.
        
        Args:
            sample_id (str): sample id.
        
        Returns:
            func: coordinate mapping function.
        """
        raster = rasterio.open(
            self.dataset.get_cache_field(sample_id, 'gis-wsi-img'),
        )
        return raster.index
    
    def get_sample_img_shape(self, sample_id):
        """Get sample image shape.
        
        Args:
            sample_id (str): sample id.
        
        Returns:
            tuple: height and width.
        """
        raster = rasterio.open(
            self.dataset.get_cache_field(sample_id, 'gis-wsi-img'),
        )
        return raster.shape

    def get_gene_max(self, sample_id, selected_cell_gene_name):
        """Get cmax.
        
        Args:
            sample_id (str): sample id.
            selected_cell_gene_name (str): list of selected cell gene name.
        
        Returns:
            float: max value.
        """
        cell_gene = load_npz(
            self.dataset.get_data_field(sample_id, 'cell-gene'),
        )
        cell_gene_name = txt_to_list(
            self.dataset.get_data_field(sample_id, 'cell-gene-name'),
        )
        cell_selected_gene = select_col_from_name(
            cell_gene, cell_gene_name, selected_cell_gene_name,
        )
        return float(cell_selected_gene.max())

    def get_cell_adata_and_img(self, sample_id):
        """Get cell adata.
        
        Args:
            sample_id (str): sample id.
        
        Returns:
            anndata.AnnData: cell adata.
            numpy.ndarray: image.
        """
        try:
            cell_adata = sc.read_h5ad(
                self.dataset.get_data_field(sample_id, 'cell'),
            )
        except FileNotFoundError:
            cell_adata = None
            pass
        img = cv2.imread(
            self.dataset.get_data_field(sample_id, 'wsi-img'),
        )
        return cell_adata, img
    
    def get_factor(self, gis_img_path, actual_distance=1e-6):
        """Get factor for converting pixel distance to actual distance.
        
        Args:
            gis_img_path (str): path to the GIS image.
            actual_distance (float): actual distance in micrometer.
            
        Returns:
            float: factor.
        """
        with rasterio.open(gis_img_path) as src:        
            lat1, lon1 = src.xy(0, 0)
            lat2, lon2 = src.xy(0, 1)
        g = Geod(ellps='clrk66') 
        _, _, dist = g.inv(lon1, lat1, lon2, lat2)
        factor = actual_distance / dist * 10 ** 6  # first convert to meter then convert to micrometer
        return factor

    def process_data(self, sample_id, height_width=None, img_path=None, mask_path=None, adata_path=None):
        """Process data.
        
        Args:
            sample_id (str): sample id.
            height_width (tuple): height and width.
            img_path (str): image path.
            mask_path (str): mask path.
            adata_path (str): adata path.
        
        Returns:
            tuple: height and width.
        """
        max_val = 10000
        if height_width is not None:
            height = height_width[0]
            width = height_width[1]
            max_dim = max(height, width)
            resize_factor = max_val / max_dim
        if img_path is not None:
            img = cv2.imread(img_path)
            height, width, _ = img.shape
            max_dim = max(height, width)
            if max_dim < max_val:
                return height, width
            sample_id = 'temp' + '-' + sample_id 
            # calculate resize factor
            resize_factor = max_val / max_dim
            resized_img = cv2.resize(img, (int(width * resize_factor), int(height * resize_factor)))
            # write data
            cv2.imwrite(self.dataset.get_data_field(sample_id, 'wsi-img'), resized_img)
        if max_dim < max_val:
            return height, width
        else:
            if mask_path is not None:
                # read mask
                mask = np.array(load_npz(mask_path).todense())
                original_shape = mask.shape
                target_shape = (int(height * resize_factor), int(width * resize_factor))  # (height, width) in numpy format

                # resize mask array
                resized_mask = ndimage.zoom(
                    mask, (target_shape[0] / original_shape[0], target_shape[1] / original_shape[1]), order=0,
                )
                # write data
                scipy.sparse.save_npz(self.dataset.get_data_field(sample_id, 'cell-mask'), scipy.sparse.csr_matrix(resized_mask))
            
            if adata_path is not None:
                # resize adata
                adata = sc.read_h5ad(adata_path)
                adata.obsm['spatial'] = adata.obsm['spatial'] * resize_factor
                
                # write data
                adata.write(self.dataset.get_data_field(sample_id, 'cell-temp'))
            return height, width
