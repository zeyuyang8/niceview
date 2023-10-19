"""Leaflet."""

import dash_leaflet as dl
from niceview.utils.tools import CMAX, CMIN, get_hex_values


def create_leaflet_map(
    map_id,
    base_client,
    base_layer,
    list_of_layers,
    cmax=CMAX,
    classes=None,
):
    """Create leaflet map.
    
    Args:
        map_id (str): Map ID.
        base_client (TileClient): Base client.
        base_layer (TileLayer): Base layer.
        list_of_layers (list[tuple]): List of layers.
        cmax (int, optional): Max value.
        classes (int): List of classes.
        
    Returns:
        Map
    """
    # viewport
    default_center = [base_client.center()[0], base_client.center()[1]]
    max_zoom = base_client.max_zoom
    
    # viewport bounds
    expanded_bounds = list(base_layer.bounds)
    expanded_bounds[0] = list(expanded_bounds[0])
    expanded_bounds[1] = list(expanded_bounds[1])
    vert_dst = (expanded_bounds[1][1] - expanded_bounds[0][1]) 
    hori_dst = (expanded_bounds[1][0] - expanded_bounds[0][0])
    expanded_bounds[0][1] -= vert_dst
    expanded_bounds[1][0] += vert_dst
    expanded_bounds[0][0] -= hori_dst
    expanded_bounds[1][1] += hori_dst

    # zoom factor
    zoom_factor = 1.5 * ((vert_dst + hori_dst) / 2) / 0.0085
    default_zoom = base_client.default_zoom + zoom_factor
    
    # overlay
    overlay_layers = []
    for arg_layer, arg_name in list_of_layers:
        layer = dl.Overlay(
            dl.TileLayer(
                opacity=1,
                url=arg_layer.url,
                maxZoom=max_zoom,
                minZoom=default_zoom,
            ),
            name=arg_name,
        )
        overlay_layers.append(layer)
    
    # create map
    width, height = 20, 200
    thor_map = dl.Map(
        id=map_id,
        children=[
            dl.TileLayer(
                url=base_layer.url,
                maxZoom=max_zoom,
                minZoom=default_zoom,
            ),
            dl.Colorbar(
                colorscale=get_hex_values('jet'), width=width, height=height, min=CMIN, max=cmax, position='bottomleft', classes=classes,
            ),
            dl.LayersControl(
                overlay_layers, hideSingleBase=True,
            ),
            dl.EasyButton(icon='fa-home', id='btn_home'),
            dl.FullScreenControl(),
            dl.FeatureGroup(
                [
                    dl.EditControl(
                        id='editControl',
                        draw={
                            'polyline': False,
                            'polygon': True,
                            'rectangle': True,
                            'circle': False,
                            'circlemarker': False,
                            'marker': False,
                        },
                    ),
                ],
            ),
            dl.EasyButton(icon='fa-save', id='btn_save'),
            dl.ScaleControl(imperial=False),
        ],
        center=default_center,
        zoom=default_zoom,
        style={'height': '850px', 'margin': 'auto', 'display': 'block', 'background': 'black'},
        attributionControl=False,
        trackViewport=True,
        maxBounds=expanded_bounds,
    )
    return thor_map
