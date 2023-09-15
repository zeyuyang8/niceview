"""Leaflet."""

import dash_leaflet as dl


def create_leaflet_map(
    map_id,
    base_client,
    base_layer,
    list_of_layers,
):
    """Create leaflet map.
    
    Args:
        map_id (str): Map ID.
        base_client (TileClient): Base client.
        base_layer (TileLayer): Base layer.
        list_of_layers (list[tuple]): List of layers.
        
    Returns:
        Map
    """
    # viewport
    default_center = [base_client.center()[0], base_client.center()[1]]
    max_zoom = base_client.max_zoom
    default_zoom = base_client.default_zoom + 0.5
    
    # overlay
    overlay_layers = []
    for arg_layer, arg_name in list_of_layers:
        layer = dl.Overlay(
            dl.TileLayer(
                url=arg_layer.url,
                maxZoom=max_zoom,
                minZoom=default_zoom,
            ),
            name=arg_name,
        )
        overlay_layers.append(layer)
    
    # create map
    thor_map = dl.Map(
        id=map_id,
        children=[
            dl.TileLayer(
                url=base_layer.url,
                maxZoom=max_zoom,
                minZoom=default_zoom,
            ),
            dl.LayersControl(
                overlay_layers, hideSingleBase=True, position='bottomright',
            ),
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
        ],
        center=default_center,
        zoom=default_zoom,
        style={'height': '70vh', 'margin': 'auto', 'display': 'block'},
        attributionControl=False,
        trackViewport=True,
    )
    return thor_map
