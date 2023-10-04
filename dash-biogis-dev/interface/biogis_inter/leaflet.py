"""Leaflet."""


import matplotlib.pyplot as plt

CMIN = 0
CMAX = 255
external_css = "https://cdnjs.cloudflare.com/ajax/libs/font-awesome/4.7.0/css/font-awesome.min.css"


def get_hex_values(colormap_name):
    """Get hex values.
    
    Args:
        colormap_name (str): Colormap name.
    
    Returns:
        list[str]: List of hex values.
    """
    cmap = plt.get_cmap(colormap_name)
    hex_values = []
    for i in range(cmap.N):
        rgba = cmap(i)
        hex_color = '#{:02X}{:02X}{:02X}'.format(int(rgba[0] * CMAX), int(rgba[1] * CMAX), int(rgba[2] * CMAX))
        hex_values.append(hex_color)
    return hex_values


def create_leaflet_map(
    map_id,
    base_client,
    base_layer,
    list_of_layers,
    cmax=CMAX,
):
    """Create leaflet map.
    
    Args:
        map_id (str): Map ID.
        base_client (TileClient): Base client.
        base_layer (TileLayer): Base layer.
        list_of_layers (list[tuple]): List of layers.
        cmax (int, optional): Max value.
        
    Returns:
        Map
    """
    #import interface.dash_leaflet.dash_leaflet as dl
    import dash_leaflet as dl
    # viewport
    default_center = [base_client.center()[0], base_client.center()[1]]
    max_zoom = base_client.max_zoom
    default_zoom = base_client.default_zoom + 0.5
    
    # overlay
    overlay_layers = []
    for arg_layer, arg_name in list_of_layers:
        layer = dl.Overla(
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
                colorscale=get_hex_values('jet'), width=width, height=height, min=CMIN, max=cmax, position='bottomleft',
            ),
            dl.LayersControl(
                overlay_layers, hideSingleBase=True,
            ),
            dl.EasyButton(icon="fa-home", id="btn_home"),
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
            dl.EasyButton(icon="fa-save", id="btn_save"),
            dl.ScaleControl(imperial=False),
        ],
        center=default_center,
        zoom=default_zoom,
        style={'height': '850px', 'margin': 'auto', 'display': 'block', 'background': 'black'},
        attributionControl=False,
        trackViewport=True,
    )
    return thor_map
