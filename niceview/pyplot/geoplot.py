"""Tile server."""

from localtileserver import TileClient, get_leaflet_tile_layer
import plotly.graph_objects as go


def scatter_overlay_georef_img(georef_img_path, coord_x, coord_y, plot_config):
    """Make scatter plot with overlay of georeferenced image.
    
    Args:
        georef_img_path (str): path to image.
        coord_x (np.ndarray): x coordinates.
        coord_y (np.ndarray): y coordinates.
        plot_config (dict): plot configuration.
    
    Returns:
        fig (plotly.graph_objects.Figure): Figure object.
    """
    marker_size = plot_config['marker_size']
    window_width = plot_config['window_width']
    window_height = plot_config['window_height']
    
    map_client = TileClient(georef_img_path, cors_all=True)
    map_tile_layer = get_leaflet_tile_layer(map_client)
    fig = go.Figure(
        data=(
            go.Scattermapbox(
                lon=coord_x,
                lat=coord_y,
                mode='markers',
                marker=go.scattermapbox.Marker(
                    size=marker_size,
                    opacity=0.0,
                ),
                text='center',
                hoverinfo='lat+lon+text',
            ),
        ),
        layout=go.Layout(
            mapbox={
                'style': 'white-bg',
                'center': {
                    'lon': map_client.center()[1],
                    'lat': map_client.center()[0],
                },
                'zoom': map_client.default_zoom,
                'pitch': 0,
                'layers': [
                    {
                        'below': 'traces',
                        'sourcetype': 'raster',
                        'source': [map_tile_layer.url],
                    },
                ],
            },
            margin={
                'l': 30,
                'r': 30,
                'b': 30,
                't': 30,
            },
            autosize=False,
            width=window_width,
            height=window_height,
            paper_bgcolor='LightSteelBlue',
        ),
    )
    return fig
