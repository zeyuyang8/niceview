"""App."""

import dash
from dash import html
from dash import dcc
from dash.dependencies import Input, Output, State


# Initialize app
app = dash.Dash(
    __name__,
    meta_tags=[
        {'name': 'viewport'},
    ],
    suppress_callback_exceptions=True,
)

app.title = 'BioGIS'
server = app.server


if __name__ == '__main__':
    app.run_server(host='localhost', debug=True, dev_tools_hot_reload=True)
