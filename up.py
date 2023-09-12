"""upload."""

import dash
import dash_html_components as html
import dash_uploader as du

app = dash.Dash(__name__)

# 1) configure the upload folder
du.configure_upload(app, './upload/')

# 2) Use the Upload component
app.layout = html.Div([
    du.Upload(max_file_size=1024 * 10),
    du.Upload(max_file_size=1024 * 0.1),
])

if __name__ == '__main__':
    app.run_server(debug=True)
