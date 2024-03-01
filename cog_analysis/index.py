"""
Entry point for the application.
"""
import json
import dash
import dash_bootstrap_components as dbc

import layouts.cog_layout as cogl

app = dash.Dash(__name__, suppress_callback_exceptions=True, external_stylesheets=[dbc.themes.BOOTSTRAP])
server = app.server

with open("data/cog_category.json", "r", encoding="utf-8") as file:
    category_dict = json.load(file)

app.layout = cogl.cog_layout(category_dict)

if __name__ == '__main__':
    app.run_server(debug=True)
