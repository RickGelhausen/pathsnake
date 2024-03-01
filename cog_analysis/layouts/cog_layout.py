"""
Layout for the overview page.
"""

from dash import html, dcc, dash_table
import dash_bootstrap_components as dbc
from plotly.colors import DEFAULT_PLOTLY_COLORS


from lib.cog_plots import plot_category_size_graph, plot_category_start_graph, plot_category_stop_graph

COG_CATEGORIES = {
    "A": "RNA processing and modification",
    "B": "Chromatin structure and dynamics",
    "C": "Energy production and conversion",
    "D": "Cell cycle control, cell division, chromosome partitioning",
    "E": "Amino acid transport and metabolism",
    "F": "Nucleotide transport and metabolism",
    "G": "Carbohydrate transport and metabolism",
    "H": "Coenzyme transport and metabolism",
    "I": "Lipid transport and metabolism",
    "J": "Translation, ribosomal structure and biogenesis",
    "K": "Transcription",
    "L": "Replication, recombination and repair",
    "M": "Cell wall/membrane/envelope biogenesis",
    "N": "Cell motility",
    "O": "Post-translational modification, protein turnover, chaperones",
    "P": "Inorganic ion transport and metabolism",
    "Q": "Secondary metabolites biosynthesis, transport, and catabolism",
    "R": "General function prediction only",
    "S": "Function unknown",
    "T": "Signal transduction mechanisms",
    "U": "Intracellular trafficking, secretion, and vesicular transport",
    "V": "Defense mechanisms",
    "W": "Extracellular structures",
    "X": "Mobilome: prophages, transposons",
    "Y": "Nuclear structure",
    "Z": "Cytoskeleton",
}


def cog_layout(category_dict):
    """
    Data for the cog overview page.
    """

    return html.Div(
        children=[

            dbc.Container(
                fluid=True,
                children=[
                    html.H1(
                        "Cluster of orthologous groups (COG) analysis",
                        className="my-3 text-center",
                    ),
                    dbc.Row(
                        className="mt-3",
                        children=[
                            dbc.Col(
                                className="col-6",
                                children=[
                                    dcc.Graph(id="cog-start-graph", figure=plot_category_start_graph(category_dict), config={
                                                    "toImageButtonOptions": {
                                                        "format": "svg",
                                                        "filename": "start_category_graph",
                                                    }
                                                })
                                ],
                            ),
                            dbc.Col(
                                className="col-6",
                                children=[dcc.Graph(id="cog-stop-graph", figure=plot_category_stop_graph(category_dict), config={
                                                    "toImageButtonOptions": {
                                                        "format": "svg",
                                                        "filename": "stop_category_graph",
                                                    }

                                })],
                            ),
                        ],
                    ),
                    dbc.Row(
                        children=[
                            dbc.Col(
                                className="col-6",
                                children=[dcc.Graph(id="cog-size-graph", figure=plot_category_size_graph(category_dict), config={
                                                    "toImageButtonOptions": {
                                                        "format": "svg",
                                                        "filename": "size_category_graph",
                                                    }
                                 })],
                            ),
                            dbc.Col(
                                className="col-6",
                                children=[
                                    dash_table.DataTable(
                                        id="cog-description-table",
                                        columns=[
                                            {"name": "Category", "id": "Category"},
                                            {
                                                "name": "Description",
                                                "id": "Description",
                                            },
                                        ],
                                        data=[
                                            {"Category": key, "Description": value}
                                            for key, value in COG_CATEGORIES.items()
                                        ],
                                        cell_selectable=False,

                                        style_cell={
                                            "font-family": "Roboto",
                                            "color": "black",
                                            "text-align": "left",
                                            "font-size": "14px",
                                        },
                                        style_header={
                                            "font-family": "Roboto",
                                            "font-weight": "bold",
                                        },
                                        style_header_conditional=[
                                            {
                                                "if": {"column_id": "Category"},
                                                "textAlign": "center",
                                            }
                                        ],
                                        style_data_conditional=[
                                            {
                                                "if": {"column_id": "Category"},
                                                "textAlign": "center",
                                            }
                                        ],
                                    )
                                ],
                            ),
                        ]
                    ),
                ],
            ),
        ]
    )