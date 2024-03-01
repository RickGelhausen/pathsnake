"""
Module for plotting functions.
"""

import pandas as pd
import plotly.graph_objs as go
from plotly.colors import DEFAULT_PLOTLY_COLORS
from plotly.subplots import make_subplots


# Set a modern style for the plots
PLOT_STYLE = {
    'font_family': 'Roboto, sans-serif',
    'title_font_size': 16,
    'axis_label_font_size': 14,
    'legend_title_font_size': 12,
    'grid_color': '#e9ecef',
    'axis_color': '#343a40'
}


def plot_category_size_graph(category_dict):
    """
    Plot the category size graph as a subplot with pie chart and bar plot.
    """

    # Calculate the size of each category
    category_size = {}
    for category, codons_dict in category_dict["stop"].items():
        category_size[category] = sum(codons_dict.values())

    # Sorting data
    data = zip(category_size.keys(), category_size.values())
    data = sorted(data, key=lambda x: x[1], reverse=True)

    # Color list
    color_list = ['#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#000000', '#006400', '#d7f9fa', '#8c7088']

    # Unpacking sorted data
    category, size = zip(*data)

    # Create subplots: 1 row, 2 columns
    fig = make_subplots(rows=1, cols=2, specs=[[{'type':'domain'}, {'type':'bar'}]])

    # Adding Pie chart
    fig.add_trace(go.Pie(labels=category, values=size, marker=dict(colors=color_list), textinfo='value'), row=1, col=1)

    # Adding Bar chart
    fig.add_trace(go.Bar(x=category, y=size, marker_color=color_list, showlegend=False), row=1, col=2)

    fig.update_xaxes(tickangle=0, tick0=0, dtick=1, row=1, col=2)
    # Update the layout with the custom style
    fig.update_layout(
        font_family=PLOT_STYLE['font_family'],
        title_text="Category Size Analysis",
        title_font_size=PLOT_STYLE['title_font_size'],
        legend_title_font_size=PLOT_STYLE['legend_title_font_size'],
        plot_bgcolor='white',
        paper_bgcolor='white',
        xaxis=dict(
            showgrid=True,
            gridcolor=PLOT_STYLE['grid_color'],
            color=PLOT_STYLE['axis_color'],
            title_font_size=PLOT_STYLE['axis_label_font_size']
        ),
        yaxis=dict(
            showgrid=True,
            gridcolor=PLOT_STYLE['grid_color'],
            color=PLOT_STYLE['axis_color'],
            title_font_size=PLOT_STYLE['axis_label_font_size']
        )
    )
    # Update layout if needed
    fig.update_layout(title_text="Category Size Analysis")

    return fig

def plot_category_start_graph(category_dict):
    """
    Plot the category start graph showing for each category the amount of start codons of each type
    """


    processed_data = {}

    for category, codons in category_dict['start'].items():
        total = sum(codons.values())
        processed_data[category] = {
            'ATG': (codons.get('ATG', 0) / total) * 100,
            'TTG': (codons.get('TTG', 0) / total) * 100,
            'GTG': (codons.get('GTG', 0) / total) * 100,
            'Others': sum((val / total) * 100 for codon, val in codons.items() if codon not in ['ATG', 'TTG', 'GTG'])
        }

    # Calculate the size of each category
    category_size = {}
    for category, codons_dict in category_dict["start"].items():
        category_size[category] = sum(codons_dict.values())

    # Convert to DataFrame for easier plotting
    df = pd.DataFrame.from_dict(processed_data, orient='index')

    # Add a new column for category sizes
    df['Size'] = df.index.map(category_size)

    # Sort the DataFrame based on the Size column
    df = df.sort_values(by='Size', ascending=False)
    df = df.drop(columns=['Size'])

    fig = go.Figure()
    # Add traces for each codon type
    for codon in df.columns:
        fig.add_trace(go.Bar(
            name=codon,
            x=df.index,
            y=df[codon],
            hoverinfo='y',
            text=df[codon].round(1),
            textposition='inside'
        ))

    # Update layout
    fig.update_layout(
        barmode='stack',
        title='Percentage of Start Codons in Each Category',
        xaxis_title='Category',
        yaxis_title='Start Codon Count (%)',
        yaxis=dict(
            ticksuffix="%"
        ),
        legend_title_text='Codon',
    )

    fig.update_layout(
        font_family=PLOT_STYLE['font_family'],
        title_font_size=PLOT_STYLE['title_font_size'],
        legend_title_font_size=PLOT_STYLE['legend_title_font_size'],
        plot_bgcolor='white',
        paper_bgcolor='white',
        xaxis=dict(
            showgrid=True,
            gridcolor=PLOT_STYLE['grid_color'],
            color=PLOT_STYLE['axis_color'],
            title_font_size=PLOT_STYLE['axis_label_font_size']
        ),
        yaxis=dict(
            showgrid=True,
            gridcolor=PLOT_STYLE['grid_color'],
            color=PLOT_STYLE['axis_color'],
            title_font_size=PLOT_STYLE['axis_label_font_size']
        )
    )
    return fig


def plot_category_stop_graph(category_dict):
    """
    Plot the category start graph showing for each category the amount of start codons of each type
    """


    processed_data = {}

    for category, codons in category_dict['stop'].items():
        total = sum(codons.values())
        processed_data[category] = {
            'TAA': (codons.get('TAA', 0) / total) * 100,
            'TAG': (codons.get('TAG', 0) / total) * 100,
            'TGA': (codons.get('TGA', 0) / total) * 100,
            'Others': sum((val / total) * 100 for codon, val in codons.items() if codon not in ['TAA', 'TAG', 'TGA'])
        }

    # Calculate the size of each category
    category_size = {}
    for category, codons_dict in category_dict["stop"].items():
        category_size[category] = sum(codons_dict.values())

    # Convert to DataFrame for easier plotting
    df = pd.DataFrame.from_dict(processed_data, orient='index')

    # Add a new column for category sizes
    df['Size'] = df.index.map(category_size)

    # Sort the DataFrame based on the Size column
    df = df.sort_values(by='Size', ascending=False)
    df = df.drop(columns=['Size'])

    fig = go.Figure()

    # Add traces for each codon type
    for codon in df.columns:
        fig.add_trace(go.Bar(
            name=codon,
            x=df.index,
            y=df[codon],
            hoverinfo='y',
            text=df[codon].round(1),
            textposition='inside'
        ))

    # Update layout
    fig.update_layout(
        barmode='stack',
        title='Percentage of Stop Codons in Each Category',
        xaxis_title='Category',
        yaxis_title='Stop Codon Count (%)',
        yaxis=dict(
            ticksuffix="%"
        ),
        legend_title_text='Codon'
    )

    fig.update_layout(
        font_family=PLOT_STYLE['font_family'],
        title_font_size=PLOT_STYLE['title_font_size'],
        legend_title_font_size=PLOT_STYLE['legend_title_font_size'],
        plot_bgcolor='white',
        paper_bgcolor='white',
        xaxis=dict(
            showgrid=True,
            gridcolor=PLOT_STYLE['grid_color'],
            color=PLOT_STYLE['axis_color'],
            title_font_size=PLOT_STYLE['axis_label_font_size']
        ),
        yaxis=dict(
            showgrid=True,
            gridcolor=PLOT_STYLE['grid_color'],
            color=PLOT_STYLE['axis_color'],
            title_font_size=PLOT_STYLE['axis_label_font_size']
        )
    )

    return fig