import plotly.express as px
from keggx.keggx import KEGG
import matplotlib
import networkx as nx
import numpy as np
from collections import defaultdict
import networkx as nx
import plotly.graph_objects as go
from typing import Any, List, Dict, Tuple, Union, Callable, Set
from strokeDTI.target_identification.pathway_merge import (
    get_kegg_list,
    get_whole_graph,
    generate_graph_list,
)

import plotly.io as pio
import warnings
import pandas as pd

# Suppress specific pandas warnings
warnings.simplefilter(action="ignore", category=pd.errors.SettingWithCopyWarning)
warnings.simplefilter(action="ignore", category=FutureWarning)

DEFAULT_COLORS = [
    "#9B5DE5",
    "#F15BB5",
    "#FEE440",
    "#00BBF9",
    "#00F5D4",
    "#F7EDE2",
    "#F5CAC3",
    "#84A59D",
    "#F28482",
]
DUPLICATE_COLOR = "#FF9F1C"  # Updated duplicate color


def make_polar_plot(inputDf):
    fig = px.bar_polar(
        inputDf,
        r="log2(FC)",
        theta="Gene name",
        color="pathway",
        template="plotly_white",
        color_discrete_sequence=DEFAULT_COLORS,
        width=1200,
        height=1200,
    )

    fig.update_traces(marker=dict(line=dict(width=1, color="black")))
    fig.update_layout(showlegend=True)

    # fig.show()

    return fig


# Implementation modified from: https://gist.github.com/mogproject/50668d3ca60188c50e6ef3f5f3ace101
Vertex = Any
Edge = Tuple[Vertex, Vertex]
Num = Union[int, float]


class GraphVisualization:
    def __init__(
        self,
        G: nx.Graph,
        pos: Dict[Vertex, Union[Tuple[Num, Num], Tuple[Num, Num, Num]]],
        node_text: Union[Dict[Vertex, str], Callable] = None,
        node_text_position: Union[Dict[Vertex, str], Callable, str] = None,
        node_text_font_color: Union[Dict[Vertex, str], Callable, str] = None,
        node_text_font_family: Union[Dict[Vertex, str], Callable, str] = None,
        node_text_font_size: Union[Dict[Vertex, Num], Callable, str] = None,
        node_size: Union[Dict[Vertex, Num], Callable, Num] = None,
        node_color: Union[
            Dict[Vertex, Union[str, Num]], Callable, Union[str, Num]
        ] = None,
        node_border_width: Union[Dict[Vertex, Num], Callable, Num] = None,
        node_border_color: Union[Dict[Vertex, str], Callable, str] = None,
        node_opacity: Num = None,
        edge_width: Union[Dict[Edge, Num], Callable, Num] = None,
        edge_color: Union[Dict[Edge, str], Callable, str] = None,
        edge_opacity: Num = None,
    ):
        # check dimensions
        if all(len(pos.get(v, [])) == 2 for v in G):
            self.is_3d = False
        elif all(len(pos.get(v, [])) == 3 for v in G):
            self.is_3d = True
        else:
            raise ValueError

        # default settings
        self.default_settings = dict(
            node_text=str,  # show node label
            node_text_position="middle center",
            node_text_font_color="#FFFFFF",  #'#000000'
            node_text_font_family="Arial",
            node_text_font_size=1,  # 14
            node_size=8 if self.is_3d else 10,  # else 18
            node_color="#fcfcfc",
            node_border_width=2,
            node_border_color="#404040",
            node_opacity=1.0,  # 0.8
            edge_width=4 if self.is_3d else 1,  # 2
            edge_color="#9AA7B1",
            edge_opacity=0.8,  # 0.8
        )

        # save settings
        self.G = G
        self.pos = pos
        self.node_text = node_text
        self.node_text_position = node_text_position
        self.node_text_font_color = node_text_font_color
        self.node_text_font_family = node_text_font_family
        self.node_text_font_size = node_text_font_size
        self.node_size = node_size
        self.node_color = node_color
        self.node_border_width = node_border_width
        self.node_border_color = node_border_color
        self.node_opacity = node_opacity
        self.edge_width = edge_width
        self.edge_color = edge_color
        self.edge_opacity = edge_opacity

    def _get_edge_traces(self) -> List[Union[go.Scatter, go.Scatter3d]]:
        # group all edges by (color, width)
        groups = defaultdict(list)

        for edge in self.G.edges():
            color = self._get_setting("edge_color", edge)
            width = self._get_setting("edge_width", edge)
            groups[(color, width)] += [edge]

        # process each group
        traces = []
        for (color, width), edges in groups.items():
            x, y, z = [], [], []
            for v, u in edges:
                x += [self.pos[v][0], self.pos[u][0], None]
                y += [self.pos[v][1], self.pos[u][1], None]
                if self.is_3d:
                    z += [self.pos[v][2], self.pos[u][2], None]

            params = dict(
                x=x,
                y=y,
                mode="lines",
                hoverinfo="none",
                line=dict(color=color, width=width, shape="spline"),
                opacity=self._get_setting("edge_opacity"),
            )

            traces += [
                go.Scatter3d(z=z, **params) if self.is_3d else go.Scatter(**params)
            ]

        return traces

    def _get_node_trace(
        self, showlabel, colorscale, showscale, colorbar_title, reversescale
    ) -> Union[go.Scatter, go.Scatter3d]:
        x, y, z = [], [], []
        for v in self.G.nodes():
            x += [self.pos[v][0]]
            y += [self.pos[v][1]]
            if self.is_3d:
                z += [self.pos[v][2]]

        params = dict(
            x=x,
            y=y,
            mode="markers" + ("+text" if showlabel else ""),
            hoverinfo="text",
            marker=dict(
                showscale=showscale,
                colorscale=colorscale,
                reversescale=reversescale,
                color=self._get_setting("node_color"),
                size=self._get_setting("node_size"),
                line_width=self._get_setting("node_border_width"),
                line_color=self._get_setting("node_border_color"),
                colorbar=dict(
                    thickness=15,
                    title=colorbar_title,
                    xanchor="left",
                    titleside="right",
                ),
            ),
            text=self._get_setting("node_text"),
            textfont=dict(
                color=self._get_setting("node_text_font_color"),
                family=self._get_setting("node_text_font_family"),
                size=self._get_setting("node_text_font_size"),
            ),
            textposition=self._get_setting("node_text_position"),
            opacity=self._get_setting("node_opacity"),
        )

        trace = go.Scatter3d(z=z, **params) if self.is_3d else go.Scatter(**params)
        return trace

    def _get_setting(self, setting_name, edge=None):
        default_setting = self.default_settings.get(setting_name)
        def_func = (
            default_setting if callable(default_setting) else lambda x: default_setting
        )
        setting = self.__dict__.get(setting_name)

        if edge is None:  # vertex-specific
            if setting is None:  # default is used
                if callable(default_setting):  # default is a function
                    return [def_func(v) for v in self.G.nodes()]
                else:  # default is a constant
                    return default_setting
            elif callable(setting):  # setting is a function
                return [setting(v) for v in self.G.nodes()]
            elif isinstance(setting, dict):  # setting is a dict
                return [setting.get(v, def_func(v)) for v in self.G.nodes()]
            else:  # setting is a constant
                return setting
        else:  # edge-specific
            if setting is None:  # default is used
                return def_func(edge)
            elif callable(setting):  # setting is a function
                return setting(edge)
            elif isinstance(setting, dict):  # setting is a dict
                return setting.get(edge, def_func(edge))
            else:  # setting is a constant
                return setting

    def create_figure(
        self,
        showlabel=True,
        colorscale="YlGnBu",
        showscale=False,
        colorbar_title="",
        reversescale=False,
        **params
    ) -> go.Figure:
        axis_settings = dict(
            autorange=True,
            showgrid=False,
            zeroline=False,
            showline=False,
            visible=False,
            ticks="",
            showticklabels=False,
        )
        scene = dict(
            xaxis=axis_settings,
            yaxis=axis_settings,
            zaxis=axis_settings,
        )

        layout_params = dict(
            paper_bgcolor="rgba(255,255,255,255)",  # white
            plot_bgcolor="rgba(0,0,0,0)",  # transparent
            autosize=False,
            height=400,
            width=450 if showscale else 375,
            title="",
            titlefont_size=16,
            showlegend=False,
            hovermode="closest",
            margin=dict(b=5, l=0, r=0, t=20),
            annotations=[],
            xaxis=axis_settings,
            yaxis=axis_settings,
            scene=scene,
        )

        # override with the given parameters
        layout_params.update(params)

        # create figure
        fig = go.Figure(layout=go.Layout(**layout_params))
        fig.add_traces(self._get_edge_traces())
        fig.add_trace(
            self._get_node_trace(
                showlabel, colorscale, showscale, colorbar_title, reversescale
            )
        )
        return fig


# Extend color_scale if there are more pathways than colors
def get_color_scale(num_colors: int) -> List[str]:
    if num_colors <= len(DEFAULT_COLORS):
        return DEFAULT_COLORS[:num_colors]
    else:
        # Repeat colors if there are more pathways than colors
        return (DEFAULT_COLORS * (num_colors // len(DEFAULT_COLORS) + 1))[:num_colors]


def get_kegg_data(kegg_path_list: List[str], species: str) -> List:
    """Fetch KEGG list based on species."""
    return get_kegg_list(kegg_path_list, species)


def identify_duplicates(
    graphs: List[nx.Graph],
) -> Tuple[Set[str], Set[Tuple[str, str]]]:
    """Identify duplicate nodes and edges across multiple graphs."""
    unique_nodes = set()
    duplicate_nodes = set()
    unique_edges = set()
    duplicate_edges = set()

    for graph in graphs:
        # Nodes
        for node in graph.nodes:
            if node in unique_nodes:
                duplicate_nodes.add(node)
            else:
                unique_nodes.add(node)

        # Edges
        for edge in graph.edges:
            if edge in unique_edges:
                duplicate_edges.add(edge)
            else:
                unique_edges.add(edge)

    return duplicate_nodes, duplicate_edges


def compose_whole_graph(kegg_path_list: List[str], species: str) -> nx.Graph:
    """Compose a whole graph from the list of KEGG paths."""
    whole_graphs = [
        get_whole_graph(lst) for lst in get_kegg_data(kegg_path_list, species)
    ]
    return nx.compose_all(whole_graphs)


def sort_nodes(whole_graph: nx.Graph, subgraphs: List[nx.Graph]) -> List[str]:
    """Sort nodes by separating those not in subgraphs and then appending the rest."""
    subgraph_nodes = set(node for graph in subgraphs for node in graph.nodes())
    sorted_nodes = [node for node in whole_graph.nodes() if node not in subgraph_nodes]
    sorted_nodes.extend(subgraph_nodes)
    return sorted_nodes


def assign_node_attributes(
    graphs: List[nx.Graph], color_scale: List[str], duplicate_nodes: Set[str]
) -> Tuple[Dict[str, str], Dict[str, float], Dict[str, int], Dict[str, str]]:
    """Assign colors, sizes, and label attributes to nodes."""
    node_color = {}
    node_size = {}
    node_label_font = {}
    node_label_color = {}

    for graph, color in zip(graphs, color_scale):
        for node in graph.nodes:
            node_color[node] = color
            node_size[node] = 45 + graph.nodes[node].get("degree_centrality", 0) * 80
            node_label_font[node] = 14
            node_label_color[node] = "#000000"

    # Update colors for duplicate nodes
    for node in duplicate_nodes:
        node_color[node] = DUPLICATE_COLOR

    return node_color, node_size, node_label_font, node_label_color


def assign_edge_attributes(
    graphs: List[nx.Graph],
    color_scale: List[str],
    duplicate_edges: Set[Tuple[str, str]],
) -> Tuple[
    Dict[Tuple[str, str], str],
    Dict[Tuple[str, str], float],
    Dict[Tuple[str, str], float],
]:
    """Assign colors, widths, and opacity to edges."""
    edge_color = {}
    edge_width = {}
    edge_opacity = {}

    for graph, color in zip(graphs, color_scale):
        for edge in graph.edges:
            edge_color[edge] = color
            edge_width[edge] = 6
            edge_opacity[edge] = 1.0

    # Update colors for duplicate edges
    for edge in duplicate_edges:
        edge_color[edge] = DUPLICATE_COLOR

    return edge_color, edge_width, edge_opacity


def create_visualization(
    graph: nx.Graph,
    pos: Dict[str, Tuple[float, float]],
    node_size: Dict[str, float],
    node_color: Dict[str, str],
    node_label_font: Dict[str, int],
    node_label_color: Dict[str, str],
    edge_color: Dict[Tuple[str, str], str],
    edge_width: Dict[Tuple[str, str], float],
    edge_opacity: float = 1.0,
    node_opacity: float = 1.0,
    border_width: float = 1.5,
) -> go.Figure:
    """Create and return a Plotly figure for the graph."""
    vis = GraphVisualization(
        graph,
        pos,
        node_size=node_size,
        node_border_width=border_width,
        node_text_font_size=node_label_font,
        node_text_font_color=node_label_color,
        edge_color=edge_color,
        edge_width=edge_width,
        edge_opacity=edge_opacity,
        node_color=node_color,
        node_opacity=node_opacity,
    )
    return vis.create_figure(height=1000, width=1200, showlabel=True)


def plot_graph(
    root_path: str,
    kegg_path_list: List[str],
    species: str,
    sequence_path: str,
    color_scale: List[str] = DEFAULT_COLORS,
):
    """Main function to plot the graph based on KEGG paths and species."""
    # Generate graphs and pathway names
    graphs, plst_name = generate_graph_list(
        root_path, kegg_path_list, species, sequence_path
    )

    # Update color_scale based on number of pathways
    color_scale = get_color_scale(len(graphs))

    # Identify duplicates
    duplicate_nodes, duplicate_edges = identify_duplicates(graphs)

    # Compose the whole graph
    whole_graph = compose_whole_graph(kegg_path_list, species)

    # Sort nodes
    sorted_nodes = sort_nodes(whole_graph, graphs)
    sorted_subgraph = whole_graph.subgraph(sorted_nodes).copy()

    # Assign node attributes
    node_color, node_size, node_label_font, node_label_color = assign_node_attributes(
        graphs, color_scale, duplicate_nodes
    )

    # Assign edge attributes
    edge_color, edge_width, edge_opacity = assign_edge_attributes(
        graphs, color_scale, duplicate_edges
    )

    # Define layout position
    pos = nx.spring_layout(sorted_subgraph, k=0.75)

    # Create visualization
    fig = create_visualization(
        sorted_subgraph,
        pos,
        node_size,
        node_color,
        node_label_font,
        node_label_color,
        edge_color,
        edge_width,
    )

    for trace in fig.data:
        trace.showlegend = False

    # Add dummy traces for legend
    for pathway_name, color in zip(plst_name, color_scale):
        fig.add_trace(
            go.Scatter(
                x=[None],
                y=[None],
                mode="markers",
                marker=dict(size=10, color=color),
                name=pathway_name,
            )
        )

    # Optionally, handle duplicates in the legend
    fig.add_trace(
        go.Scatter(
            x=[None],
            y=[None],
            mode="markers",
            marker=dict(size=10, color=DUPLICATE_COLOR),
            name="Overlapping Nodes/Edges",
        )
    )

    # Update layout to show legend
    fig.update_layout(showlegend=True)

    # fig.show()

    image_path = root_path + "graph2.png"
    pio.write_image(fig, image_path, scale=2)

    cellDeathGraph = nx.compose_all(graphs)
    nx.write_graphml(cellDeathGraph, root_path + "graph.net")
