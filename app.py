import os
import logging

import yaml

from dash import Dash, dcc, html, Input, Output
from dash.exceptions import PreventUpdate

import scipy.sparse
import anndata

import plotly.express as px
import pandas as pd
import numpy as np


class DataCollection:
    """Collection of datasets"""

    def __init__(self, data: dict):
        self.data = data

    def gene_counts_df(self, dataset_name, gene_id):
        data = self.data[dataset_name]
        count_data = data[:,gene_id].to_df().rename(columns={gene_id: "Expression"})
        return count_data

    def gene_counts_dummy_df(self, dataset_name):
        data = self.data[dataset_name]
        dummy =  pd.DataFrame(np.zeros(len(data.obs_names)), index=data.obs_names, columns=["Expression"])
        return dummy 

    def embedding_df(self, dataset_name, embedding_name):
        data = self.data[dataset_name]
        plot_data = data.obsm[embedding_name]
        if isinstance(plot_data, scipy.sparse.spmatrix):
            plot_data = plot_data.to_dense()
        plot_data = pd.DataFrame(plot_data[:,0:2], index=data.obs_names, columns=["x", "y"])
        return plot_data

    def grouping_df(self, dataset_name, group_vars):
        data = self.data[dataset_name]
        group_data = data.obs[group_vars]
        return group_data

    def available_gene_ids(self, dataset_name):
        return self.data[dataset_name].var_names

    def available_embedding_keys(self, dataset_name):
        data = self.data[dataset_name]
        return list(data.obsm.keys())

    def available_group_vars(self, dataset_name):
        data = self.data[dataset_name]
        catvars = data.obs.columns[data.obs.dtypes == "category"]
        catvar_counts = pd.Series([len(data.obs[v].cat.categories) for v in catvars])
        catvars = catvars[(catvar_counts < 24) & (catvar_counts > 1)]
        return catvars

    def keys(self):
        return list(self.data.keys())

if "DASH_DEBUG" in os.environ:
    logging.basicConfig(level=logging.DEBUG)
else:
    logging.basicConfig(level=logging.INFO)

if "DASH_DATA_DIR" in os.environ:
    DASH_DATA_DIR = os.environ["DASH_DATA_DIR"]
else:
    DASH_DATA_DIR = "data"
DASH_META_FILE = f"{DASH_DATA_DIR}/datasets.yaml"

logging.info(f"reading files from `{DASH_DATA_DIR}`")

DATASETS = list(filter(lambda x: x.endswith(".h5ad"), os.listdir(DASH_DATA_DIR)))

logging.info(f"found {len(DATASETS)} files:")
for line in DATASETS:
    logging.info(f"File `{line}`")

DATA = DataCollection({dataset.split(".")[0]: anndata.read_h5ad(f"{DASH_DATA_DIR}/{dataset}", backed="r") for dataset in DATASETS})

logging.info(f"loaded {len(DATA.keys())} datasets.")

META_DEFAULT = {
    "filename": "",
    "reference": "",
    "display": "A290 Molecular Neurobiology Data Viewer",
}

if os.path.exists("data/datasets.yaml"):
    with open("data/datasets.yaml") as connection:
        DATAMETA = yaml.safe_load(connection)
    DATAMETA = {d["file"].split(".")[0]: d for d in DATAMETA}
    logging.info(f"Including metadata from {DASH_META_FILE}")
else:
    DATAMETA = {d: META_DEFAULT for d in DATA.keys()}
    logging.info(f"Created default metadata")


APP = Dash(name=__name__, server=True)

APP.layout = html.Div([
    html.Header([
        html.H2("A290 Molecular Neurobiology Data Viewer", id="data-name"),
        html.Div("", id="data-info"),
        html.Div([
            dcc.Dropdown(DATA.keys(), DATA.keys()[0], id="dataset-name", className="", placeholder="Dataset..."), 
        ], className=""),
    ], className="row"),
    html.Section([
        html.Div([
            html.Div([
                dcc.Dropdown(id="selected-embedding", className="flex-child", placeholder="Embedding..."),
                dcc.Dropdown(["gene", "categorical"], id="selected-embedding-source", className="flex-child", placeholder="Embedding colour source..."),
                dcc.Dropdown(id="selected-embedding-var", className="flex-child", placeholder="Embedding colour variable...")
            ], className="flex-container"),
            dcc.Graph(id="graph-umap", className="graph-output"),
        ], className="flex-child"), 
        html.Div([
            html.Div([
                dcc.Dropdown(id="selected-grouping-var", className="flex-child", placeholder="Grouping variable...", multi=True),
                dcc.Dropdown(id="selected-gene-id", className="flex-child", placeholder="Gene Name..."),
            ], className="flex-container"),
            dcc.Graph(id="graph-boxes", className="graph-output"),
        ], className="flex-child"), 
    ], className="flex-container"),
    html.Footer(["DKFZ/A290 Dash App"], className="row"),
], className="box")

@APP.callback(Output("data-info", "children"),
              Input("dataset-name", "value"))
def update_info(dataset_name):
    content = DATAMETA[dataset_name]["reference"]
    display = DATAMETA[dataset_name]["display"]
    return html.P([
        f"{display}: ",
        html.A(f"{content}", href=content),
    ])

@APP.callback(Output("selected-gene-id", "options"),
              Input("dataset-name", "value"),
              Input("selected-gene-id", "search_value"))
def update_gene_options(dataset_name, search_value):
    if not search_value:
        raise PreventUpdate
    options = DATA.available_gene_ids(dataset_name)
    return [o for o in options if search_value in o]

@APP.callback(Output("selected-embedding", "options"),
              Input("dataset-name", "value"))
def update_embedding_options(dataset_name):
    return DATA.available_embedding_keys(dataset_name)

@APP.callback(Output("selected-grouping-var", "options"),
              Input("dataset-name", "value"))
def update_group_options(dataset_name):
    return DATA.available_group_vars(dataset_name)

@APP.callback(Output("selected-embedding-var", "options"),
              Input("dataset-name", "value"))
def update_embedding_cat_options(dataset_name):
    return DATA.available_group_vars(dataset_name)

@APP.callback(Output("graph-umap", "figure"),
              Input("dataset-name", "value"),
              Input("selected-gene-id", "value"),
              Input("selected-embedding", "value"),
              Input("selected-embedding-var", "value"),
              Input("selected-embedding-source", "value"))
def update_umap(dataset_name, gene_id, embedding_name, group_var, colour_source):
    if not dataset_name or not embedding_name:
        return px.scatter(template="simple_white")
    plot_data = DATA.embedding_df(dataset_name, embedding_name)
    if colour_source == "gene" and gene_id:
        plot_data = plot_data.join(DATA.gene_counts_df(dataset_name, gene_id))
        fig = px.scatter(x="x", y="y", color="Expression" if gene_id else None, 
                         data_frame=plot_data, template="simple_white")
        fig.update_xaxes(title_text="", showticklabels=False, tickvals=[])
        fig.update_yaxes(title_text="", showticklabels=False, tickvals=[],
                         scaleanchor="x", scaleratio=1)
    elif colour_source == "categorical" and group_var:
        plot_data = plot_data.join(DATA.grouping_df(dataset_name, group_var))
        fig = px.scatter(x="x", y="y", color=group_var if group_var else None, 
                         data_frame=plot_data, template="simple_white")
    else:
        fig =  px.scatter(x="x", y="y", data_frame=plot_data,
                          template="simple_white")
    fig.update_xaxes(title_text="", showticklabels=False, tickvals=[])
    fig.update_yaxes(title_text="", showticklabels=False, tickvals=[],
                     scaleanchor="x", scaleratio=1)
    return fig

@APP.callback(Output("graph-boxes", "figure"),
              Input("dataset-name", "value"),
              Input("selected-gene-id", "value"),
              Input("selected-grouping-var", "value"))
def update_boxplot(dataset_name, gene_id, group_vars):
    if not dataset_name or not gene_id or not group_vars:
        return px.scatter(template="simple_white")
    if len(group_vars) > 3:
        group_vars = group_vars[0:3]
    group_vars = [v for v in group_vars if v in DATA.available_group_vars(dataset_name)]
    plot_data = DATA.grouping_df(dataset_name, group_vars)
    if gene_id:
        plot_data = plot_data.join(DATA.gene_counts_df(dataset_name, gene_id))
    else:
        plot_data = plot_data.join(DATA.gene_counts_dummy_df(dataset_name))
    fig = px.box(y="Expression", 
                 x=group_vars[0], 
                 color=group_vars[1] if len(group_vars) > 1 else None,
                 facet_col=group_vars[2] if len(group_vars) > 2 else None,
                 data_frame=plot_data, 
                 template="simple_white")
    return fig

print("App ready!")

if __name__ == "__main__":
    APP.run(debug=True if "DASH_DEBUG" in os.environ else False)

SERVER = APP.server
