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

    meta_default = {
        "reference_uri": "",
        "reference_text": "",
        "display": "",
    }

    def __init__(self, dir):
        meta_file = f"{dir}/datasets.yaml"
        self.load_datasets(dir)
        self.load_metadata(meta_file)

    def load_datasets(self, dir):
        logging.info(f"reading files from `{DASH_DATA_DIR}`")
        datasets = list(filter(lambda x: x.endswith(".h5ad"), os.listdir(dir)))
        logging.info(f"found {len(datasets)} files:")
        for line in datasets:
            logging.info(f"File `{line}`")
        self.data = {dataset.split(".")[0]: anndata.read_h5ad(f"{dir}/{dataset}", backed="r") for dataset in datasets}
        logging.info(f"loaded {len(self.data.keys())} datasets.")

    def load_metadata(self, meta_file):
        logging.info(f"reading meta file from `{DASH_DATA_DIR}`")
        if os.path.exists(meta_file):
            logging.info(f"Including metadata from {meta_file}")
            with open(meta_file) as connection:
                meta = yaml.safe_load(connection)
            meta = {d["file"].split(".")[0]: d for d in meta}
        else:
            meta = {}
        self.meta = meta

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
    
    def meta_value(self, dataset_name, key):
        meta = self.meta.get(dataset_name, self.meta_default)
        return meta.get(key, "")

    def available_gene_ids(self, dataset_name):
        return self.data[dataset_name].var_names

    def available_embedding_keys(self, dataset_name):
        data = self.data[dataset_name]
        return list(data.obsm.keys())

    def available_group_vars(self, dataset_name):
        data = self.data[dataset_name]
        catvars = data.obs.columns[data.obs.dtypes == "category"]
        catvar_counts = pd.Series([len(data.obs[v].cat.categories) for v in catvars])
        catvars = list(catvars[(catvar_counts < 24) & (catvar_counts > 1)])
        boolvars = list(data.obs.columns[data.obs.dtypes == "bool"])
        return [x for x in data.obs.columns if x in (catvars + boolvars)]

    def keys(self):
        return list(self.data.keys())

    def keys_dicts(self):
        return [{"label": self.meta_value(x, "display"), "value": x} for x in self.meta.keys()]


if "DASH_DEBUG" in os.environ:
    logging.basicConfig(level=logging.DEBUG)
else:
    logging.basicConfig(level=logging.INFO)

if "DASH_DATA_DIR" in os.environ:
    DASH_DATA_DIR = os.environ["DASH_DATA_DIR"]
else:
    DASH_DATA_DIR = "data"
if "DASH_TITLE" in os.environ:
    DASH_TITLE = os.environ["DASH_TITLE"]
else:
    DASH_TITLE = "Martin-Villalba Lab Data Explorer"

if "DASH_DEBUG" in os.environ:
    DASH_TITLE = "(DEBUG) " + DASH_TITLE

DATA = DataCollection(DASH_DATA_DIR)

APP = Dash(name=__name__, server=True, title=DASH_TITLE)

APP.layout = html.Div([
    html.Header([
        html.H2(DASH_TITLE, id="data-name"),
        html.Div("", id="data-info"),
        html.Div([
            dcc.Dropdown(DATA.keys_dicts(), DATA.keys()[0], id="dataset-name", className="", placeholder="Dataset..."), 
        ], className=""),
    ], className="flex-col-top"),
    html.Section([
        html.Div([
            html.Div([
                dcc.Dropdown(id="selected-embedding", placeholder="Embedding..."),
                dcc.Dropdown(["gene", "categorical"], "gene", id="selected-embedding-source", placeholder="Embedding colour source..."),
                dcc.Dropdown(id="selected-embedding-var", placeholder="Embedding colour variable...")
            ], className="flex-col-top flex-row-container three-elem-input"),
            dcc.Graph(id="graph-umap", className="flex-col-mid", responsive=True),
        ], className="flex-row-element flex-col-container"), 
        html.Div([
            html.Div([
                dcc.Dropdown(id="selected-grouping-var", placeholder="Grouping variable...", multi=True),
                dcc.Dropdown(id="selected-gene-id", placeholder="Gene Name..."),
            ], className="flex-col-top flex-row-container two-elem-input"),
            dcc.Graph(id="graph-boxes", className="flex-col-mid", responsive=True),
        ], className="flex-row-element flex-col-container"), 
    ], className="flex-col-mid flex-row-container"),
    html.Footer(["DKFZ/A290 Dash App"], className="flex-col-bot"),
], className="flex-col-container root")

# Callback to update dataset info from metadata.
@APP.callback(Output("data-info", "children"),
              Input("dataset-name", "value"))
def update_info(dataset_name):
    ref_uri = DATA.meta_value(dataset_name, "reference_uri")
    ref_text = DATA.meta_value(dataset_name, "reference_text")
    display = DATA.meta_value(dataset_name, "display")
    return html.P([
        f"{display} | ",
        html.A(f"{ref_text}", href=ref_uri),
    ])

# Callbacks to update dropdown values if the dataset changed.
# Check if the available values are in the new dataset, if not, set something else
# For some reason this is not necessary for the grouping variable dropdown, perhaps because it allows multiple selections?
@APP.callback(Output("selected-embedding-var", "value"),
              Input("dataset-name", "value"),
              Input("selected-embedding-var", "value"))
def update_embedding_var(dataset_name, selected):
    if selected in DATA.available_group_vars(dataset_name):
        return selected
    return ""

@APP.callback(Output("selected-gene-id", "value"),
              Input("dataset-name", "value"),
              Input("selected-gene-id", "value"))
def update_gene_id(dataset_name, selected):
    if selected in DATA.available_gene_ids(dataset_name):
        return selected
    return ""

@APP.callback(Output("selected-embedding", "value"),
              Input("dataset-name", "value"),
              Input("selected-embedding", "value"))
def update_embedding_value(dataset_name, selected):
    if selected in DATA.available_embedding_keys(dataset_name):
        return selected
    return ""

# Callbacks to update dropdown options if the dataset changed.
@APP.callback(Output("selected-gene-id", "options"),
              Input("dataset-name", "value"),
              Input("selected-gene-id", "search_value"))
def update_gene_options(dataset_name, search_value):
    if not search_value:
        raise PreventUpdate
    options = DATA.available_gene_ids(dataset_name)
    return [o for o in options if search_value.lower() in o.lower()]

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

# Callbacks to update the plots if necessary fields changed.
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
