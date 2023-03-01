from dash import Dash, dcc, html, Input, Output
from dash.exceptions import PreventUpdate

import scipy.sparse
import anndata

import plotly.express as px
import pandas as pd
import numpy as np

data = anndata.read_h5ad("data/data.h5ad")

class DataCollection:
    """Collection of datasets"""

    def __init__(self, data):
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

    def grouping_df(self, dataset_name, group_var):
        data = self.data[dataset_name]
        group_data = data.obs[[group_var]].rename(columns={group_var: "x"})
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
        catvars = catvars[(catvar_counts < 16) & (catvar_counts > 1)]
        return catvars

    def keys(self):
        return list(self.data.keys())


DATA = DataCollection({"SmartSeq3GLAST": data})

app = Dash(name=__name__, server=True)

app.layout = html.Div([
    html.H6("Data"),
    html.Div([
        dcc.Dropdown(DATA.keys(), DATA.keys()[0], id="dataset-name"),
        dcc.Dropdown(id="selected-gene-id"),
        dcc.Dropdown(id="selected-embedding"),
        dcc.Dropdown(id="selected-grouping-var"),
    ]),
    html.Br(),
    dcc.Graph(id="graph-umap"),
    dcc.Graph(id="graph-boxes"),
])

def search_decorator(search_function):
    def inner_function(function):
        def callback_function(dataset_name, search_value):
            if not search_value:
                raise PreventUpdate
            options = search_function(dataset_name)
            return function(search_value, options)
        return callback_function
    return inner_function



@app.callback(Output("selected-gene-id", "options"),
              Input("dataset-name", "value"),
              Input("selected-gene-id", "search_value"))
@search_decorator(lambda x: DATA.available_gene_ids(x))
def update_gene_options(search_value, options):
    return [o for o in options if search_value in o]

@app.callback(Output("selected-embedding", "options"),
              Input("dataset-name", "value"),
              Input("selected-embedding", "search_value"))
@search_decorator(lambda x: DATA.available_embedding_keys(x))
def update_embedding_options(search_value, options):
    return [o for o in options if search_value in o]

@app.callback(Output("selected-grouping-var", "options"),
              Input("dataset-name", "value"),
              Input("selected-grouping-var", "search_value"))
@search_decorator(lambda x: DATA.available_group_vars(x))
def update_group_options(search_value, options):
    return [o for o in options if search_value in o]

@app.callback(Output("graph-umap", "figure"),
              Input("dataset-name", "value"),
              Input("selected-gene-id", "value"),
              Input("selected-embedding", "value"))
def update_umap(dataset_name, gene_id, embedding_name):
    plot_data = DATA.embedding_df(dataset_name, embedding_name)
    if gene_id:
        plot_data = plot_data.join(DATA.gene_counts_df(dataset_name, gene_id))
    fig = px.scatter(x="x", y="y", color="Expression" if gene_id else None, data_frame=plot_data, template="simple_white")
    fig.update_xaxes(title_text="", showticklabels=False, tickvals=[])
    fig.update_yaxes(title_text="", showticklabels=False, tickvals=[],
                     scaleanchor="x", scaleratio=1)
    return fig

@app.callback(Output("graph-boxes", "figure"),
              Input("dataset-name", "value"),
              Input("selected-gene-id", "value"),
              Input("selected-grouping-var", "value"))
def update_boxplot(dataset_name, gene_id, group_var):
    plot_data = DATA.grouping_df(dataset_name, group_var)
    if gene_id:
        plot_data = plot_data.join(DATA.gene_counts_df(dataset_name, gene_id))
    else:
        plot_data = plot_data.join(DATA.gene_counts_dummy_df(dataset_name))
    fig = px.box(x="x", y="Expression", data_frame=plot_data, template="simple_white")
    fig.update_xaxes(title_text=group_var)
    return fig

if __name__ == "__main__":
    app.run(debug=True)
