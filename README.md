# Martin-Villalba Lab Dash App

This is the Martin-Villalba Lab's public-facing dash app source code.

# Configuration

The app is configured by the data it is provided and two environment variables:

- `DASH_TITLE` sets the title of the dash app and `<h2>`.
- `DASH_DEBUG` sets dash into debug mode. Only use this for development purposes.
- `DASH_DATA_DIR` is the path of the directory that contains the data and an optional `datasets.yaml`

## Data

The app reads all available `h5ad` files in `DASH_DATA_DIR` in a disk-backed mode.
This is done once at startup, if new data is added, the app currently needs to be restarted.

Additionally, it reads the contents of `datasets.yaml` to provide some pretty information about each dataset like reference and title.
If a dataset doesn't have metadata specified, a default title will be shown.

An example yaml looks like this:

```
- display: "10x"dataset
  file: "data.h5ad"
  reference_text: "People et al. (2015)"
  reference_uri: "https://doi.org/doi/doi"
- display: "Some other dataset"
  file: "cool_data.h5ad"
  reference_uri: "https://example.com"
  reference_text: "Example et al. (2010)"
```

This provides information on the data in `cool_data.h5ad` and `data.h5ad`.

Currently, only the keys `display`, `file` and `reference` are read for each dataset.
Other entries are ignored.
The name provided in `file` needs to match the filename in the `DASH_DATA_DIR`.

## Preparing data

To keep memory consumption low, I strongly advise only keeping the necessary parts of the `h5ad` files.
This currently only includes:

- `/X` (for expression values)
- Columns with a categorical datatype in `/obs/`
- The first two columns of all matrices in `/obsm/`
- `/var_names` and `/obs_names` (actually part of `/obs` and `/var` which are mandatory in `h5ad`)

Entries that should be deleted:

- `/layers`
- `/obsp`
- `/uns`
- `/varp`
- `/varm`
- All columns in `/var`
- All columns in `/obs` that aren't categorical.

Matrices with more than two columns should be deleted or chopped of at two columns.

## Data available

When data is loaded into the app, all categorical columns from `/obs` are available as grouping variables in the group plot and as colour variables for the embedding plot.
All embeddings in `/obsm` are available for the embedding plot. This plot takes the first two columns of any matrix in `/obsm`.
Further, all genes in `/X` can be used for colouring the embedding plot or plotting in the group plot.

# Running

## With containers

Start the app:

```
podman run -v /srv/dash-data:/data:ro,z -p 3838:8000 -dit --rm ghcr.io/lysogeny/a290dash:latest
```

- `-v` is necessary to map the directory containing your data (here `/srv/dash-data` to the container).
- `z` in volume is necessary for linux systems with an enforcing SELinux (e.g. CentOS).
- `-p` maps the local port to the container port. The first part is what your reverse proxy should be pointing at, the second part shouldn't change.

## Manual

Alternatively, if necessary requirements are present, the app can also be run with `gunicorn`.

```
gunicorn app:SERVER -b :8000
```

For development, `poetry` is used:

```
poetry run python app.py
DASH_DEBUG=1 poetry run python app.py # debug mode
```

## System requirements

Currently the app doesn't use much memory (<500MB), but this varies depending on the data it is provided.
For python dependencies see either the `pyproject.toml` or `requirements.txt`.
The docker container only needs a working docker.

## Adding data

With the docker setup described above, copy extra data into `/srv/dash-data`.
Optionally you may want to add some information about that data in the `datasets.yaml` file.
