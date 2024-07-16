# Static Site Setup
These instructions are for deploying a local version of the static PHG website:

## Retrieve `mkdocs-material` from conda

```
# create optional environment
conda create --name site-generator
conda activate site-generator

# install mkdocs-material
conda install conda-forge::mkdocs-material
```

## Launch site

```
# ensure you are in root of project directory
cd phg_v2

mkdocs serve
```

