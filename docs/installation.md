# PHGv2 - Installation

In this document, we will discuss the overall steps on how
to download and set up the PHGv2 package.

## Quick start
* Run on Linux or macOS (_Windows not currently tested_)
* Make sure you have $\geq$ Java 17
* Make sure you have [miniConda](https://docs.conda.io/projects/miniconda/en/latest/index.html#quick-command-line-install) installed
* Make sure you have the [libmamba]() solver installed:
  ```shell
  conda update -n base conda
  conda install -n base conda-libmamba-solver
  conda config --set solver libmamba
  ```
* Download and untar the [latest release](https://github.com/maize-genetics/phg_v2/releases/latest) of the PHGv2 package:
  ```shell
  curl -s https://api.github.com/repos/maize-genetics/phg_v2/releases/latest \
  | awk -F': ' '/browser_download_url/ && /\.tar/ {gsub(/"/, "", $(NF)); system("curl -LO " $(NF))}'
  ```
* Untar the package:
  ```shell
  tar -xvf <PHGv2_Release>.tar
  ```
  
* Navigate into uncompressed PHGv2 package:
  ```shell
  cd phg/bin
  ```
* Invoke PHGv2 through the `phg` wrapper script:
  ```shell
  ./phg --version
  ```
* Basic syntax is:
  ```shell
  phg [<options>] <command> [<args>]...
  ```
