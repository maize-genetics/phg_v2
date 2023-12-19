# PHGv2 - Installation

In this document, we will discuss the overall steps on how
to download and set up the PHGv2 package.

## Quick start
* Run on Linux (_Windows not currently tested_)
* Make sure you have $\geq$ Java 17
* Make sure you have [Miniconda](https://docs.conda.io/projects/miniconda/en/latest/index.html#quick-command-line-install) installed
* Make sure you have the [libmamba]() solver installed:
  ```shell
  conda update -n base conda
  conda install -n base conda-libmamba-solver
  conda config --set solver libmamba
  ```
* Download the [latest release](https://github.com/maize-genetics/phg_v2/releases/latest) of the PHGv2 package:
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

## Requirements
PHGv2 requires basic software components: a Unix-style operating 
system and Java version 17 or higher. PHGv2 
also relies on external programs for alignment and storage, including 
[AnchorWave](https://github.com/baoxingsong/AnchorWave) and 
[TileDB-VCF](https://docs.tiledb.com/main/integrations-and-extensions/genomics/population-genomics). 
To facilitate this, we strongly recommend using the 
[Conda](https://en.wikipedia.org/wiki/Conda_(package_manager)) 
environment management system, with a focus on the lightweight Conda 
package manager, [Miniconda](https://conda.io/miniconda.html).

> [!NOTE]
> This has currently been tested on _Fedora_- and 
> _Debian_-derived Unix systems

> [!NOTE]
> AnchorWave is currently not supported on Windows. 
> [See AnchorWave documentation for further details](https://github.com/baoxingsong/AnchorWave#installation).


## Get PHGv2
You can download the latest version of PHGv2 
[here](https://github.com/maize-genetics/phg_v2/releases/latest). 
Assuming you have downloaded PHGv2 locally, these instructions 
presume you will run the program directly. Obtain the `.tar` file 
manually from the provided link or use the following `curl` and `awk` 
commands to retrieve the latest release:

```shell
curl -s https://api.github.com/repos/maize-genetics/phg_v2/releases/latest \
| awk -F': ' '/browser_download_url/ && /\.tar/ {gsub(/"/, "", $(NF)); system("curl -LO " $(NF))}'
```

Once downloaded, untar the release using:
```
tar -xvf <PHGv2_release>.tar
```
...where `<PHGv2_release>.tar` is the
downloaded PHGv2 package. After the source has been decompressed,
we can remove the initial tar file using: 
```
rm <PHGv2_release>.tar
```


## "Installation"
No traditional installation is required, as the precompiled jar 
files are designed to function on any 
[POSIX](https://en.wikipedia.org/wiki/POSIX) platform meeting the 
specified requirements. Just open the downloaded package and place 
the folder containing the jar files and launch script in a preferred 
directory on your hard drive or server filesystem.

To run PHGv2, you can manually enter into the package directory and
run the `phg` wrapper script from the `bin` directory. Another
option is to add the wrapper script to your `PATH` variable. If you
are using the `bash` terminal shell, the classic syntax is:

```shell
export PATH="/path/to/phgv2-package/:$PATH"
```

...where `/path/to/phgv2-package/` is the path to the location of the
`phg` executable wrapper script.

> [!NOTE]
> The above path example must be the path to the `bin` subdirectory
> found in the `phg` directory.
 
> [!NOTE]
> The Java JAR files (`.jar`) in the `lib` subdirectory
> must remain in the same directory as `phg` for it to work.

> [!NOTE]
> Be sure to include the final `/` in your path.


## Test that PHGv2 works
To test that you can successfully run the `phg` executable. Run
the following command:

```shell
./phg --help
```

> [!NOTE]
> This assumes that you have added `phg` to your `PATH` using the
> above example, or you are within the `bin` subdirectory.

This should output summary text to the terminal including syntax
help and a list of subcommands and descriptions.


## Setting memory
The amount of data you wish to process will affect the amount of
computational resources that you will need. Since PHGv2 leverages
a [Java virtual machine](https://en.wikipedia.org/wiki/Java_virtual_machine) 
(JVM) for a majority of its tasks, we can manually alter the maximum 
amount of memory allocated to the JVM using the following command 
prompt:

```shell
export JAVA_OPTS="-Xmx<memory_amount>"
```

...where `<memory_amount>` is a specified unit of memory. For 
example, if I want to allocate a maximum of 50 gigabytes (GB) of 
memory for my operations, I would use the input `"-Xmx50g"`, where `g`
stands for GB:

```shell
export JAVA_OPTS="-Xmx50g"
```

> [!NOTE]
> In order for memory to properly be set, you must set this 
> **before** running any of the PHGv2 commands.

> [!NOTE]
> Setting JVM memory will only affect JVM-intensive commands. Since
> PGHv2 utilizes several external pieces of software several commands
> will **not** be affected by this. Currently, these are:
>
> * `setup-environment`
> * `initdb`
> * `align-assemblies`
> * `agc-compress`
> 
> ...which rely on conda, TileDB, AnchorWave, and AGC, respectively.

