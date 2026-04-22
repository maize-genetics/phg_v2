!!! warning "Legacy Documentation - PHG Version 1"

    This section contains documentation for **PHG version 1**, which is
    no longer actively developed. It is preserved here for archival and
    historical reference only. If you are looking to use the Practical
    Haplotype Graph, please refer to the [PHG v2 documentation](../../index.md),
    which reflects the current version of the software.

# Overview #
This shell script will extract out the haplotypes from the PHG DB and will then run minimap2 to index the pangenome.  This is needed for the next step [FindPathsMinimap2.sh](find_path_minimap2.md)  The outputs of this step is a fasta file containing the haplotypes for the graph and a .mmi index file created by minimap indexing said fasta file.

# This only needs to be run 1 time for a given graph.

## Pipeline Steps ##
1. The Script will first check to see if the fasta file already is located in the path specified.
    1. If not, Pull out the fasta using WriteFastaFromGraphPlugin
2. The Script will then check to see if the fasta file has already been indexed with the parameters specified.
    2. If not, Index using minimap2.

## Example Run Command ##

```
#!bash

IndexPangenome.sh [BASE_HAPLOTYPE_NAME] [CONFIG_FILE_NAME] [HAPLOTYPE_METHOD] [NUM_BASES_LOADED] [MINIMIZER_SIZE] [WINDOW_SIZE]

```

## Command Line Flags ##

```
#!bash
BASE_HAPLOTYPE_NAME: This is the base of the name of the haplotype fasta file.
CONFIG_FILE_NAME: This is the path to the config.txt file used to create the DB.  All the needed DB connection information should be in here.
HAPLOTYPE_METHOD: This is the HAPLOTYPE_METHOD_NAME used when creating the PHG.  This is used to pull out the correct haplotypes from the database.  This Name can also be a list of methods to be pulled at the same time.
NUM_BASES_LOADED: This is a Minimap2 parameter for how many Base Pairs should be loaded into the database for each batch.  Typically you should set this to be larger than the number of Base Pairs in your Pangenome.  If this is less, Minimap2 is inconsistent in its mapping. We have used 90G for the assembly DB in maize.
MINIMIZER_SIZE: This is the kmer size used by Minimap2 to do alignments.  We suggest you use a k of 21.
WINDOW_SIZE: This is the number of consecutive Minimizer windows used by Minimap.  We suggest using 11.
```


## Full Command Example

```
#!bash
DB=/workdir/user/DockerTuningTests/DockerOutput/phgTestMaizeDB.db
OUTPUT_DIR=/workdir/user/DockerTuningTests/DockerOutput/PangenomeFasta/
CONFIG_FILE=/workdir/user/DockerTuningTests/configSQLite.txt
docker run --name index_pangenome_container --rm \
                    -w / \
                    -v ${OUTPUT_DIR}:/tempFileDir/outputDir/pangenome/ \
                    -v ${DB}:/tempFileDir/outputDir/phgTestMaizeDB.db \
                    -v ${CONFIG_FILE}:/tempFileDir/data/configSQLite.txt \
                    -t maizegenetics/phg \
                    /IndexPangenome.sh pangenome_fasta configSQLite.txt CONSENSUS 4G 15 10
```

This command will extract the pan genome to a fasta named pangenome_fasta.fa and will also create an index file named pangenome_fasta.mmi for use in [FindPathMinimap2](find_path_minimap2.md). It will be saved to the directory specified by OUTPUT_DIR.

The --name parameter provides a name for the container.  This is optional.

The --rm parameter indicates the container should be deleted when the program finishes executing.  This is optional.

The -v directives are used to mount data from the user machine into the Docker.  The path preceding the ":" is the path on the user machine.  The directory path following the ":" are the paths inside the Docker where the user home directories will be mounted.  NOTE:  If you are using postgreSQL instead of SQLite db, you will not need to include the -V for the database file.

The "-t" directive indicates the Docker image of which this container will be an instance.  The last line tells the Docker container to run the IndexPangenome.sh script which is found in the root directory.  The items following are the parameters to the IndexPangenome.sh script.