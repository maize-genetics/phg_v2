!!! warning "Legacy Documentation - PHG Version 1"

    This section contains documentation for **PHG version 1**, which is
    no longer actively developed. It is preserved here for archival and
    historical reference only. If you are looking to use the Practical
    Haplotype Graph, please refer to the [PHG v2 documentation](../../index.md),
    which reflects the current version of the software.

# Extract a pangenome fasta file

## Quick Start

1. Change `param1` to match file paths on your computer.
2. Run `phg extractPangenome [config.txt]`

## Details

This step extracts haplotypes from the PHG database and writes these haplotypes to a pangenome fasta file which is then indexed. Both the fasta file and the .mmi minimap index file are output. The fasta file will be formatted so that the header of each sequence is named with the haplotype ID from the database. 

**This step only needs to be executed once**

## Kitchen Sink

This shell script will extract out the haplotypes from the PHG database and will then run minimap2 to index the pangenome. This is needed for subsequent steps in the [findPaths](impute_with_phg_main.md) pipeline. The outputs of this step are a fasta file containing the haplotypes for the graph and a .mmi index file created by minimap indexing said fasta file.

he fasta file will be formatted so that the header of each sequence is named with the haplotype ID from the database. Every haplotype from every reference range will be included in this pangenome fasta file, so the total number of sequences will be equal to the number of haplotypes * the number of reference ranges in your database.

## Pipeline Steps ##
1. The Script will first check to see if the fasta file already is located in the path specified.
    1. If not, Pull out the fasta using WriteFastaFromGraphPlugin
2. The Script will then check to see if the fasta file has already been indexed with the parameters specified.
    2. If not, Index using minimap2.

There are 5 parameters that can be set in this step.

* BASE_HAPLOTYPE_NAME: This is the base of the name of the haplotype fasta file.
* CONFIG_FILE_NAME: This is the path to the config.txt file used to create the DB. All the needed DB connection information should be in here.
* HAPLOTYPE_METHOD: This is the HAPLOTYPE_METHOD_NAME used when creating the PHG. This is used to pull out the correct haplotypes from the database. This name can also be a list of methods to be pulled at the same time.
* NUM_BASES_LOADED: This is a Minimap2 parameter for how many base pairs should be loaded into the database for each batch. Typically you should set this to be larger than the number of Base Pairs in your Pangenome. If this is less, Minimap2 is inconsistent in its mapping. We have used 90G for the assembly DB in maize.
* MINIMIZER_SIZE: This is the kmer size used by Minimap2 to do alignments. We suggest you use a k of 21.
* WINDOW_SIZE: This is the number of consecutive Minimizer windows used by Minimap.  We suggest using 11.

### *Details on running this step with wrapper scripts*

This step can be run directly from the command line or from a wrapper script. When running this step with a wrapper script, all file paths and parameters are set in the config file and the only call that needs to be run in the terminal is `phg extractPangenome /path/to/config.txt`. 

If you would like to overwrite the parameters set in the config file, you can do that by setting the parameters on the command line directly.

For example, to ignore the config file base haplotype name and set one directly, you could run:
```
phg extractPangenome -CONFIG_FILE_NAME /path/to/config.txt -BASE_HAPLOTYPE_NAME Pangenome1
```

You can also run this script directly from the command line using the IndexPangenome.sh bash script. An example command is below:
```
#!bash

IndexPangenome.sh [BASE_HAPLOTYPE_NAME] [CONFIG_FILE_NAME] [HAPLOTYPE_METHOD] [NUM_BASES_LOADED] [MINIMIZER_SIZE] [WINDOW_SIZE]

```

### *Details on running this step through docker*

Mount points for use with the PHG Docker:

* Mount ${OUTPUT_DIR} (localMachine:/pathToOutputDir/) to docker:/tempFileDir/outputDir/. If using sqlite, the database mentioned in the config file should be in this directory.
* Mount ${DB} (localMachine:/pathToDBDir/DB) to docker:/tempFileDir/outputDir/phgTestMaizeDB.db. This should be your PHG database file. 
* Mount ${CONFIG_FILE} (localMachine:/pathToDataDir/) to docker:/tempFileDir/data/configSQLite.txt. This contains the config file.

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

This command will extract the pan genome to a fasta named pangenome_fasta.fa and will also create an index file named pangenome_fasta.mmi for use in subesquent pathfinding steps. It will be saved to the directory specified by OUTPUT_DIR.

The --name parameter provides a name for the container.  This is optional.

The --rm parameter indicates the container should be deleted when the program finishes executing.  This is optional.

The -v directives are used to mount data from the user machine into the Docker.  The path preceding the ":" is the path on the user machine.  The directory path following the ":" are the paths inside the Docker where the user home directories will be mounted.  NOTE:  If you are using postgreSQL instead of SQLite db, you will not need to include the -V for the database file.

The -t directive indicates the Docker image of which this container will be an instance.  The last line tells the Docker container to run the IndexPangenome.sh script which is found in the root directory.  The items following are the parameters to the IndexPangenome.sh script.

Also note the -w / parameter. This is needed to guarantee that the script will run correctly. When running a normal docker, this is likely not needed, but if running on a system like cbsu, the working directory needs to be set to the root directory.

### *Files*

**Config file**

An example can be found here: Master config file

### *Steps and Plugins*

#### **WriteFastaFromGraphPlugin**

The WriteFastaFromGraphPlugin method pulls sequence information from the PHG database and writes the sequence to a fasta file, with each entry in the fasta file corresponding to one haplotype from one reference range. The fasta entries are named with the haplotype_id from the database. The plugin will also BWA index the fasta file if possible.

The parameters to this plugin are:

* -outputFile <Output file> Output filename prefix for pangenome fasta. (REQUIRED)

The plugin requires a haplotype graph, which can be acquired by chaining the HaplotypeGraphBuilderPlugin with the WriteFastaFromGraphPlugin. Here is an example command to call this pipeline with TASSEL:
```
perl /tassel-5-standalone/run_pipeline.pl $fullXmx -debug \
-HaplotypeGraphBuilderPlugin \
	-configFile $CONFIGFILE \
	-methods $HAPLOTYPE_METHOD \
	-includeVariantContexts false \
	-endPlugin \
-WriteFastaFromGraphPlugin \
	-outputFile $HAPLOTYPE_FASTA \
	-endPlugin
```

#### **Index fasta with Minimap2**

If the haplotype fasta index file does not exist, you can create one with minimap2. This is needed for the next step in the pipeline.
```
time /minimap2/minimap2 -d ${HAPLOTYPE_INDEX} -k ${MINIMIZER_SIZE} -I ${NUM_BASES_LOADED} -w ${WINDOW_SIZE} ${HAPLOTYPE_FASTA}
```

## Troubleshooting

1. If you seem to be getting poor-quality mappings, try increasing the value for the NUM_BASES_LOADED parameter. This is a Minimap2 parameter for how many base pairs should be loaded into the database for each batch. Typically you should set this to be larger than the number of base pairs in your pangenome. If this is less than the number of base pairs in your pangenome, Minimap2 is inconsistent in its mapping. We have used 90G for the assembly DB in maize, but it may need to be set larger for larger genomes.


[Return to Step 3 pipeline](impute_with_phg_main.md)

[Return to Wiki Home](../home.md)