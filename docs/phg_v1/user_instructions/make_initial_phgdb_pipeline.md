!!! warning "Legacy Documentation - PHG Version 1"

    This section contains documentation for **PHG version 1**, which is
    no longer actively developed. It is preserved here for archival and
    historical reference only. If you are looking to use the Practical
    Haplotype Graph, please refer to the [PHG v2 documentation](../../index.md),
    which reflects the current version of the software.

# Make Initial PHG DB Pipeline Plugin

MakeInitialPHGDBPipelinePlugin is the first required step to build a PHG.  It will run the following plugins:

```
GetDBConnectionPlugin
LoadAllIntervalsToPHGdbPlugin
LiquibaseUpdatePlugin
```

GetDBConnectionPlugin will create a new DB with the required Schema and LoadAllIntervalsToPHGdbPlugin will then begin to populate the DB with Reference Range Information.  It will also load reference haplotypes into the PHG for the provided Reference Ranges.  LiquibaseUpdatePlugin is then run to mark the DB version as compatible with the current PHG software.   More information is in [Details](#markdown-header-Details). 

## Quick Start
Note this can be run outside of the docker as well.

You will need a valid intervals file for this step.  For details on the interval file contents and creation, click [here](create_phg_step1_bedfile.md) .

You will need to set WORKING_DIR for things to work correctly.  This should be the location where [MakeDefaultDirectory](make_default_directory.md) put all the PHG files.

```
WORKING_DIR=local/directory/where/MakeDefaultDirectory/was/run/
DOCKER_CONFIG_FILE=/phg/config.txt

docker run --name create_initial_db --rm \
    -v ${WORKING_DIR}/:/phg/ \
    -t maizegenetics/phg:latest \
    /tassel-5-standalone/run_pipeline.pl -Xmx100G -debug -configParameters ${DOCKER_CONFIG_FILE} \
    -MakeInitialPHGDBPipelinePlugin -endPlugin

```

Note we are assuming that you used MakeDefaultDirectoryPlugin and are using a config file named config.txt.  If you use something else, replace the file name in DOCKER_CONFIG_FILE.

## Details

As mentioned above, this plugin will first run GetDBConnectionPlugin which will create the DB schema for the PHG.  Then LoadAllIntervalsToPHGdbPlugin is run which will take in a BED file containing the coordinates for the Reference Ranges and will begin loading in that information.  Once the ReferenceRanges have been loaded, the plugin will also extract out sequence information for each Reference Range from the Reference Fasta file and will create a haplotype entry in the DB.  These haplotypes can be used for later analysis pipelines. 

Once the Reference Haplotypes are created, the LiquibaseUpdatePlugin is run to mark the DB as up-to-date, indicating that all existing database updates have been applied. This allows the DB to be updated via the Liquibase migration tool should future database changes occur. 
Further information can be found here:

[LoadAllIntervalsToPHGdbPlugin](create_phg_step1_load_genome_intervals.md)

[LiquibaseUpdatePlugin](update_phg_schema.md)

### Config File Parameters for this step
The MakeInitialPHGDBPipelinePlugin makes use of parameters stored in a config file. The relevant config file parameters are as follows. Change the values to match your configuration.  The DBType must be either "sqlite" or "postgres".  The host, user, password and DB should match your own system configuration.  Change other file names to match the files you have defined.  The value for "DB" should be a docker-relevant value if you are running this plugins via docker.  When running via docker, you will have mounted your local folder to something inside docker.  



```
host=localHost
user=sqlite
password=sqlite
DB=/phg/phg_db_name.db
DBtype=sqlite
# Load genome intervals parameters
referenceFasta=/phg/inputDir/reference/Ref.fa
anchors=/phg/anchors.bed
genomeData=/phg/inputDir/reference/load_genome_data.txt
refServerPath=irods:/ibl/home/assemblies/
#liquibase results output directory, general output directory
outputDir=/phg/outputDir
liquibaseOutdir=/phg/outputDir
```
This example config file follows the directory structure defined by [MakeDefaultDirectoryPlugin](make_default_directory.md) but is the Docker specific version of the paths.  If running outside of the docker, please make the necessary adjustments.


When running a postgresl data base, the value for "DB" should be just the db name without a path.  

Note the "refServerPath" is a field that stores where future users may find the source genome Fasta used for populating the reference genome. It can be any string the user wants.  This example shows a directory on an i iRods server.

[Return to Step 1 pipeline version 0.0.40 or earlier](create_phg_step1_2_main.md)

[Return to Step 1 pipeline version 0.1.0 or later](create_phg_step1_create_db_load_ref.md)