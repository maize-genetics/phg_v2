!!! warning "Legacy Documentation - PHG Version 1"

    This section contains documentation for **PHG version 1**, which is
    no longer actively developed. It is preserved here for archival and
    historical reference only. If you are looking to use the Practical
    Haplotype Graph, please refer to the [PHG v2 documentation](../../index.md),
    which reflects the current version of the software.

The LoadGenomeIntervals.sh script reads reference genome sequence, genome intervals (anchor bed file) information and populates the initial database tables with reference genome data.  The database can be either postgres or SQLite.  The type of database and the means to access it are provided in a config file.

The db loading functionality is performed by a call to [LoadGenomeIntervalsToPHGdbPlugin](load_genome_intervals_to_db_plugin.md).  When LoadGenomeIntervals.sh is run as part of a Docker container script, the Docker script expects the following directory mount points:

* Mount localMachine:/pathToOutputDir/ to docker:/tempFileDir/outputDir/.  (If using sqlite, db mentioned in the config file should be in here.)
* Mount localMachine:/pathToReferenceDir/ to docker:/tempFileDir/data/reference (holds the reference fasta)
* Mount localMachine:/pathToDataDir/ to docker:/tempFileDir/data/ (config file and genome data file are here)
* Mount localMachine:/pathToAnswerDir/ to docker:/tempFileDir/answer/ (anchor bed file must live here)

The parameters to this script are:

* Config File:  DB Config File containing properties host,user,password,DB and DBtype where DBtype is either sqlite or postgres.  The shell script assumes the configuration file lives in the mounted /tempFileDir/data/ directory.  If using sqlite, the sqlite database should be written to the directory mounted to /tempFileDir/outputDir/. A sample config file can be found here: [Config File Wiki](https://bitbucket.org/bucklerlab/practicalhaplotypegraph/wiki/DockerPipeline/ConfigFile).  Only the first 5 entries in the config file deal with database connection.  The example config file contains additional parameters that are used in other steps of the pipeline.
* Reference Fasta File: The reference fasta file, which should exist in the folder mounted to the Docker /tempFileDir/data/reference directory.
* Anchor bed file:  A BED-formatted, tab-delimited file containing the chromosome, start-position, end-position for the genome intervals that will comprise the reference ranges.  This file is expected to live in the folder mounted to the Docker /tempFileDir/answer/ folder.
* Genome data file:  This is a tab-delimited file containing specific data related to the reference genome.  It must have the following columns:  Genotype Hapnumber Dataline Ploidy Reference GenePhased ChromPhased Confidence Method MethodDetails.  This file is expected to live in the folder mounted to the Docker /tempFileDir/data directory.  See this example genome data file: [Sample Genome Data File](../files/sample_load_data.txt)
* create_new boolean:  either "true" or "false", indicating whether a new database instance should be created.  If the value is "true", a new db as identified in the config file will be created to hold the reference anchor data.  if "false", the data created will be added to an existing database specified by the config file.  The boolean is "false" and there is no database of the specified name, an error is thrown and processing stops.  The "false" value allows users to add multiple reference versions to a single database.  We recommend against this, and recommend this value be "true" to facilitate creation of a new database.  WARNING:  "true" will delete any existing database of the specified name before recreating it.

The genome data file contents are described below:

* Genotype:  the name of the line as you want it to appear in the db genotypes table "name" column, e.g. "B104" or "B104_haplotype_caller".  These names must be unique for each line.  If your reference line is B73 and you also have a B73 haplotype you want to process from other data, then make the names distinct, e.g. B73Ref and B73.
* Hapnumber:  The 0-based chromosome number.  For inbreds, this is always 0.
* Dataline: Text description defining the line. This data will be stored in the "description" field of the genotypes table.
* Ploidy:  Number of chromosomes for this genome, the value to be stored in the "ploidy" field of the genome_intervals table.  Should be 1 for the reference.
* GenePhased: a boolean field indicating if the genes are phased.  Should be false for the reference.
* ChromPhased:  a boolean field indicating if the chromosomes are phased.  Should be false for the reference.
* Confidence:  a float field indicating confidence in phasing.  Generally 1 when no phasing.
* Method:  a mathod name by which this version of the reference data can be identified.
* MethodDetails:  Text description defining the method used to create the reference (or indicating the AGP version)

An example Docker script to run the LoadGenomesIntervals.sh script is:

```
#!python

# This script assumes the config.txt file lives in the directory mounted to /tempFileDir/data
# It assumes the reference fasta (in this case, Zea_mays.AGPv4.dna.toplevelMtPtv3.fa) lives
# in the directory mounted to /tempFileDir/data/reference.
# It assumes the genome intervals file (in the example below, maizeRefAnchor_intervals_bed) lives
# in the directory mounted to /tempFileDir/answer.
# It assumes the genomeData file describing the reference (in the example below, B&3Ref_load_data.txt) lives
# in the directory mounted to /tempFileDir/data
# It assumes your sqlite database lives in the directory mounted to /tempFileDir/outputDir/  This in only relevant when running an SQLite database.  This path shows up in the config file, parameter "db".

# You must change "/workdir/user/DockerTuningTests/..." to match your own directory paths
docker run --name load_phg_container --rm \
        -v /workdir/user/DockerTuningTests/DockerOutput/:/tempFileDir/outputDir/ \
        -v /workdir/user/DockerTuningTests/InputFiles/Reference/:/tempFileDir/data/reference/ \
        -v /workdir/user/DockerTuningTests/DataFolders/LoadRefDataDocker/:/tempFileDir/data/ \
        -v /workdir/user/DockerTuningTests/DataFolders/LoadRefDataDocker/:/tempFileDir/answer/ \
        -t maizegenetics/phg:latest \
        /LoadGenomeIntervals.sh config.txt Zea_mays.AGPv4.dna.toplevelMtPtv3.fa maizeRefAnchor_intervals.bed B73Ref_load_data.txt true
```

The --name parameter provides a name for the container.  This is optional.

The --rm parameter indicates the container should be deleted when the program finishes executing.  This is optional.

The -v directives are used to mount data from the user machine into the Docker.  The path preceding the ":" is the path on the user machine.  The directory path following the ":" are the paths inside the Docker where the user home directories will be mounted.

The "-t" directive indicates the Docker image of which this container will be an instance.  

The last line tells the Docker container to run the LoadGenomeIntervals.sh script which is found in the root directory.  The items following are the parameters to the LoadGenomeIntervals.sh script.