!!! warning "Legacy Documentation - PHG Version 1"

    This section contains documentation for **PHG version 1**, which is
    no longer actively developed. It is preserved here for archival and
    historical reference only. If you are looking to use the Practical
    Haplotype Graph, please refer to the [PHG v2 documentation](../../index.md),
    which reflects the current version of the software.

# Load genome intervals to PHG


## Details

This step populates the PHG database with a reference genome and creates reference intervals to be used in later steps of the pipeline. These tend to be conserved regions identified from your reference genome, but can be identified and defined any way you choose. Some options include defining gene regions from a gff file, exon regions, or arbitrary 1000 bps regions. Any region not specified in the reference intervals file will not be included in as reference range intervals in the database. Haplotype data is not added for any region that is not represented in the intervals file. After running this step you will have a database containing the reference genome and any reference range groups you have defined.

## Kitchen Sink

The LoadGenomeIntervals.sh script reads reference genome sequence, genome intervals (anchor bed file) information and populates the initial database tables with reference genome data.  The database can be either postgres or SQLite.  The type of database and the means to access it are provided in the config file.

There are 6 parameters used by this step:

* configFile: DB Config File containing properties host,user,password,DB and DBtype where DBtype is either sqlite or postgres. The shell script assumes the configuration file lives in the mounted /tempFileDir/data/ directory. If using sqlite, the sqlite database should be written to the directory mounted to /tempFileDir/outputDir/. The master config file can be found here: Config File. Only the first 5 entries in the config file deal with database connection. The example config file contains additional parameters that are used in other steps of the pipeline.
* referenceFastaFile: The reference genome fasta file. The path to this file should be included in the config file. 
* anchorBedFile: A [BED-formatted](create_phg_step1_bedfile.md), tab-delimited file containing the chromosome, interval start position, interval end position and "type" for the genome intervals that will comprise the reference ranges. The "type" column indicates the group to which the interval belongs.  The default groups, when creating this file from the CreateValidIntervalsFilePlugin,  are "FocusRegion" and "FocusComplement" which can be thought of as "genic" and "inter-genic" if the user were separating their regions by gene space.  The user may use any group names, but this field must be populated for each reference range in the anchors file.  The path to this file should be included in the config file.   NOTE: only ranges represented in the anchorBedFile will be included in the DB. This script will no longer automatically created inter-regions for those regions not in the anchor bed file. You may use the CreateValidIntervalsFilePlugin plugin to create a bedfile containing inter-regions.
* genomeDataFile: This is a tab-delimited file containing specific data related to the reference genome. It must have the following columns:  "Genotype" "Hapnumber" "Dataline" "Ploidy" "Reference" "GenesPhased" "ChromsPhased" "Confidence" "Method" "MethodDetails".  If running with PHG version 1.0 or greater, it additionally needs the column "gvcfServerPath".  The location of this file should be stored in the config file.  See this example genome data file: [Sample Genome Data File](../files/sample_load_data.txt).  Note this example genome data file contains the "gvcfServerPath" column needed when running PHG versions 1.0 or greater.
* refServerPath: This parameter provides information on where, in the real world, the genome fasta used to populate the reference haplotypes is stored.  This reflects a more permanent location, e.g. an irods directory.  Or could be a directory on AWS, or at some institution. It is just a string. 
* create_new boolean: either "true" or "false", indicating whether a new database instance should be created. If the value is "true", a new db as identified in the config file will be created to hold the reference anchor data. If "false", the data created will be added to an existing database specified by the config file. If the boolean is "false" and there is no database of the specified name, an error is thrown and processing stops.  The "false" value allows users to add multiple reference versions to a single database. We recommend against this, and recommend this value be "true" to facilitate creation of a new database. WARNING:  "true" will delete any existing database of the specified name before recreating it.

### *Details on running this step through docker*

When LoadGenomeIntervals.sh is run as part of a Docker container script, the Docker script expects the following directory mount points:

* Mount localMachine:/pathToOutputDir/ to docker:/phg/outputDir/. If using sqlite, the database mentioned in the config file should be in this directory.
* Mount localMachine:/pathToReferenceDir/ to docker:/phg/data/reference. This directory holds the reference genome fasta file
* Mount localMachine:/pathToDataDir/ to docker:/phg/data/. This contains the config file and genome data file
* Mount localMachine:/pathToAnswerDir/ to docker:/tempFileDir/answer/. This directory much contain the anchor bed file.

An example Docker script to run the LoadGenomesIntervals.sh script is:

```
#!python

# This script assumes the user is using the default directory structure created
# when the PHG MakeDefaultDirectoryPlugin TASSEL plugin has been run.

# It assumes the config.txt file lives in the top level working directory
# It assumes the reference fasta (in this case, Zea_mays.AGPv4.dna.toplevelMtPtv3.fa) lives
# in workingDir/inputDir/reference/
# It assumes the genome intervals file (in the example below, maizeRefAnchor_intervals_bed) lives
# in workingDir/inputDir/
# It assumes the genomeData file describing the reference genome lies in workingDir/inputDir/reference/ 
# It assumes your sqlite database lives in the workingDir/outputDir/  This in only relevant when running an SQLite database.  This path shows up in the config file, parameter "db".

# You must change "/workdir/user/DockerTuningTests/..." to match your own directory paths
docker run --name load_phg_container --rm \
        -v /workdir/user/MyGenomePHG/:/phg/ \
        -t maizegenetics/phg:latest \
        /LoadGenomeIntervals.sh config.txt Zea_mays.AGPv4.dna.toplevelMtPtv3.fa maizeRefAnchor_intervals.bed B73Ref_load_data.txt  "irods:/ibl/assemblies" true

```


Or run using the MakeInitialPHGDBPipelinePlugin  from the docker:


```
#!python
docker run --name small_seq_test_container --rm \
    -v /Users/lcj34/temp/phgSmallSeq/dockerBaseDir/:/phg/ \
    -t phgrepository_test:latest \ 
    /tassel-5-standalone/run_pipeline.pl -Xmx10G -debug -configParameters /phg/configSQLiteDocker.txt -MakeInitialPHGDBPipelinePlugin -endPlugin

```

The --name parameter provides a name for the container.  This is optional.

The --rm parameter indicates the container should be deleted when the program finishes executing.  This is optional.

The -v directives are used to mount data from the user machine into the Docker.  The path preceding the ":" is the path on the user machine.  The directory path following the ":" are the paths inside the Docker where the user home directories will be mounted.

The -t directive indicates the Docker image of which this container will be an instance.  

The last line tells the Docker container which script to run, either LoadGenomeIntervals.sh, or tassel with plugin MakeInitialPHGDBPipelinePlugin. Both tassel-t-standalone and LoadGenomeIntervals.sh script are found in the root directory.  The items following are the parameters to the LoadGenomeIntervals.sh script.

### *Files*

**Config file**

An example can be found here: Master config file

**Reference Fasta file**

Fasta file for the reference genome you want to use in the database

**Anchor bed file**

Bed file with selected reference range intervals. Information on how to create this file can be found here: [Create reference ranges](create_phg_step1_bedfile.md)

**Genome Data File**

An example of the genome data file can be found here: [Sample Genome Data File](../files/sample_load_data.txt)

The contents of the genome data file are described below:

* Genotype: The name of the line as you want it to appear in the database genotypes table "name" column (e.g. "B104" or "B104_haplotype_caller"). These names must be unique for each line. If your reference genome is B73 and you also have a B73 haplotype that you want to process from other data, then you must make the names distinct (e.g. B73Ref and B73).
* Hapnumber: The 0-based chromosome number. For inbreds, this is always 0. 
* Dataline: Text description defining the taxon. This data will be stored in the "description" field of the genotypes table.
* Ploidy: The number of chromosomes for this genome, the value to be stored in the "ploidy" field of the genome_intervals table. Should be 1 for the reference genome.
* GenesPhased: A boolean field indicating if the genes are phased. This should be false for the reference genome.
* ChromsPhased: A boolean field indicating if the chromosomes are phased. This should be false for the reference genome.
* Confidence: A float field indicating confidence in phasing. Leave this value as 1 when there is no phasing.
* Method: A method name by which this version of the reference data can be identified.
* MethodDetails: Text description defining the method used to create the reference (or indicating the AGP version)
* gvcfServerPath: The server and path where the gvcf files will be stored outside of the database.  This should be the remove server path, not the local path where these file are downloaded for processing, unless the remote and local server are the same.

### *Plugins*

#### **LoadAllIntervalsToPHGdbPlugin**

The database loading functionality is performed by a call to the LoadAllIntervalsToPHGdbPlugin in TASSEL, which is described below. 

This plugin takes a reference genome fasta, a genome data file and an anchors file as described above. The plugin code grabs sequence from the reference genome fasta based on the coordinates from the anchors file and loads the genome interval sequence to the specified PHG database. The genome data file contains genome specific details used when adding to the genotypes, gametes and method tables.
 
When finished loading the user defined regions, the intervals will be grouped based on the "name" field from the user provided intervals file. 

Database tables populated via this method are:

* genotypes
* gametes
* gamete_groups
* methods
* genome_interval_versions
* genome_intervals
* haplotypes
* gamete_haplotypes
* ref_range_ref_range_method

The parameters to this plugin are:

* -ref <Reference Genome File>  Reference Genome File for aligning against. (Default is null) (REQUIRED)
* -anchors <Anchors File> CSV file containing chrom, anchor start positions, anchor end position, gene start, gene end. The positions are physical positions (1-based, inclusive/inclusive).  (Default is null) (REQUIRED)
* -genomeData <Genome Data File> A tab-delimited file containing genome specific data with columns as described above. (Default is null) (REQUIRED)
* -refServerPath <External location of the reference genome> The external server/path where users may find the genome fasta used as the reference for this PHG db instance.

This plugin expects a database connection as input. This can be obtained by chaining the GetDBConnectionPlugin with the LoadGenomeIntervalsToPHGdbPlugin when called from the command line with the TASSEL run_pipeline.pl script.

An example of chaining these plugins is below:

```
#!python

/tassel-5-standalone/run_pipeline.pl -Xmx10G -debug \
-GetDBConnectionPlugin \
	-config dbConfigFile \
	-create true  \
	-endPlugin \
-LoadAllIntervalsToPHGdbPlugin \
	-ref reference \
	-anchors rangesFile \
	-genomeData genomeDataFile  \
        -refServerPath refServerString \
	-outputDir <pathToOutPut> \
	-endPlugin
```

## Troubleshooting

1. The reference genome name must be unique - if you have duplicate taxa, make sure they have unique labels when added to the database. 
2. Check that the genomeDataFile file is tab-delimited and that all fields are present. Check that there is not a trailing tab after the last field and that there are no extra columns in the file.



[Return to Step 1 pipeline version 0.0.40 or earlier](create_phg_step1_2_main.md)

[Return to Step 1 pipeline version 1.0 or later](create_phg_step1_create_db_load_ref.md)

[Return to Wiki Home](../home.md)
