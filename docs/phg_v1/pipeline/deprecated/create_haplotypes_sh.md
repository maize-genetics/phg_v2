!!! warning "Legacy Documentation - PHG Version 1"

    This section contains documentation for **PHG version 1**, which is
    no longer actively developed. It is preserved here for archival and
    historical reference only. If you are looking to use the Practical
    Haplotype Graph, please refer to the [PHG v2 documentation](../../../index.md),
    which reflects the current version of the software.

# THIS SCRIPT IS DEPRECATED

Use the Groovy Scripts instead.

The CreateHaplotypes.sh script will create haplotype sequences for each anchor and load them into the database.  The script can either start with raw WGS fastq files and a reference fasta or can operate on pre-aligned BAM files created with BWA mem.  Either fastq files or BAM files must be passed into the Docker container using volume mounting.  If neither are supplied, no haplotypes will be created.

# Mount Points For use with the PHG Docker#

* Mount localMachine:/pathToInputs/refFolder to docker://tempFileDir/data/reference/ (required)
* Mount localMaching:/pathToInputs/config.txt to docker://tempFileDir/data/config.txt (required)
* Mount localMachine:/pathToInputs/fastq to docker://tempFileDir/data/fastq/ (optional if WGS sequence needs to be aligned, if not declared, BAM files must be input)
* Mount localMachine:/pathToInputs/bam to docker://tempFileDir/data/bam/externalBams/ (optional to skip alignment, if not declared, Fastq files must be input.  BAMs must be in format TAXON_otherInformation.bam. We split on the "_" and extract the first token as the taxon/sample name)
* Mount localMachine:/pathToGVCFOutput/ to docker://tempFileDir/data/outputs/gvcfs/ (optional to get the output GVCF file created by HaplotypeCaller)
* Mount localMachine:/pathToGVCFFilteredOutput/ to docker://tempFileDir/data/filteredGvcfsAndFasta/Filtered/ (optional to get out the filtered gvcf, useful to see if your filtering parameters are too stringent)

# CreateHaplotypes.sh parameters #

* configFile: Path to config file containing DB parameters host, user, password, DB, type.  Used for making the database connection.  Type must be wither "sqlite" or "postgres" to identify db type for connection.  This file will also contain filtering criteria for both the BAM file filtering and the GVCF filtering.  A sample config file can be found here:[Config File Wiki](https://bitbucket.org/bucklerlab/practicalhaplotypegraph/wiki/DockerPipeline/ConfigFile)
* currentTaxon :  Name of taxon currently being processed.
* operationMode:  paired or single (for aligning via BWA)
* methodName: Name used to store this run of haplotypes into the database.  Later on it makes it easier to pull a specific set of haplotypes from the database when making the graph.
* fastqList: Optional list of fastq file names separated by commas to be aligned with BWA.  If paired, the length of this list must be even.  If this is blank, the script will assume you have BAM files to import.

# Pipeline Steps #
1. Check /tempFileDir/data/reference/ for bwa and .fai indices.  If they do not exist create new ones
2. Pull the Anchor Bed file from the Database
3. If fastq files are stored in /tempFileDir/data/fastq/ and they are specified as an argument, execute BWA MEM in either single or paired end mode(also specified by command line argument)
3a. Dedup the resulting BAM files
4. Apply a MapQ filter to remove any reads which map to multiple places on the genome(Default: MapQ>48)
5. Create GVCF files. If sentieon_license is specified in config.txt, the script will run Sentieon.  Else run GATK 4s implementation.
6. Apply a GVCF filter to the gvcf based on parameters in config.txt.  This is generally used to filter out abnormally high depth regions(repeats) or heterozygous regions. This step uses the [FilterGVCFPlugin](../filter_gvcf_plugin.md).
7. Load the final GVCF file to the DB.  This process will split the GVCF by anchor regions and will store both the VariantContext records(think VCF records) and the fasta sequence in the db.  This step uses the [LoadHapSequencesToDBPlugin](../load_hap_sequences_to_db_plugin.md).

# Relevant Config File Parameters. Must be in the form param=value#
A sample config file can be found here:[Config File Wiki](https://bitbucket.org/bucklerlab/practicalhaplotypegraph/wiki/DockerPipeline/ConfigFile)
### DB Parameters ###
* host - host name of the db
* user - db user name
* password - password for db
* DB - name of the db
* DBtype - sqlite or postgres

### Java Arguments ###
* Xmx - max heap space for the pipeline.  This is similar to Java's -Xmx argument.

### BAM Filtering Arguments ###
* mapQ - minimum mapQ filtering default it 48

### GVCF Filtering Arguments ###
#### These are optional and can be used in combination with one another ####
* DP_poisson_min - between 0 and 1. If used the pipeline will create a poisson distribution of the depths and filter out any depths below this poisson bound.
* DP_poisson_max - between 0 and 1. If used the pipeline will create a poisson distribution of the depths and filter out any depths above this poisson bound.
* DP_min - minimum Depth threshold.  Cannot be used in conjunction with DP_poisson_min or DP_poisson_max.
* DP_max - maximum Depth threshold.  Cannot be used in conjunction with DP_poisson_min or DP_poisson_max.
* GQ_min - minimum Genotype Quality(GQ annotation in VCF) threshold.  
* GQ_max - maximum Genotype Quality(GQ annotation in VCF) threshold.
* QUAL_min - minimum base pair quality(QUAL column in VCF) threshold. **WARNING** This parameter is heavily dependent on Depth.  We recommend to not use this option.
* QUAL_max - maximum base pair quality(QUAL column in VCF) threshold. **WARNING** This parameter is heavily dependent on Depth.  We recommend to not use this option.
* filterHets - true or false.  This will remove any VCF records which suggest a het based on Allele Depth(AD).  For simplicity it will remove any records where two or more alleles have greater than 0 depth.  In other words if the AD values look like this: 0,10,0 or 1,0,0 the record will stay.  If the values are 1,2,0 or 1,1,0 or 1,10,0 it will be filtered out.
* exclusionString - this will override all previous parameters with a bcftools valid exclusion string.  **This option is not recommended unless you have a very complicated or specific filter which needs to be applied to the gvcf**.

### Other Optional Parameters ###
* numThreads - number of threads which are available for HaplotypeCaller to run. (Default: 10).
* extendedWindowSize - When filtering BAM files, we will keep some flanking regions around the anchor in case HaplotypeCaller requires some additional information for the active window.  Default is 1000 bp on both sides of the anchor.
* sentieon_license - location of Sentieon Licensing server.  If this is not defined, the Pipeline will use GATK's implementation of HaplotypeCaller.  If it is, it will attempt to use Sentieon.


# Example Run Scripts #

An example shell script that runs a Docker container in a loop to process a list of taxon is here:
The --name parameter provides a name for the container.  This is optional.

The --rm parameter indicates the container should be deleted when the program finishes executing.  This is optional.

The -v directives are used to mount data from the user machine into the Docker.  The path preceding the ":" is the path on the user machine.  The directory path following the ":" are the paths inside the Docker where the user home directories will be mounted.

The "-t" directive indicates the Docker image of which this container will be an instance.  The last line tells the Docker container to run the CreateHaplotypes.sh script which is found in the root directory.  The items following are the parameters to the CreateHaplotypes.sh script.

If you are running on Cornell CBSU you need to change docker to docker1.

### This script is running with just BAMs ###
```
#!bash

#Set these properties
TAXON=B73
REF_DIR=/workdir/user/DockerTuningTests/InputFiles/Reference/
BAM_DIR=/workdir/user/DockerTuningTests/InputFiles/WGSBams/${TAXON}/
DB=/workdir/user/DockerTuningTests/DockerOutput/phgTestMaizeDB.db
CONFIG_FILE=/workdir/user/DockerTuningTests/DataFolders/LoadRefDataDocker/config.txt
#Optional debug parameters
GVCF_OUTPUT_DIR=/workdir/user/DockerTuningTests/DockerOutput/gvcfOut/${TAXON}/
GVCF_FILTERED_DIR=/workdir/user/DockerTuningTests/DockerOutput/gvcfOutFilter/${TAXON}/

docker run --name cbsu_phg_container_${TAXON} --rm \
        -v ${REF_DIR}:/tempFileDir/data/reference/ \
        -v ${BAM_DIR}:/tempFileDir/data/bam/${TAXON}/DedupBAMs/ \
        -v ${DB}:/tempFileDir/outputDir/phgTestMaizeDB.db \
        -v ${GVCF_OUTPUT_DIR}:/tempFileDir/data/outputs/gvcfs/ \
        -v ${GVCF_FILTERED_DIR}:/tempFileDir/data/filteredGvcfsAndFasta/Filtered/ \
        -v ${CONFIG_FILE}:/tempFileDir/data/config.txt \
        -t phgrepository_test:latest \
        /CreateHaplotypes.sh /tempFileDir/data/config.txt \
                          ${TAXON} \
                          single \
                          GATK_PIPELINE
```

### This script is running with fastq files in Single End Mode ###

```
#!bash

#Set these properties
TAXON=B73
REF_DIR=/workdir/user/DockerTuningTests/InputFiles/Reference/
FASTQ_DIR=/workdir/user/DockerTuningTests/InputFiles/Fastas/${TAXON}/
FASTQ_LIST=B73_batch1.fastq,B73_batch2.fastq,B73_batch3.fastq
DB=/workdir/user/DockerTuningTests/DockerOutput/phgTestMaizeDB.db
CONFIG_FILE=/workdir/user/DockerTuningTests/DataFolders/LoadRefDataDocker/config.txt
#Optional debug parameters
GVCF_OUTPUT_DIR=/workdir/user/DockerTuningTests/DockerOutput/gvcfOut/${TAXON}/
GVCF_FILTERED_DIR=/workdir/user/DockerTuningTests/DockerOutput/gvcfOutFilter/${TAXON}/

docker run --name cbsu_phg_container_${TAXON} --rm \
        -v ${REF_DIR}:/tempFileDir/data/reference/ \
        -v ${FASTQ_DIR}:/tempFileDir/data/fastq/ \
        -v ${DB}:/tempFileDir/outputDir/phgTestMaizeDB.db \
        -v ${GVCF_OUTPUT_DIR}:/tempFileDir/data/outputs/gvcfs/ \
        -v ${GVCF_FILTERED_DIR}:/tempFileDir/data/filteredGvcfsAndFasta/Filtered/ \
        -v ${CONFIG_FILE}:/tempFileDir/data/config.txt \
        -t phgrepository_test:latest \
        /CreateHaplotypes.sh /tempFileDir/data/config.txt \
                          ${TAXON} \
                          single \
                          GATK_PIPELINE \
                          $FASTQ_LIST
```

### This script is running with fastq files in Paired End Mode ###


```
#!bash

#Set these properties
TAXON=B73
REF_DIR=/workdir/user/DockerTuningTests/InputFiles/Reference/
FASTQ_DIR=/workdir/user/DockerTuningTests/InputFiles/Fastas/${TAXON}/
FASTQ_LIST=B73_batch1_1.fastq,B73_batch1_2.fastq,B73_batch2_1.fastq,B73_batch2_2.fastq,B73_batch3_1.fastq,B73_batch3_2.fastq
DB=/workdir/user/DockerTuningTests/DockerOutput/phgTestMaizeDB.db
CONFIG_FILE=/workdir/user/DockerTuningTests/DataFolders/LoadRefDataDocker/config.txt
#Optional debug parameters
GVCF_OUTPUT_DIR=/workdir/user/DockerTuningTests/DockerOutput/gvcfOut/${TAXON}/
GVCF_FILTERED_DIR=/workdir/user/DockerTuningTests/DockerOutput/gvcfOutFilter/${TAXON}/

docker run --name cbsu_phg_container_${TAXON} --rm \
        -v ${REF_DIR}:/tempFileDir/data/reference/ \
        -v ${FASTQ_DIR}:/tempFileDir/data/fastq/ \
        -v ${DB}:/tempFileDir/outputDir/phgTestMaizeDB.db \
        -v ${GVCF_OUTPUT_DIR}:/tempFileDir/data/outputs/gvcfs/ \
        -v ${GVCF_FILTERED_DIR}:/tempFileDir/data/filteredGvcfsAndFasta/Filtered/ \
        -v ${CONFIG_FILE}:/tempFileDir/data/config.txt \
        -t phgrepository_test:latest \
        /CreateHaplotypes.sh /tempFileDir/data/config.txt \
                          ${TAXON} \
                          paired \
                          GATK_PIPELINE \
                          $FASTQ_LIST
```

### This script is running running over a list of BAM files for multiple taxon ###


```
#!bash

taxonList=(B73 A632 B14 B37 B97 CO125 LH74 Ms71 Oh43 OH7B W22)

REF_DIR=/workdir/user/DockerTuningTests/InputFiles/Reference/
BAM_DIR=/workdir/user/DockerTuningTests/InputFiles/WGSBams/
DB=/workdir/user/DockerTuningTests/DockerOutput/phgTestMaizeDB.db
CONFIG_FILE=/workdir/user/DockerTuningTests/DataFolders/LoadRefDataDocker/config.txt
#Optional debug parameters
GVCF_OUTPUT_DIR=/workdir/user/DockerTuningTests/DockerOutput/gvcfOut/
GVCF_FILTERED_DIR=/workdir/user/DockerTuningTests/DockerOutput/gvcfOutFilter/

for TAXON in "${taxonList[@]}"
do
mkdir -p /workdir/user/DockerTuningTests/DockerOutput/gvcfOut/${TAXON}/
mkdir -p /workdir/user/DockerTuningTests/DockerOutput/gvcfOutFilter/${TAXON}/

                          
docker run --name cbsu_phg_container_${TAXON} --rm \
        -v ${REF_DIR}:/tempFileDir/data/reference/ \
        -v ${BAM_DIR}/${TAXON}/:/tempFileDir/data/bam/${TAXON}/DedupBAMs/ \
        -v ${DB}:/tempFileDir/outputDir/phgTestMaizeDB.db \
        -v ${GVCF_OUTPUT_DIR}/{TAXON}/:/tempFileDir/data/outputs/gvcfs/ \
        -v ${GVCF_FILTERED_DIR}/${TAXON}/:/tempFileDir/data/filteredGvcfsAndFasta/Filtered/ \
        -v ${CONFIG_FILE}:/tempFileDir/data/config.txt \
        -t phgrepository_test:latest \
        /CreateHaplotypes.sh /tempFileDir/data/config.txt \
                          ${TAXON} \
                          single \
                          GATK_PIPELINE

done
```