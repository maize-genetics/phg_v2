!!! warning "Legacy Documentation - PHG Version 1"

    This section contains documentation for **PHG version 1**, which is
    no longer actively developed. It is preserved here for archival and
    historical reference only. If you are looking to use the Practical
    Haplotype Graph, please refer to the [PHG v2 documentation](../../index.md),
    which reflects the current version of the software.

# Separate genome intervals into groups

## Quick Start

1. Change `param1`, `param2`, `param3`, and `param4` to match file paths on your computer.
2. Run 
```
#!bash

CreateHaplotypesFromBAM.groovy -[hd] -config [dbConfigFile]
```

**Command line flags**
```
#!bash

[-h, -help] 'Show usage information'
[-d, -description] 'Show information about this Pipeline'
[-config] 'DB Config File(required)'
```

## Details

This step is part of a pipeline that allows you to load Haplotypes into the Database. Currently we support inputs from WGS fastq files, BAM files of WGS reads aligned to a reference and a GVCF file. This step adds data from a BAM file.

Steps to the full pipeline:

1. Check to see if the reference is BWA and fai indexed
  * If not, index using bwa, samtools, and picard
2. Align WGS fastq files to the reference using bwa mem
3. Sort and MarkDuplicates in the output BAM file
4. Filter the BAM file based on Minimum MapQ
5. Run HaplotypeCaller using GATK or Sentieon
6. Filter the GVCF file
7. Extract out the reference ranges and load the haplotypes to the database

This page describes the CreateHaplotypesFromBAM.groovy script, which runs steps 4-7.

If you have data from multiple formats, we suggest running CreateHaplotypesFromFastq.groovy first on all fastq files, then CreateHaplotypesFromBAM.groovy on all BAM files, and finally CreateHaplotypesFromGVCF.groovy on all GVCF files.

## Kitchen Sink

This Groovy Script will create haplotypes originating from BAM Files specified within a [keyFile](create_phg_step2_haplotype_keyfile.md).  This groovy script can be run by itself, but it is also called from [CreateHaplotypesFromFastq.groovy](create_phg_step2_add_haps_from_fastq.md).  This expects BAM files where WGS short reads have been aligned to a reference using BWA MEM.

Note only records with type BAM and GVCF will be processed by this script.  BAM records will be run though HaplotypeCaller and filtered.  Then the GVCFs coming out of HaplotypeCaller will be added to the DB along with the GVCF records in the keyfile.

This Groovy Script will open up the [keyFile](create_phg_step2_haplotype_keyfile.md)(location is specified as the LoadHaplotypesFromGVCFPlugin.keyFile entry in the config file) and will loop over all FASTQ records and run BWA MEM to create BAM files. Then all BAM records in the keyfile and those output by BWA MEM are run through GATK/Sention's HaplotypeCaller. If multiple BAM files are stored for a given taxon, they will all be input for that taxon's HaplotypeCaller run. The resulting GVCFs will then be filtered and uploaded along with GVCF keyfile records to the DB. If multiple GVCFs are stored for a given taxon, each GVCF will be treated as an independent haplotype.

This CreateHaplotypes Script will need some configuration parameters set in order to work correctly. You must set referenceFasta, LoadHaplotypesFromGVCFPlugin.keyFile, LoadHaplotypesFromGVCFPlugin.gvcfDir, LoadHaplotypesFromGVCFPlugin.referenceFasta, LoadHaplotypesFromGVCFPlugin.haplotypeMethodName.

#### File Directories.
* gvcfFileDir=/tempFileDir/data/gvcfs/
    * The **Optional** Output GVCF file Directory.  If you need to keep these files, mount a local drive to this location
* tempFileDir=/tempFileDir/data/bam/temp/
    * The Temp file directory.  This Location is used to make any intermediate files when processing.
* filteredBamDir=/tempFileDir/data/bam/filteredBAMs/
    * The **Optional** Output filtered file Directory. After Filtering the BAM files by MapQ, the files are written here.  Mount a local drive to this location if you need the bams.
* dedupedBamDir=/tempFileDir/data/bam/DedupBAMs/
    * After filtering and Deduplication, the BAM files are stored here. This will be where the script is expecting the input BAM files.

#### TASSEL parameters
* Xmx=10G
    * Max Java Heap Space used when running TASSEL code.
* tasselLocation=/tassel-5-standalone/run_pipeline.pl
    * Location of TASSEL on machine.  If using PHG docker this is the correct location


#### PHG CreateHaplotypes Parameters
* referenceFasta
    * Reference fasta file location.  **This is Required.** Note, if using Docker, this needs to be the Docker specific path to the reference file.  You will need to mount these files to the Docker in order for it to work.
* LoadHaplotypesFromGVCFPlugin.keyFile
    * Location of the [keyfile](create_phg_step2_haplotype_keyfile.md). **This is Required.** Note, if using Docker, this needs to be the Docker specific path to the keyfile.  You will need to mount this file to the Docker in order for it to work.
* LoadHaplotypesFromGVCFPlugin.gvcfDir
    * Directory of the GVCF files you wish to upload. **This is Required.** Note, if using Docker, this needs to be the Docker specific path to the gvcfs.  You will need to mount these files to the Docker in order for it to work. **Note: This needs to match gvcfFileDir in the config file.**
* LoadHaplotypesFromGVCFPlugin.referenceFasta
    * This should be set to the same location as referenceFasta. **This is Required.**
* LoadHaplotypesFromGVCFPlugin.haplotypeMethodName=GATK_PIPELINE
    * This is the method name which you are going to write to the DB.  **This is Required.** If you attempt to upload the same gvcfs and do not change this an error will be thrown.
* LoadHaplotypesFromGVCFPlugin.haplotypeMethodDescription
    * This is the description for the haplotypeMethodName.  If this is not set and empty description will be written to the DB.  It is **strongly** suggested that you put something in this field.
* extendedWindowSize=1000
    * This Script will extract out the Focused reference ranges from the Database into a BED file.  When doing the BAM file filtering, we allow the user to add flanking regions around each of the regions.
* mapQ=48
    * Minimum MapQ value allowed in the BAM file filtering.  This is to reduce the number of reads which map to multiple locations across the genome.

#### GATK and Sentieon Parameters
* gatkPath=/gatk/gatk
    * Location of GATK on this system.  This is the location within the PHG docker
* numThreads=10
    * Number of threads requested to run HaplotypeCaller.  
* sentieon_license
    * If you have access to Sentieon and wish to use it, provide the sentieon license here.  The Script will automatically attempt to use Sentieon if this parameter is specified.
* sentieonPath=/sentieon/bin/sentieon
    * This is the Docker specific path to sentieon.  If sentieon is installed somewhere else on your system, change this parameter.

### *Details on running this step with wrapper scripts*

When running this step on the command line, all file paths and parameters are set in the config file. If you would like to overwrite the parameters set in the config file, you can do that by setting the parameters on the command line directly.

### *Details on running this step through docker*

An example Docker script to run the CreateHaplotypesFromFastq.groovy script is:

```
#!bash
REF_DIR=/workdir/user/DockerTuningTests/InputFiles/Reference/
DB=/workdir/user/DockerTuningTests/DockerOutput/phgTestMaizeDB.db
CONFIG_FILE=/workdir/user/DockerTuningTests/DataFolders/LoadRefDataDocker/config.txt
CONFIG_FILE_IN_DOCKER=/tempFileDir/data/config.txt
DEDUP_BAM_DIR=/workdir/user/DockerTuningTests/InputFiles/WGSBams/
FILTERED_BAM_DIR=/workdir/user/DockerTuningTests/InputFiles/WGSBams_Filtered/
GVCF_OUTPUT_DIR=/workdir/user/DockerTuningTests/DockerOutput/gvcfOut/
KEY_FILE=/workdir/user/DockerTuningTests/DataFolders/LoadRefDataDocker/keyfile.txt
KEY_FILE_IN_DOCKER=/tempFileDir/data/keyFile.txt

docker run --name small_seq_test_container --rm \
        -w / \
        -v ${REF_DIR}:/tempFileDir/data/reference/ \
        -v ${DB}:/tempFileDir/outputDir/phgSmallSeq.db \
        -v ${CONFIG_FILE}:${CONFIG_FILE_IN_DOCKER} \
        -v ${GVCF_OUTPUT_DIR}:/tempFileDir/data/gvcfs/ \
        -v ${DEDUP_BAM_DIR}:/tempFileDir/data/bam/DedupBAMs/ \
        -v ${FILTERED_BAM_DIR}:/tempFileDir/data/bam/filteredBAMs/ \
        -v ${KEY_FILE}:${KEY_FILE_IN_DOCKER} \
        -t maizegenetics/phg /CreateHaplotypesFromBAM.groovy -config ${CONFIG_FILE_IN_DOCKER}

```

Note that /tempFileDir/data/gvcfs/ and /tempFileDir/data/bam/filteredBAMs/ are locations for intermediate files to be stored.  These are not required, but if you wish to keep the GVCFs and/or the filtered BAM files you should mount these locations otherwise they will not persist.  If you would like to change the docker locations for these files, simply provide a gvcfFileDir and filteredBamDir config parameter in the config file and the script will write the files to that location.

Also note the -w / parameter. This is needed to guarantee that the script will run correctly. When running a normal docker, this is likely not needed, but if running on a system like cbsu, the working directory needs to be set to the root directory.

### *Files*

**Haplotype keyfile**

Information on the keyfile format is here: [Haplotype keyfile](create_phg_step2_haplotype_keyfile.md) 

**Reference fasta**

The reference genome in fasta file format

### *Steps*

#### **Filter BAMs**

Extract the BED file from the database for filtering and use with HaplotypeCaller. A BED file with extended anchors will be used to filter the BAM file.

Filter BAM files based on user-defined MapQ score and intervals contained in BED file. Re-index the filtered file.

#### **Run GATK HaplotypeCaller or Sentieon HaplotypeCaller**

Run HaplotypeCaller with Sentieon if Sentieon license is present in config file, otherwise with GATK. The BED file without extended anchors is used with this step.

#### **Run GVCF filtering**

Reheader the GVCF to make sure it matches the name of the taxon, then filter.

#### **Pass on to [CreateHaplotypesFromGVCF.groovy](create_phg_step2_add_haps_from_gvcf.md)**

## Troubleshooting

1. Make sure you have write permissions for all directories.



[Return to Step 2 pipeline, version 0.0.40-](create_phg_step1_2_main.md)

[Return to Step 2 pipeline, version 1.0+](create_phg_step2_assembly_and_wgs_haplotypes.md)

[Return to Wiki Home](../home.md)
