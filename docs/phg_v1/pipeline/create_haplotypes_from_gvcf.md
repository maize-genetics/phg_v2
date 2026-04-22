!!! warning "Legacy Documentation - PHG Version 1"

    This section contains documentation for **PHG version 1**, which is
    no longer actively developed. It is preserved here for archival and
    historical reference only. If you are looking to use the Practical
    Haplotype Graph, please refer to the [PHG v2 documentation](../../index.md),
    which reflects the current version of the software.

# Overview
This Groovy Script will create haplotypes originating from a set of filtered GVCF files specified in a [keyfile](https://bitbucket.org/bucklerlab/practicalhaplotypegraph/wiki/DockerPipeline/CreateHaplotypeKeyFile).  This groovy script can be run by itself, but it is also called from [CreateHaplotypesFromBAM.groovy](create_haplotypes_from_bam.md) ([source](https://bitbucket.org/bucklerlab/practicalhaplotypegraph/src/master/docker/buildFiles/scripts/groovy/CreateHaplotypesFromBAM.groovy)).  Basically if you already have GVCFs that are filtered to your liking and are ready to be uploaded, you should use this script.

### Note, as of PHG version 0.0.20, the keyFile plugin parameter has been renamed wgsKeyFile.

## More information about the keyfile format is here: 

https://bitbucket.org/bucklerlab/practicalhaplotypegraph/wiki/DockerPipeline/CreateHaplotypeKeyFile.

Note only records with type GVCF will be processed by this plugin.

If you need to filter, either use an external tool like [bcftools](https://samtools.github.io/bcftools/bcftools.html) or [vcftools](https://vcftools.github.io/index.html).  Alternatively you can use [FilterGVCFSingleFilePlugin](https://bitbucket.org/bucklerlab/practicalhaplotypegraph/wiki/DockerPipeline/FilterGVCFPlugin) in PHG. The Documentation for FilterGVCFPlugin describes the possible filters applied to the GVCF file.  

For Maize, we set mapQ=48, DP_poisson_min=0.01, DP_poisson_max=0.99, GQ_min=50 and filterHets=t. **NOTE THESE ARE FOR MAIZE WGS WITH GOOD COVERAGE >5X.**  These can be used as starting points, but you should not follow these values blindly.  

## **We strongly suggest you use an external tool so you know exactly what filters are applied to your GVCF file.**

## Pipeline Steps run by CreateHaplotypesFromGVCF
1. Extract out the Reference Ranges and load the haplotypes to the DB


## Example Run Command

```
#!bash

CreateHaplotypesFromGVCF.groovy -[hd] -config [dbConfigFile]
```

## Command Line Flags (note some parameters have both a short and long name)
```
#!bash

[-h, -help] 'Show usage information'
[-d, -description] 'Show information about this Pipeline'
[-config] 'DB Config File(required)'
```

This Groovy Script will open up the [keyFile](https://bitbucket.org/bucklerlab/practicalhaplotypegraph/wiki/DockerPipeline/CreateHaplotypeKeyFile)(location is specified as the LoadHaplotypesFromGVCFPlugin.keyFile entry in the config file) and will loop over all GVCF records and upload to the DB.  If multiple GVCFs are stored for a given taxon, each GVCF will be treated as an independent haplotype.

## Sample Script:


```
#!bash
#Set these properties
REF_DIR=/workdir/user/DockerTuningTests/InputFiles/Reference/
GVCF_DIR=/workdir/user/DockerTuningTests/InputFiles/
GVCF_DIR_IN_DOCKER=/tempFileDir/data/outputs/gvcfs/
DB=/workdir/user/DockerTuningTests/DockerOutput/phgTestMaizeDB.db
CONFIG_FILE=/workdir/user/DockerTuningTests/DataFolders/LoadRefDataDocker/config.txt
CONFIG_FILE_IN_DOCKER=/tempFileDir/data/config.txt
KEY_FILE=/workdir/user/DockerTuningTests/DataFolders/LoadRefDataDocker/keyfile.txt
KEY_FILE_IN_DOCKER=/tempFileDir/data/keyFile.txt

docker run --name small_seq_test_container --rm \
        -v ${REF_DIR}:/tempFileDir/data/reference/ \
        -v ${GVCF_DIR}:${GVCF_DIR_IN_DOCKER} \
        -v ${DB}:/tempFileDir/outputDir/phgTestMaizeDB.db \
        -v ${CONFIG_FILE}:${CONFIG_FILE_IN_DOCKER} \
        -v ${KEY_FILE}:${KEY_FILE_IN_DOCKER} \
        -t maizegenetics/phg \
        /CreateHaplotypesFromGVCF.groovy \
            -config ${CONFIG_FILE_IN_DOCKER}

```

Note the GVCF_DIR_IN_DOCKER needs to also be specified in the config file under the Parameter:LoadHaplotypesFromGVCFPlugin.gvcfDir

## Config File Parameters
This CreateHaplotypes Script will need some configuration parameters set in order to work correctly.  These must go in you config.txt file. You must set referenceFasta, LoadHaplotypesFromGVCFPlugin.keyFile, LoadHaplotypesFromGVCFPlugin.gvcfDir, LoadHaplotypesFromGVCFPlugin.referenceFasta, LoadHaplotypesFromGVCFPlugin.haplotypeMethodName.

#### File Directories.
* tempFileDir=/tempFileDir/data/bam/temp/
    * The Temp file directory.  This Location is used to make any intermediate files when processing.  For this script, only the BED File will be extracted containing the Focused Reference Ranges setup in a previous step.

#### TASSEL parameters
* Xmx=10G
    * Max Java Heap Space used when running TASSEL code.
* tasselLocation=/tassel-5-standalone/run_pipeline.pl
    * Location of TASSEL on machine.  If using PHG docker this is the correct location


#### PHG CreateHaplotypes Parameters
* referenceFasta
    * Reference fasta file location.  **This is Required.** Note, if using Docker, this needs to be the Docker specific path to the reference file.  You will need to mount these files to the Docker in order for it to work.
* LoadHaplotypesFromGVCFPlugin.wgsKeyFile
    * Location of the [keyfile](https://bitbucket.org/bucklerlab/practicalhaplotypegraph/wiki/DockerPipeline/CreateHaplotypeKeyFile). **This is Required.** Note, if using Docker, this needs to be the Docker specific path to the keyfile.  You will need to mount this file to the Docker in order for it to work.  Note before PHG version 0.0.20, this parameter was named keyFile.
* LoadHaplotypesFromGVCFPlugin.gvcfDir
    * Directory of the GVCF files you wish to upload. **This is Required.** Note, if using Docker, this needs to be the Docker specific path to the gvcfs.  You will need to mount these files to the Docker in order for it to work. 
* LoadHaplotypesFromGVCFPlugin.referenceFasta
    * This should be set to the same location as referenceFasta. **This is Required.**
* LoadHaplotypesFromGVCFPlugin.haplotypeMethodName=GATK_PIPELINE
    * This is the method name which you are going to write to the DB.  **This is Required.** If you attempt to upload the same gvcfs and do not change this an error will be thrown.
* LoadHaplotypesFromGVCFPlugin.haplotypeMethodDescription
    * This is the description for the haplotypeMethodName.  **This is Required.**
* extendedWindowSize=1000
    * This Script will extract out the Focused reference ranges from the Database into a BED file.  When doing the BAM file filtering, we allow the user to add flanking regions around each of the regions.
* numThreads=10
    * Number of threads requested to run HaplotypeCaller.  
* refRangeMethods=refRegionGroup
    * This is used to extract a BED file out of the DB before the GVCF file is processed.  The BED file is then used to extract out regions of the GVCF used to become the haplotypes.  Typically, refRegionGroup refers to the anchor Reference ranges.  If "refRegionGroup,refInterRegionGroup" is used it will create a BED file representing both anchors and inter anchors.  **We strongly suggest not setting this parameter in the Config File**