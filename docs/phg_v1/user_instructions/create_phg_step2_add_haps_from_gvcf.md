!!! warning "Legacy Documentation - PHG Version 1"

    This section contains documentation for **PHG version 1**, which is
    no longer actively developed. It is preserved here for archival and
    historical reference only. If you are looking to use the Practical
    Haplotype Graph, please refer to the [PHG v2 documentation](../../index.md),
    which reflects the current version of the software.

# Filter GVCF and add variants to database

## Quick Start

1. Change `param1`, `param2`, `param3`, and `param4` to match file paths on your computer.
2. Run 
```
#!bash

CreateHaplotypesFromGVCF.groovy -[hd] -config [dbConfigFile]
```

**Command line flags**
```
#!bash

[-h, -help] 'Show usage information'
[-d, -description] 'Show information about this Pipeline'
[-config] 'DB Config File(required)'
```

## Details

This step is part of a pipeline that allows you to load Haplotypes into the Database. Currently we support inputs from WGS fastq files, BAM files of WGS reads aligned to a reference and a GVCF file. This step adds data from a GVCF file.

Steps to the full pipeline:

1. Check to see if the reference is BWA and fai indexed
  * If not, index using bwa, samtools, and picard
2. Align WGS fastq files to the reference using bwa mem
3. Sort and MarkDuplicates in the output BAM file
4. Filter the BAM file based on Minimum MapQ
5. Run HaplotypeCaller using GATK or Sentieon
6. Filter the GVCF file
7. Extract out the reference ranges and load the haplotypes to the database

This page describes the CreateHaplotypesFromBAM.groovy script, which runs step 7.

If you have data from multiple formats, we suggest running CreateHaplotypesFromFastq.groovy first on all fastq files, then CreateHaplotypesFromBAM.groovy on all BAM files, and finally CreateHaplotypesFromGVCF.groovy on all GVCF files.

## Kitchen Sink

This Groovy Script will create haplotypes originating from a set of filtered GVCF files specified in a [keyfile](create_phg_step2_haplotype_keyfile.md). This groovy script can be run by itself, but it is also called from [CreateHaplotypesFromBAM.groovy](create_phg_step2_add_haps_from_bam.md) ([source](https://bitbucket.org/bucklerlab/practicalhaplotypegraph/src/master/docker/buildFiles/scripts/groovy/CreateHaplotypesFromBAM.groovy)).  Basically if you already have GVCFs that are filtered to your liking and are ready to be uploaded, you should use this script.

A GVCF file contains more information than a regular VCF file. We use the [GATK GVCF format](https://gatk.broadinstitute.org/hc/en-us/articles/360035531812-GVCF-Genomic-Variant-Call-Format) for the PHG.

Note only records with type GVCF will be processed by this plugin.

If you need to filter, either use an external tool like [bcftools](https://samtools.github.io/bcftools/bcftools.html) or [vcftools](https://vcftools.github.io/index.html).  Alternatively you can use the [FilterGVCFSingleFilePlugin](create_phg_step2_filter_gvcf_single_file_plugin.md) in the PHG. The documentation for FilterGVCFPlugin describes the possible filters applied to the GVCF file.  

For Maize, we set mapQ=48, DP_poisson_min=0.01, DP_poisson_max=0.99, GQ_min=50 and filterHets=t. **NOTE THESE ARE FOR MAIZE WGS WITH GOOD COVERAGE >5X.**  These can be used as starting points, but you should not follow these values blindly.  

### **We strongly suggest you use an external tool so you know exactly what filters are applied to your GVCF file.**

This Groovy Script will open up the [keyFile](create_phg_step2_haplotype_keyfile.md)(location is specified as the LoadHaplotypesFromGVCFPlugin.keyFile entry in the config file) and will loop over all GVCF records and upload to the DB. If multiple GVCFs are stored for a given taxon, each GVCF will be treated as an independent haplotype.

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
* LoadHaplotypesFromGVCFPlugin.keyFile
    * Location of the [keyfile](create_phg_step2_haplotype_keyfile.md). **This is Required.** Note, if using Docker, this needs to be the Docker specific path to the keyfile.  You will need to mount this file to the Docker in order for it to work.
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

### *Details on running this step with wrapper scripts*

When running this step on the command line, all file paths and parameters are set in the config file. If you would like to overwrite the parameters set in the config file, you can do that by setting the parameters on the command line directly.

### *Details on running this step through docker*

An example Docker script to run the CreateHaplotypesFromFastq.groovy script is:

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

### *Files*

**Haplotype keyfile**

Information on the keyfile format is here: [Haplotype keyfile](create_phg_step2_haplotype_keyfile.md) 

**Reference fasta**

The reference genome in fasta file format

### *Steps*

#### **Upload GVCFs to Database**

Pull out the BED file from the database if this script is not called from CreateHaplotypesFromBAM.groovy.

Upload variants to the database for all intervals within the BAM file.

## Troubleshooting
1. If you are having problems with this step, first verify that your GVCF files are formatted correctly and that you have only a single sample per GVCF file.


[Return to Step 2 pipeline, version 0.0.40](create_phg_step1_2_main.md)

[Return to Step 2 pipeline, version 0.1.0+](create_phg_step2_assembly_and_wgs_haplotypes.md)

[Return to Wiki Home](../home.md)
