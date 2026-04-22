!!! warning "Legacy Documentation - PHG Version 1"

    This section contains documentation for **PHG version 1**, which is
    no longer actively developed. It is preserved here for archival and
    historical reference only. If you are looking to use the Practical
    Haplotype Graph, please refer to the [PHG v2 documentation](../../index.md),
    which reflects the current version of the software.

# Quick Start

The simplest way to load haplotypes from whole genome sequence (WGS : fastq, BAM, or GVCF) is to use the [PopulatePHGDBPipelinePlugin](populate_phgdb_pipeline.md). 
The parameter descriptions in the details section provide information about how the parameters are used and which ones you will need to set.
If your life is not complicated enough or you need more control over the process, you can run the scripts used by this plugin on your own as described in Details. 

# Details

The scripts used by the plugin and distributed as part of the PHG Docker are CreateHaplotypesFromFastq.groovy, 
CreateHaplotypesFromBAM.groovy, and CreateHaplotypesFromGVCF.groovy. These scripts can also be found in the PHG source code
repository.

To align sequence from fastq files to the reference genome and load the resulting haplotypes to the database, use the script
CreateHaplotypesFromFastq.groovy. This script aligns the reads to created BAM files then calls CreateHaplotypesFromBAM.groovy, which in turn calls CreateHaplotypesFromGVCF.groovy.
To start with BAM files run CreateHaplotypesFromBAM.groovy. To start with GVCF files, run CreateHaplotypesFromGVCF.groovy

Because the BAM files can be quite large, you may wish to delete them after the GVCFs are created. Saving the GVCFs is worthwhile
should it be necessary to rebuild the database.

## Sentieon

Sentieon is a commercial version of GATK that runs much faster. By default, GATK HaplotypeCaller is used to call
variants from BAM files. If a Sentieon license is provided in the config file, Sentieon will be used.

## Parameters

The following are parameters that can or be assigned in the config file. Some of the parameters are required and do not have default values. 
The default paths will work with the Docker PHG when the MakeDefaultDirectory plugin is used to create a directory structure.
Running CreateHaplotypesFromFastq.groovy uses all of the parameters. Running CreateHaplotypesFromBAM.groovy uses parameters for
that script plus CreateHaplotypesFromGVCF.groovy. Running CreateHaplotypesFromGVCF.groovy uses only CreateHaplotypesFromGVCF parameters.

### CreateHaplotypesFromFastq.groovy
| Name | Default | Description |
| ---- | ------- | ----------- |
| referenceFasta | required  | full path to the reference fasta file |
| wgsKeyFile | required | the key file listing the fastq's to be processed |
| numThreads | 10 | the maximum number of threads to be used |
| picardPath | /picard.jar | The path to picard.jar |
| fastqFileDir | /phg/inputDir/loadDB/fastq/ | the directory holding the fastq files |
| tempFileDir | /phg/inputDir/loadDB/temp/ | the storage location for temporary files created by the script |
| dedupedBamDir | /phg/inputDir/loadDB/bam/dedup/ | the stroace location for the deduped BAM files created by the script |

### CreateHaplotypesFromBAM.groovy

| Name | Default | Description |
| ---- | ------- | ----------- |
| referenceFasta | required | full path to the reference fasta file |
| numThreads | 10 | the maximum number of threads to be used |
| gvcfFileDir | /phg/inputDir/loadDB/gvcf/ | directory where gvcf files will be written |
| tempFileDir | /phg/inputDir/loadDB/bam/temp/ | directory where temporary files will be written |
| filteredBamDir | /phg/inputDir/loadDB/bam/mapqFiltered/ | directory where filtered BAMs will be written |
| dedupedBamDir | /phg/inputDir/loadDB/bam/dedup/ | directory where deduped BAMs will be written |
| Xmx | 10G | Max Java Heap Space used when running TASSEL code |
| extendedWindowSize | 1000 | adds flanking regions (base pairs) to reference ranges for BAM filtering |
| tasselLocation | /tassel-5-standalone/run_pipeline.pl | location of TASSEL, default works for docker |
| mapQ | 48 | minimum MapQ value used in the BAM file filtering |
| sentieon_license | optional | the sentieon license if you have one |
| sentieonPath | /sentieon/bin/sentieon | the path to the Sentieon executable |
| gatkPath | /gatk/gatk | the path to the gatk executable |

### CreateHaplotypesFromGVCF.groovy

| Name | Default | Description |
| ---- | ------- | ----------- |
| refRangeMethods | refRegionGroup |  |
| tempFileDir | /phg/inputDir/loadDB/bam/temp/ |  |
| tasselLocation | /tassel-5-standalone/run_pipeline.pl |  |
| extendedWindowSize | 1000 |  |
| Xmx | 10G |  |
