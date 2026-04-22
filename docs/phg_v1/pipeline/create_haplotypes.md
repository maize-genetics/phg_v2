!!! warning "Legacy Documentation - PHG Version 1"

    This section contains documentation for **PHG version 1**, which is
    no longer actively developed. It is preserved here for archival and
    historical reference only. If you are looking to use the Practical
    Haplotype Graph, please refer to the [PHG v2 documentation](../../index.md),
    which reflects the current version of the software.

# Create Haplotypes

## Overview
This step allows you to load Haplotypes into the Database.  Currently we support inputs from WGS fastq files, BAM files of WGS reads aligned to a reference and a GVCF file.
Depending on the inputs, you will need to run a different Script.

### Note, as of PHG version 0.0.20, the keyFile plugin parameter has been renamed wgsKeyFile.

### Scripts

* [CreateHaplotypesFromFastq.groovy](https://bitbucket.org/bucklerlab/practicalhaplotypegraph/wiki/DockerPipeline/CreateHaplotypesFromFastq)
* [CreateHaplotypesFromBAM.groovy](https://bitbucket.org/bucklerlab/practicalhaplotypegraph/wiki/DockerPipeline/CreateHaplotypesFromBAM)
* [CreateHaplotypesFromGVCF.groovy](https://bitbucket.org/bucklerlab/practicalhaplotypegraph/wiki/DockerPipeline/CreateHaplotypesFromGVCF)

Please refer to each scripts documentation for more detailed instructions.

### Steps of the Pipeline
1. Check to see if the reference is BWA and fai indexed
    1. If not index using bwa, samtools and picard
2. Align WGS fastqs to the reference using bwa mem
3. Sort and MarkDuplicates in the output BAM file
4. Filter the BAM file based on Minimum MapQ
5. Run HaplotypeCaller using GATK or Sentieon
6. Filter the GVCF file
7. Extract out the Reference Ranges and load the haplotypes to the DB

Steps 1-7 are completed by CreateHaplotypesFromFastq.groovy, Steps 4-7 are completed by CreateHaplotypesFromBAM.groovy and only step 7 is completed by CreateHaplotypesFromGVCF.groovy.

If you have data from multiple different formats, we suggest running CreateHaplotypesFromFastq.groovy first on all your fastq files, then CreateHaplotypesFromBAM.groovy on your BAMs and finally CreateHaplotypesFromGVCF.groovy on your GVCF files.


Deprecated Script:

* [CreateHaplotypes.sh](https://bitbucket.org/bucklerlab/practicalhaplotypegraph/wiki/DockerPipeline/DeprecatedScripts/CreateHaplotypes_sh) - This script is no longer maintained.