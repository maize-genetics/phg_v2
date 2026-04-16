!!! warning "Legacy Documentation - PHG Version 1"

    This section contains documentation for **PHG version 1**, which is
    no longer actively developed. It is preserved here for archival and
    historical reference only. If you are looking to use the Practical
    Haplotype Graph, please refer to the [PHG v2 documentation](../../index.md),
    which reflects the current version of the software.

# ConfigFile Documentation #

Most steps of the Practical Haplotype Graph require a Configuration File.  For the most part, users can create one config file and add parameters when needed, reusing it along the way.  


### Sample Configuration File ###

Here is a sample configuration file.  This file contains the required parameters and some additional optional ones.  Please refer to each individual pipeline's documentation page in this wiki to get a full list of available parameters and acceptable values for each.


```
#!bash

#database config parameters
host=localHost
user=sqlite
password=sqlite
DB=/tempFileDir/outputDir/DBFile.db
DBtype=sqlite

#Java arguments
Xmx=500G

#CreateHaplotype Params
referenceFasta=/path/to/reference/ref.fasta
LoadHaplotypesFromGVCFPlugin.keyFile=/path/to/keyfile/keyFile.txt
LoadHaplotypesFromGVCFPlugin.gvcfDir=/path/to/gvcfs/
LoadHaplotypesFromGVCFPlugin.referenceFasta=/path/to/reference/ref.fasta
LoadHaplotypesFromGVCFPlugin.haplotypeMethodName=GATK_PIPELINE
LoadHaplotypesFromGVCFPlugin.haplotypeMethodDescription=GVCF_DESCRIPTION

#Haplotype filtering
mapQ=48
DP_poisson_min=.01
DP_poisson_max=.99
GQ_min=50
filterHets=t

#sentieon license
sentieon_license=cbsulogin2.tc.cornell.edu:8990

#Consensus parameters
#Optional argument to get out merged VCF files for debugging consensus
includeVariants=true
mxDiv=.001
maxError=0.2

#FindPath Config parameters
BestHaplotypePathPlugin.maxNodes=30
BestHaplotypePathPlugin.minTaxa=5
BestHaplotypePathPlugin.minReads=1
BestHaplotypePathPlugin.maxReads=100
BestHaplotypePathPlugin.minTransitionProb=0.001
BestHaplotypePathPlugin.probCorrect=0.99
BestHaplotypePathPlugin.splitNodes=true
BestHaplotypePathPlugin.splitProb=0.99
```