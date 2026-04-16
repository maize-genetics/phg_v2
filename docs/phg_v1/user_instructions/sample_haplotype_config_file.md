!!! warning "Legacy Documentation - PHG Version 1"

    This section contains documentation for **PHG version 1**, which is
    no longer actively developed. It is preserved here for archival and
    historical reference only. If you are looking to use the Practical
    Haplotype Graph, please refer to the [PHG v2 documentation](../../index.md),
    which reflects the current version of the software.

# Sample Config file for Loading Haplotypes.


To Load haplotypes a config file will be needed.  Here is a sample one with most of the needed information filled out.  Anything marked with ***UNASSIGNED*** will need to be updated with a correct value and anything ***OPTIONAL*** is optional and depending on your use should be removed.


```
#!

###Example config file. 
### Anything marked with UNASSIGNED needs to be set for at least one of the steps
### If it is marked as OPTIONAL, it will only need to be set if you want to run specific steps. 
host=localHost
user=sqlite
password=sqlite
DB=/phg/smallSeqDB.db
DBtype=sqlite

#System parameters.  Xmx is the java heap size and numThreads will be used to set threads available for multithreading components.
Xmx=10G
numThreads=10

liquibaseOutdir=/phg/outputDir

anchors=***UNASSIGNED***
genomeData=***UNASSIGNED***

referenceFasta=***UNASSIGNED***

asmMethodName=mummer4
asmKeyFile=***OPTIONAL***

wgsMethodName=GATK_PIPELINE
wgsKeyFile=***OPTIONAL***

consensusMethodName=CONSENSUS
inputConsensusMethods=GATK_PIPELINE

fastqFileDir=/phg/inputDir/loadDB/fastq/
dedupedBamDir=/phg/inputDir/loadDB/bam/dedup/

gvcfFileDir=/phg/inputDir/loadDB/gvcf/
filteredBamDir=/phg/inputDir/loadDB/bam/mapqFiltered/

# BAM and GVCF uploading parameters
mapQ=48
refRangeMethods=FocusRegion,FocusComplement
extendedWindowSize=1000

# WGS Haplotype Filtering criteria.  These are the defaults.
GQ_min=50
QUAL_min=200
DP_poisson_min=.01
DP_poisson_max=.99
filterHets=true

## If you have a sentieon license you can set the server location here(and remove the #).  If it is set, it will use Sentieon instead of GATK
#sentieon_license= ***OPTIONAL***


##Consensus Plugin Parameters
minFreq=0.5
maxClusters=30
minSite=30
minCoverage=0.1
maxThreads=10
minTaxa=1
mxDiv=0.01

#This sets the type of clustering mode.
#Valid params are: upgma, upgma_assembly, and kmer_assembly
#The two assembly parameters are designed for assembly haplotypes and will choose a representative haplotype as the consensus instead of attempting to merge calls like with upgma.
clusteringMode=upgma

#If you want to use an assembly clusteringMode, you must have a ranking file.
#The ranking file must be a tab separated list of taxon\trankingScore where higher numbers are a better rank.  This file is used to chose the representative haplotype
rankingFile=***OPTIONAL***

##Optional if you want to use kmer_assembly as the clusteringMode. Otherwise is ignored 
kmerSize=7
distanceCalculation=Euclidean



```