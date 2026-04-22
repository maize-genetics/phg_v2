!!! warning "Legacy Documentation - PHG Version 1"

    This section contains documentation for **PHG version 1**, which is
    no longer actively developed. It is preserved here for archival and
    historical reference only. If you are looking to use the Practical
    Haplotype Graph, please refer to the [PHG v2 documentation](../../index.md),
    which reflects the current version of the software.

Steps to run Minimap2BasedPipeline without needing a docker.

#Prerequisites:
* A Database filled with haplotypes
* The most recent tassel-5-standalone from bitbucket
* The most recent phg.jar built from practicalhaplotypegraph project.  Put this in tassel-5-standalone/lib/
* The most recent version of minimap2 downloaded and installed.

#Step 1: pull out the fasta files

```
#!bash

time ./tassel-5-standalone/run_pipeline.pl -Xmx100g -debug -HaplotypeGraphBuilderPlugin -configFile [CONFIG_FILE_NAME] -methods [HAPLOTYPE_METHOD_NAME] -includeVariantContexts false -endPlugin -WriteFastaFromGraphPlugin -outputFile [OUTPUT_FASTA_NAME] -endPlugin > [LOG_FILE]
```



##Parameters:
* CONFIG_FILE_NAME - The name of the config file you are wanting to run through the pipeline.  This is used to figure out what DB to access.
* HAPLOTYPE_METHOD_NAME - The name of the Haplotype Methods used when populating the DB.  
* OUTPUT_FASTA_NAME - The name of the output fasta file you are creating by running this step.
* LOG_FILE - file to write the logs to.

#Step 2: Index the pangenome fasta using minimap2 (~30 minutes for 20 assemblies with 75k anchors)

```
#!bash
time minimap2/minimap2 -d [OUTPUT_MINIMAP_INDEX] -k 21 -w 11 -I 90G [PANGENOME_FASTA_FILE] 

```

##Parameters:
* OUTPUT_MINIMAP_INDEX - output file name for the minimap2 index file.  This should have a .mmi extension
* -k minimap2 parameter for kmer-size.  We Tested 21 to be the best for maize.
* -w minimap2 parameter for window size.  We tested 11 to be the best for maize.
* -I - minimap2 parameter may need to be dropped if running into memory issues.  This should be above the number of basepairs in the entire pangenome.  If it is too low it will batch the index and cause the mapping to be incorrect.
* PANGENOME_FASTA_FILE - the output from the last step. 

#Step 3: Run FastqDirToMappingPlugin(This one is the slowest step.  Highly dependent on read coverage.)

```
#!bash

time ./tassel-5-standalone/run_pipeline.pl -debug -Xmx200G -HaplotypeGraphBuilderPlugin -configFile [CONFIG_FILE] -methods [HAPLOTYPE_METHOD_NAME] -includeVariantContexts false -includeSequences false  -endPlugin -FastqDirToMappingPlugin -minimap2IndexFile [MINIMAP_INDEX_FILE] -fastqDir [READ_DIRECTORY]/ -mappingFileDir [OUTPUT_MAPPING_DIR]/ -paired [PAIRED] -endPlugin > [LOG_FILE]
```

##Parameters:
* CONFIG_FILE_NAME - The name of the config file you are wanting to run through the pipeline.  This is used to figure out what DB to access.
* HAPLOTYPE_METHOD_NAME - The name of the Haplotype Methods used when populating the DB.  
* MINIMAP_INDEX_FILE - Minimap2 index file made in previous step
* READ_DIRECTORY - directory holding all the reads you wish to align.  Make sure to have a "/" at the end.
* OUTPUT_MAPPING_DIR - directory where all the mapping files are stored. Make sure to have a "/" at the end. This will go away likely by fall 2019.
* PAIRED - either "true" or "false". If true, it will attempt to pair up reads based on file name and use the pair to align.  For this to work correctly the file names need to be in the form B73_1.fq and B73_2.fq or even B73_R1.fq and B73_R2.fq.  Basically everything but the last token separated by "_" must be the same for the pair.  If false it will map each file independently and make a mapping file for each.  False is generally used for GBS.

#Step 4: Find Path

```
#!bash

time ./tassel-5-standalone/run_pipeline.pl -debug -Xmx240g -HaplotypeGraphBuilderPlugin -configFile [CONFIG_FILE_NAME] -methods [HAPLOTYPE_METHOD_NAME],[REF_RANGE_METHOD_NAME] -includeVariantContexts false -includeSequences false -endPlugin -HapCountBestPathToTextPlugin -configFile [CONFIG_FILE_NAME] -inclusionFileDir [MAPPING_DIR]/ -outputDir [OUTPUT_PATH_DIR]/ -hapCountMethod [HAP_COUNT_METHOD_NAME] -pMethod [PATH_METHOD_NAME] -endPlugin > [LOG_FILE] 2>&
```

##Parameters:
* CONFIG_FILE_NAME - The name of the config file you are wanting to run through the pipeline.  This is used to figure out what DB to access.
* HAPLOTYPE_METHOD_NAME - The name of the Haplotype Methods used when populating the DB.  
* REF_RANGE_METHOD_NAME - Name of the reference ranges you wish to run on.  It is highly suggested you use "refRegionGroup" as this is the genic regions.  This approach does not work very well with the intergenes.
* MAPPING_DIR - Name of the mapping directory output from the previous step.
OUTPUT_PATH_DIR - Output directory for the path output files.
* HAP_COUNT_METHOD_NAME - name of the hap counts stored in the DB.  Just put a dummy value here for now.  It will be used in the future.
* PATH_METHOD_NAME - name of the path method to be stored in the DB.  Just put a dummy value here for now. It will be used in the future.

Depending on what you are attempting to do, you may need to tweak the config file parameters.  In particular:

* minReads - this will control the amount of imputation the HMM will do.  If it is set to 0 the HMM will attempt to find the path through all reference ranges even if we did not have any reads supporting those haplotypes.  If set to 1, we require at least one Read hitting haplotypes in the reference range.
* splitTaxa - should probably set this to true when using assemblies.  Assemblies tend to skip over some reference regions which can bottle neck the HMM.  By doing splitTaxa=true, this is solved.

Here are the parameters we typically use:

```
#!
maxNodesPerRange=30
minTaxaPerRange=1
minReads=0
maxReadsPerKB=1000
minTransitionProb=0.001
probReadMappedCorrectly=0.99
emissionMethod=allCounts
splitTaxa=true
useBF=false

```

#Step 5: Create VCF(Optional)

```
#!bash
time ./tassel-5-standalone/run_pipeline.pl -debug -Xmx100g -HaplotypeGraphBuilderPlugin -configFile [CONFIG_FILE_NAME] -methods [HAPLOTYPE_METHOD_NAME],[REF_RANGE_METHOD_NAME] -includeVariantContexts true -endPlugin -ImportHaplotypePathFilePlugin -inputFileDirectory [PATH_DIR_NAME]/ -endPlugin -PathsToVCFPlugin -positionVCF [KNOWN_SITES_VCF] -ref [REFERENCE_FILE] -outputFile [OUTPUT_VCF_FILE] -endPlugin > [LOG_FILE] 2>&1

```

##Parameters:
* CONFIG_FILE_NAME - The name of the config file you are wanting to run through the pipeline.  This is used to figure out what DB to access.
* HAPLOTYPE_METHOD_NAME - The name of the Haplotype Methods used when populating the DB.  
* REF_RANGE_METHOD_NAME - Name of the reference ranges you wish to run on.  It is highly suggested you use "refRegionGroup" as this is the genic regions.  This approach does not work very well with the intergenes.
* PATH_DIR_NAME - Directory holding the paths. This was filled in the FindPaths step.
* KNOWN_SITES_VCF - Optional parameter to force certain sites to be exported.  If you have a known SNP set which you need all the sites to compare, pass that SNP set in as a vcf file.  If not, remove the -positionVCF flag as well.
* REFERENCE_FILE - reference fasta used when creating the graph.
* OUTPUT_VCF_FILE - name of the output vcf file.