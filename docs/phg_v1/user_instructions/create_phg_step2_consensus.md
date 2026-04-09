!!! warning "Legacy Documentation - PHG Version 1"

    This section contains documentation for **PHG version 1**, which is
    no longer actively developed. It is preserved here for archival and
    historical reference only. If you are looking to use the Practical
    Haplotype Graph, please refer to the [PHG v2 documentation](../../index.md),
    which reflects the current version of the software.

# Create consensus haplotypes


## Details

The CreateConsensi.sh script is used to create Consensus Haplotypes for each genome interval using a specific set of raw haplotypes. Basically this script will create a haplotype graph and attempt to create consensus haplotypes for each anchor. Once a set of consensus haplotypes are created, the script will upload all consensus haplotypes to the Database. Once this process is done, a graph can be made using the consensus haplotypes.

## Kitchen Sink

Steps in the create consensus pipeline:

1. Create a Haplotype Graph
2. For each anchor do the following:
    1. If the graph does not have the variants, extract the variants from the DB for the raw haplotypes.
    2. Take all the GVCF variants and combine them into a single TASSEL GenotypeTable object.
    3. Remove any Positions covering an indel.
    4. Create consensus clusters and create GenotypeTables for each cluster.
    5. Reintroduce indels which were filtered out if they agree.  If they do not agree, it will follow the indelMergeRule defined in the configuration file.
3. Collect up to 1000 Consensus haplotypes at a time and upload them to the DB.

There are a number of parameters needed for this step. We highly recommend tuning the clustering parameters to match the diversity present in the database you are working with. 

**CreateConsensi parameters** 

* configFile: Path to config file containing DB parameters host, user, password, DB, DBtype. Used for making the database connection. DBtype must be either "sqlite" or "postgres" to identify db type for connection. This file will also contain parameters related to creating the consensus sequences and is described in the *Consensus Parameters* section below. A sample configuration file can be found here: Master Config File
* referenceFasta: Reference fasta file. If running on the command line, this parameter is set only in the config file. If running in the docker this is contained in the mount point /phg/inputDir/reference/.
* haplotypeMethod: Name of the raw haplotypes method from which Consensus haplotypes should be made.  These are method pairs comprised of a haplotype method name and reference range group method name. For example:  You may have  created haplotypes using a method name "gvcfHaplotypes" and want to include only those haplotypes that represent a reference range from the range group of "FocusRange".  The method name and range group would be separated by a comma and would be represented as "gvcfHaplotypes,FocusRange". Multiple method pairs are separated by colons. The range group is optional.  So you could have a method string with 2 pairs e.g. "gvcfHaplotypes,FocusRange:coreHaplotypes,FocusRange:.   You can also omit the range group name.  In this case, haplotypes from all defined reference range groups will be included.  The previous example with no reference range group specified would be "gcvfHaplotypes:coreHaplotypes"
                   
* consensusMethod: Name to be stored in the database denoting this run of CreateConsensi.  This is a user-defined method name.  It is often a name that includes parameters used to create the consensus haplotypes, e.g. "collapseMethod_NAM_CONSENSUS_mxDiv_10ToNeg4".  This name tells me that consensus was created using NAM haplotypes with a mxDiv parameter of 10 to the minus 4.  

**DB Parameters**

* host - host name of the db
* user - db user name
* password - password for db
* DB - name of the db
* DBtype - sqlite or postgres

**Java Arguments**

* Xmx - max heap space for the pipeline. This is similar to Java's -Xmx argument.

**Graph Building Parameters**

* includeVariantContexts - true. This needs to be true to work correctly. Given the new Variant Storage, the graphs variants can fit into memory much better.  Setting it true means the variant info will be included in the graph nodes.
* localGvcfFolder - Folder where the taxon gvcfs are stored.  This is required for consensus processing.  The varianats used to determine consensus will be pulled from the GVCF files.

**Consensus Parameters**

* minFreq - Minimum allele frequency: At each position, if no allele has the minimum frequency, the consensus haplotype allele will be set to missing.
* mxDiv - maximum amount of divergence allowed when clustering. The default is 0.01. It is highly recommended you test this parameter at different levels to determine how much or little haplotype collapse is best for your database. 
* maxClusters - The maximum number of clusters that will be created for a reference range. If mxDiv produces too many clusters then the cut height that produces maxClusters number of clusters will be substituted.
* minSites - minimum number of non-missing sites to be clustered. Any haplotypes with fewer number of non-missing calls will be ignored. The default is 20.
* minCoverage -  For each range, any taxon with coverage of less than this amount will not be used to generate consensus haplotypes and will not be included in any haplotype group for that range.
* minTaxa - Minimum number of taxa, default is 1.
* rankingFile - path to a ranking file. If your haplotypes were created from Assemblies, you need to produce a method ranking the taxon. This file will take the form: taxon\tscore where higher scores mean we trust that taxon more highly.  Do not include a header line. When clustering assemblies, when we have a cluster of similar haplotypes, we choose whichever taxon in that group which has the higher ranking score. To break ties, be sure to give each taxon a different score. One simple way to score things is to count the number of haplotypes covered by each taxon in the DB and use that count as the score. Any other arbitrary ranking can be used.
* clusteringMode - either upgma_assembly(default) or kmer_assembly. Clustering mode "upgma" is supported in PHG versions 0.0.40 or below, but is not supported in PHG version 1.0 or higher.   The differences between the upgma and upgma_assembly is that upgma just builds a tree based on a pairwise calculated distance matrix and then will try to merge haplotypes together.  upgma_assembly will also do a pairwise distance matrix but then will select the better haplotype in the group based on the ranking file specified by the "rankingFile" parameter.   Clustering mode "upgma" is no longer supported in PHG versions supporting variants stored in gvcf files (PHG version 1.0 or higher) as the software now stores a gvcf file id from which each haplotype may be found.  With the upgma method, a haplotype is created by merging regions of different haplotypes together to form a new haplotype.  That new haplotype does not exist in any gvcf file associated with this clustered set of genomes.
* maxThreads - The maximum number of threads to be used to create consensi
* kmerSize - size of kmers for the kmer clustering method (default = 7)
* distanceCalculation - distance calculation type, required when clusteringMode = kmer_assembly

### *Details on running this step through docker*

```
#!bash

WORKING_DIR=/workdir/lcj34/phg_newVersion/
DOCKER_CONFIG_FILE=/phg/config.txt

docker1 run --name consensus_container --rm \
    -v ${WORKING_DIR}/:/phg/ \
    -t maizegenetics/phg:1.0 \
    /tassel-5-standalone/run_pipeline.pl -Xmx100G -debug -configParameters ${DOCKER_CONFIG_FILE} \
        -HaplotypeGraphBuilderPlugin -configFile ${DOCKER_CONFIG_FILE}.txt -methods method1:method2 \ 
                -includeVariantContexts true -endPlugin \
        -RunHapConsensusPipelinePlugin \
                -referenceFasta /phg/inputDir/reference/Zm-B73-REFERENCE-NAM-5.0.fa \
                -dbConfigFile ${DOCKER_CONFIG_FILE}  \
                -collapseMethod NAM_CONSENSUS_mxDiv_10ToNeg4 \
                -collapseMethodDetails NAM_CONSENSUS_mxDiv_10ToNeg4 \
                -rankingFile /phg/rankingFile.txt \
                -mxDiv 0.0001 \
                -clusteringMode kmer_assembly -endPlugin
```

The --name parameter provides a name for the container.  This is optional.

The --rm parameter indicates the container should be deleted when the program finishes executing.  This is optional.

The -v directives are used to mount data from the user machine into the Docker.  The path preceding the ":" is the path on the user machine.  The directory path following the ":" are the paths inside the Docker where the user home directories will be mounted.

The -t directive indicates the Docker image of which this container will be an instance.  The last line tells the Docker container to run the CreateConsensi.sh script which is found in the root directory.  The items following are the parameters to the CreateConsensi.sh script.

### *Files*

**Config file**

An example can be found here: Master config file

**Reference Fasta file**

Fasta file for the reference genome you want to use in the database

**Ranking file**

If your haplotypes were created from Assemblies, you need to produce a method ranking the taxon. This file must be a tab-delimited file of the form: taxon\tscore where higher scores mean we trust that taxon more highly. Do not include a header line. When clustering assemblies, when we have a cluster of similar haplotypes, we choose whichever taxon in that group which has the higher ranking score. To break ties, be sure to give each taxon a different score. One simple way to score things is to count the number of haplotypes covered by each taxon in the DB and use that count as the score. Any other arbitrary ranking can be used.

An example ranking file is below.  In this file, B73 is the best taxon:
```
B73	4
CML103	2
Mo17	1
W22	3
```

### *Plugins*

The RunHapConsensusPlugin goes through each reference range in the database in parallel and invokes software to merge and cluster the haplotype sequence for each range.  The consensus sequences created are then loaded to the database and associated with the method name provided in the plugin parameters.


An example command to string these plugins together and run the create consensus step through TASSEL directly is below the plugin descriptions.

#### **RunHapConsensusPlugin**

The RunHapCollapsePlugin method was written to consolidate calls to the collapse pipeline.  Invocation of this plugin results in the following functionality performed:

Loop through each reference range in the graph:

* Extract the HaplotypeNodes with the VariantContexts (we assume that the user has not pulled these yet for memory reasons)
* Merge all the VariantContext records for each base pair (bp) of the reference range
* Export a GenotypeTable containing each bp
* Run HapCollapse Finding algorithm on this genotype table
* load consensus haplotypes to database haplotypes table



Example command to chain the HaplotypeGraphBuilderPlugin and the RunHapConsensusPipelinePlugin when calling them outside of a docker instance:

The ranking file is needed for both the upgma_assembly and kmer_assembly methods.  The upgma method is supported in PHG versions 0.0.40 and earlier, but not in PHG version 1.0 or later.

```
tassel-5-standalone/run_pipeline.pl -Xmx500G -debug -configParameters config.txt \
        -HaplotypeGraphBuilderPlugin -configFile config.txt -methods method1:method2 \ 
                -includeVariantContexts true -endPlugin \
        -RunHapConsensusPipelinePlugin \
                -referenceFasta /workdir/zrm22/Maize2_0/BuildConsensus/Zm-B73-REFERENCE-NAM-5.0.fa \
                -dbConfigFile config.txt \
                -collapseMethod NAM_CONSENSUS_mxDiv_10ToNeg4 \
                -collapseMethodDetails NAM_CONSENSUS_mxDiv_10ToNeg4 \
                -rankingFile rankingFile.txt \
                -mxDiv 0.0001 \
                -clusteringMode kmer_assembly -isTestMethod true -endPlugin 
```


## Troubleshooting

1. If you are not seeing as many haplotypes as you would like (too much collapse) try setting the mxDiv parameter to a smaller value.
2. If you find you are missing taxa at more reference ranges than expected, try setting the minTaxa parameter to 1. This will keep divergent taxa that do not cluster with at least one other taxon in the database.



[Return to Step 2 pipeline, version 0.0.40 and earlier](create_phg_step1_2_main.md)

[Return to Step 2 pipeline, version 1.0 and greater](create_phg_step2_assembly_and_wgs_haplotypes.md)

[Return to Wiki Home](../home.md)