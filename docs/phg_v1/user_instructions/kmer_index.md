!!! warning "Legacy Documentation - PHG Version 1"

    This section contains documentation for **PHG version 1**, which is
    no longer actively developed. It is preserved here for archival and
    historical reference only. If you are looking to use the Practical
    Haplotype Graph, please refer to the [PHG v2 documentation](../../index.md),
    which reflects the current version of the software.

# Create a kmer index

The KmerHashMapFromGraphPlugin is a tool that creates a kmer hash map with a specified value of k (in this case, k = 32). It can optionally save this map to a file. Let’s break down its functionality:
  
1. Purpose
 
+ The primary purpose of this plugin is to create a map that associates each kmer hash with the haplotypes in which it occurs.
+ It helps discriminate between different haplotypes based on observed kmers.

2. Steps
 
+ The initial step involves finding a set of kmers that are observed in only one reference range.
+ These kmers occur in at most (maxHaplotypes * number of haplotypes) haplotypes in that reference range.
+ Once a diagnostic kmer set is identified, the plugin creates a map of kmer hash values to the list of haplotypes in which they occur.

3.	Options

+	-maxHaplotypes: Specifies the maximum number of haplotypes as a proportion of total haplotypes for a kmer to be considered.
+	-hashMask: Determines how the positions of a kmer hash are masked for pre-filtering hashes.
+	-hashFilter: Defines which kmers are used based on their masked values (0 = A, 1 = C, 2 = G, 3 = T).
+	-saveFile: Specifies the name and path of the file where the map will be saved.
+	-configFile: Requires a config file containing database connection parameters.

4.	Example Command

> ./run_pipeline.pl -Xmx100G -debug -configParameters /myDir/myConfig.txt -HaplotypeGraphBuilderPlugin -methods myMethod -includeSequences true -endPlugin -KmerHashMapFromGraphPlugin -saveFile myFileName -endPlugin

### Details

Remember to replace the placeholders (/myDir/myConfig.txt, myMethod, myFileName) with actual values relevant to your use case. Setting -Xmx100G should provide enough RAM for most databases. Expect this step to take a few hours. This tool is particularly useful for mapping reads during imputation. 😊
Typing the following command produces help information: ./run_pipeliine.pl -KmerHashMapFromGraphPlugin   
This returns the following:

KmerHashMapFromGraphPlugin Description...  
This plugin creates a kmer hash map with k = 32 and, optionally, saves it to a file. The kmer map is a map of kmer hash to hapid list that associates each kmer hash with haplotypes it occurs in. The first step in creating the map is to find a set of kmers that are observed in only one reference range and that occur in at most maxHaplotypes of them, so that the kmers can discriminate between haplotypes. Once a diagnostic kmer set is found, a map is created of kmer hash value to the list of haplotypes in which they occur. The map is then saved to a file for later use in mapping reads for imputation.

Usage:
KmerHashMapFromGraphPlugin <options>  
-maxHaplotypes <Max Haplotypes> : If a kmer is in more than (maxHaplotypes * number of haplotypes) in a reference range it will not be used. (Default: 0.75)  
-hashMask <Hash Mask> : The value used to mask the positions of a kmer hash used to pre-filter hashes. A value of 3 uses the last position of the kmer, and a value of 15 uses the last two positionns (Default: 3)  
-hashFilter <Hash Filter> : A kmer is used when its masked value equals hashFilter, where 0 = A, 1 = C, 2 = G, and 3 = T. (Default: 1)  
-saveFile <Save File> : The name and path of the file to which the map will be saved.  
-configFile <Config File> : A config file containing database connection parameters. Either this value or a value for -configParameters is required.