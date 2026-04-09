!!! warning "Legacy Documentation - PHG Version 1"

    This section contains documentation for **PHG version 1**, which is
    no longer actively developed. It is preserved here for archival and
    historical reference only. If you are looking to use the Practical
    Haplotype Graph, please refer to the [PHG v2 documentation](../../index.md),
    which reflects the current version of the software.

The HapCountBestPathToTextPlugin processes haplotype node data either from inclusion files or from the haplotype_counts table of a PHG database instance.  This plugin loops through all the taxon and their inclusion/exclusion count set.  For each taxon, HapCountBestPathPlugin is called to perform the processing. HapCountBestPathPlugin uses a Hidden Markov Model (HMM) to determine the most likely path for a group of haplotype ids.  It determines the nodes on the path through the graph that is most likely to have generated the set of haplotypes.

When the DB is used to pull haplotype_counts data, all paths created will be stored in the DB.  If an inclusion file is used, paths are NOT stored to the DB.  This is because the paths table needs the haplotype_counts_id.  When an inclusion file is used, the assumption is that these haplotype counts are NOT stored in the DB.

The parameters to this plugin are below, most of which are passed to the HapCountBestPathPlugin::

* -configFile <Config File> Database Configuration file containing properties host,user,password,DB and DBtype where DBtype is either sqlite or postgres (Default=null) (REQUIRED)
* -taxa <Taxa> A comma delimited list of taxa (no spaces allowed) to include in graph. Only nodes containing these taxa will be included in the graph. If no taxa list is supplied, then all taxa in the full graph will be used. (Default=null) (OPTIONAL)
* -inclusionFileDir <Inclusion File Dir> The name of the file containing read inclusion and exclusion counts for hapids. (Default=null) (OPTIONAL)
* -target <Target> The taxon that will be used to evaluate the node list returned. (Default=null) (OPTIONAL)
* -refRangeFile <Ref Range File> The name of the file containing the reference ranges to keep.  Supplied when only a subset of referene ranges are desired for processing. (Default=null) (OPTIONAL)
* -refFileName <Ref File Name> Reference file name to allow indexing on the fly.  (Default=null) (OPTIONAL)
* -outputDir <Output Dir> Output directory for writing files. (Default=null) (REQUIRED)
* -refVersion <Ref Version> Ref version name as stored in the genome_interval_versions table.  Needed for storing data to the haplotype_counts table. (Default=null) (REQUIRED)
* -hapCountMethod <Hap Count Method> Name of method used to creates inclusion/exclusion counts in FastqToHapCountPlugin. (Default=null) (REQUIRED)
* -pMethod <P Method> Name of method to be stored that was used to create paths through the graph. (Default=null) (REQUIRED)