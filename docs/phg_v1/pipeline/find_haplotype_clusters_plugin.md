!!! warning "Legacy Documentation - PHG Version 1"

    This section contains documentation for **PHG version 1**, which is
    no longer actively developed. It is preserved here for archival and
    historical reference only. If you are looking to use the Practical
    Haplotype Graph, please refer to the [PHG v2 documentation](../../index.md),
    which reflects the current version of the software.

The FindHaplotypeClustersPlugin processes a multi-haplotype VCF file from one region, and identifies haplotype clusters to be collapsed into a consensus haplotype.  For every alignment in a genotype table, taxa are filtered by coverage to create clusters.  The clusters are determined based on number of alignments without indels, a specified minimum number of polymorphic sites to calculate identity and a maximum distance from each taxon.  The resulting clusters are filtered again based on a minimum number of taxa in each group.  Those clusters passing the filtering criteria have their sequence pulled and written to a fasta file for loading to the PHG DB haplotypes table.

The parameters for this plugin are:

* -ref <Reference> input Reference fasta for pulling reference sequence where gvcf indicates REFRANGE. (Default=null) (REQUIRED)
* -haplotypeMethod <Haplotype Method> Name of haplotype method to sue when retrieving haplotypes form the graph.  Must be a method currently stored in the PHG db methods table.  (Default=null) (REQUIRED)
* -dbConfigFile <Db Config File> File holding configuration information for connecting to the PHG DB. (Default=null) (REQUIRED)
* -consensusVCFOutputDir <Consensus VCF Output Dir> Directory in which to store the output VCFs from the consensus process. (Default=null) (REQUIRED)
* -consensusFastaOutputDir <Consensus Fasta Output Dir> Directory in which to store the output fastas from the consensus process. (Default=null) (REQUIRED)
* -refVersion <Ref Version> Version name of the reference intervals.  Must match a version name that is stored in the genome_interval_versions DB table. (Default=null) (REQUIRED)
* -collapseMethod <Collapse Method> Name to be stored in the DB for this method used to collapse the haplotypes into consensus sequences. (Default=null) (REQUIRED)
* collapseMethodDetails <Collapse Method Details> Details describing the collapse method. (Default=null) (REQUIRED)


