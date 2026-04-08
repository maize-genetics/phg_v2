!!! warning "Legacy Documentation - PHG Version 1"

    This section contains documentation for **PHG version 1**, which is
    no longer actively developed. It is preserved here for archival and
    historical reference only. If you are looking to use the Practical
    Haplotype Graph, please refer to the [PHG v2 documentation](../../index.md),
    which reflects the current version of the software.

The RunHapCollapsePlugin method was written to consolidate calls to the collapse pipeline.  Invocation of this plugin results in the following functionality performed:

Loop through each reference range in the graph:

* Extract the HaplotypeNodes with the VariantContexts (we assume that the user has not pulled these yet for memory reasons)
* Merge all the VariantContext records for each base pair (bp) of the reference range
* Export a GenotypeTable containing each bp
* Run HapCollapse Finding algorithm on this genotype table

For each VCF file exported, upload them to the DB as a consensus

The parameters to this plugin are:

* -mergedOutputDir <Merged Output Dir> Directory to store the ouput VCFs from the merge process. (Default=null)(REQUIRED)
* -ref <Ref> Input reference fasta to be passed to FindHaplotypeClusters for extracting reference sequence.
* -haplotypeMethod <Haplotype Method> Haplotype calling method used to retrieve the correct haplotypes from the graph. (Default=null) (REQUIRED)
* -dbConfigFile <Db Config File> File holding the database configuration information.
* -consensusVCFOutputDir <Consensus VCF Output Dir> Directory in which to store the output VCFs from the consensus process. (Default=null) (REQUIRED)
* -consensusFastaOutputDir <Consensus Fasta Output Dir> Directory in which to store the output fastas from the consensus process. (Default=null) (REQUIRED)
* -refVersion <Ref Version> Version name of the reference intervals.  Must match a version name that is stored in the genome_interval_versions DB table. (Default=null) (REQUIRED)
* -collapseMethod <Collapse Method> Name to be stored in the DB for this method used to collapse the haplotypes into consensus sequences. (Default=null) (REQUIRED)
* collapseMethodDetails <Collapse Method Details> Details describing the collapse method. (Default=null) (REQUIRED)

This method calls [FindHaplotypeClustersPlugin](find_haplotype_clusters_plugin.md)/[MergeGVCFPlugin](merge_gvcf_plugin.md) in parallel for each reference range to create the haplotype clusters.  The resulting consensus sequences are loaded to the haplotypes table via a call to [LoadConsensusAnchorSequencesPlugin](load_consensus_anchor_sequences_plugin.md).