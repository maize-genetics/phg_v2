!!! warning "Legacy Documentation - PHG Version 1"

    This section contains documentation for **PHG version 1**, which is
    no longer actively developed. It is preserved here for archival and
    historical reference only. If you are looking to use the Practical
    Haplotype Graph, please refer to the [PHG v2 documentation](../../index.md),
    which reflects the current version of the software.

# Haplotype Methods

A haplotype method describes a set of haplotypes. In the database, haplotypes have methods which indicate how a haplotype was created. Haplotypes created from assemblies have a user defined method name with no default. 
 Haplotypes created from WGS have a user assigned method name, which defaults to "GATK_PIPELINE". Consensus haplotypes also have a user assigned method name, such as "CONSENSUS". In the PHG database, reference ranges are assigned to groups. The default initia two groups of reference ranges are named "FocusRegion" and "FocueComplement". Additional groups can be assigned by users. A haplotype method for building a graph can be a haplotype method, a comma-separated list of haplotype methods, a refRange group:haplotype method pair, or a comma-separated list of pairs.