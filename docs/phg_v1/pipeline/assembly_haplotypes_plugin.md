!!! warning "Legacy Documentation - PHG Version 1"

    This section contains documentation for **PHG version 1**, which is
    no longer actively developed. It is preserved here for archival and
    historical reference only. If you are looking to use the Practical
    Haplotype Graph, please refer to the [PHG v2 documentation](../../index.md),
    which reflects the current version of the software.

This plugin takes a reference fasta file containing sequence for a single chromosome and an assembly fasta file with sequence for a single corresponding chromosome.  The fastas are aligned using mummer4 nucmer script with a cluster size of 250 and the --mum parameter (anchor matches are unique in both reference and query).  

The nucmer results are filtered using mummer4 delta-filter script with the -g option.  The -g option provides a 1-1 global alignment not allowing rearrangments.  Mummer script show-coords is run on both the original delta file created by nucmer, and the filtered delta created via delta-filter.  Post-processing is done to add back inversions.

Mummer script show-snps is run using the post-processed coords file.  The list of final-snps is stored in htsjdk VariantContext records.  The haplotype sequence is pulled from the VariantContext records.  For alignments that map to reference anchor regions, haplotype sequence and variant lists are created and loaded to the PHG database.

Both assembly anchor and inter-anchor regions are processed in this method.


Tables populated via this method are:

* genotypes
* gametes
* gamete_groups
* methods
* haplotypes
* gamete_haplotypes

The parameters to this plugin are:

* -ref <Reference Fasta File>  Reference Genome File for a single chromosome.  (Default is null) (REQUIRED)
* -assembly <Assembly Fasta File> Path to assembly fasta for a single chromosome from which to pull sequence for aligning.
* -outputDir <Output directory>  Path to output directory including trailing / that will be used when creating temporary output files.
* -dbConfigFile <DB Config File> Path to config file containing DB parameters host, user, password, DB, type.  Used for making the database connection.  Type must be wither "sqlite" or "postgres" to identify db type for connection.
* -version <version> Version name for the set of DB anchor/inter-anchor coorcinates as stored in the anchor_versions table (Default is null) (REQUIRED)
* -assemblyName <Assembly Name> Name of assembly taxon, to be stored as taxon name in the DB. (Default is null) (REQUIRED)
* -chrom <Chromosome Name> Name of chromosome as it appears both for the reference in the db genome_intervals table and in the fasta file idline for the assembly.  The chromosome name can be just a number, a number preceded by "chr" or a number preceded by "chromsome".  Reference and assembly need to be consistent.

An example config file looks as below:


```
#!java

host=localHost
user=sqlite
password=sqlite
DB=/tempFileDir/outputDir/phgSmallSeq.db
DBtype=sqlite
```