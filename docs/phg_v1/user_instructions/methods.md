!!! warning "Legacy Documentation - PHG Version 1"

    This section contains documentation for **PHG version 1**, which is
    no longer actively developed. It is preserved here for archival and
    historical reference only. If you are looking to use the Practical
    Haplotype Graph, please refer to the [PHG v2 documentation](../../index.md),
    which reflects the current version of the software.

## PHG Methods

The PHG methods are used to define both the source of data and how the data was processed.  Each defined method is stored in the PHG methods table which contains information on the method name, method type and a method description. The method information allows for filtering of data based on the method type stored against the data in specific tables.  Method names and descriptions are determined by the user.  This document provides guidance on how these methods are intended to be used and suggestions on defining methods names.  Ultimately, the method names are determined by the used and should be descriptive of their projects/studies.



### Method Types

The PHG methods are used to define both the source of data and how the data was processed.  Each defined method is stored in the PHG methods table which contains information on the method name, method type and a method description. The method information allows for filtering of data based on the method type stored against the data in specific tables.  Method names and descriptions are determined by the user.  This document provides guidance on how these methods are intended to be used and suggestions on defining methods names.  Ultimately, the method names are determined by the used and should be descriptive of their projects/studies.


The defined methods are:


* ANCHOR_HAPLOTYPES: indicates non-consensus haplotypes stored to the haplotypes table.  This is often the type used to store the reference haplotypes, but can also be used for haplotypes created from Whole Genome Sequencing (WGS) files.
* ASSEMBLY_HAPLOTYPES: used for haplotypes created from assembly genomes that are stored to the haplotypes table.
* CONSENSUS_ANCHOR_SEQUENCE: method used to identify haplotypes created by collapsing similar haplotypes to a single haplotype that is stored to the haplotypes table
* EDGE: currently not used
* READ_MAPPING: used to define the methods used to store data to the read_mapping table.  Reads are mapped to a pangenome using minimap2 or a kmer counting method.  The read mapping method name will often indicate the parameters used to map the reads or information about the origin of the set of reads that are mapped. All entries to the read_mapping table should have a method of type READ_MAPPING.
* PATHS: this method type is used to define a particular set of path finding and read mapping parameters that were used to create the most likely path for a taxon through the haplotype graph.  All entries to the PHG paths table should have a method of type PATHS.
* REF_RANGE_GROUP: this method is used when grouping reference ranges together.  Initially it is used when reference ranges are loaded.  From the initial reference range bed file, a “type” is defined for each reference range.  This type becomes the group associated with that range.  A single reference range may belong to multiple groups when users run the AddRefRangeGroup plugin to subset their ranges.


### Method Usage per Table

There are 4 tables in the PHG that contain a method_id field:  haplotypes table, read_mapping table, paths, and ref_range_ref_range_method table.  The usage of the method_id field for each table is defined below


#### Haplotypes table:

Each haplotype loaded to the haplotypes table is associated with a method.  When creating the haplotype graph, a list of methods is provided.  This allows the graph to contain only those haplotypes defined by a certain method(s).  This is useful, for example, when a user wants to create a graph with only consensus methods.  Or perhaps the user wants only a certain set of known good assemblies in their graph, identified by a particular set of method ids.

In addition, after a genome’s haplotypes have been loaded to the haplotypes table, the software creates a path through the graph of this genome.  This path is defined as all the haplotypes that were inserted into the haplotypes table for this genome.  The path method is given the name &lt;haplotypes_method>_path where “haplotypes_method” is the user provided method name from the plugin used to load the haplotypes.


#### Read_mapping table:

Each read mapping entry is associated with a particular taxon.  The read mappings are created when fastq files are mapped to a pangenome created from a particular set of haplotypes.  The method associated with each mapping indicates the set of reads used and the parameters defined for each run.  There may be multiple read mapping entries for a single taxon.  Each method name defined should be descriptive of what was run to create this set of read_mapping entries.  The read_mapping table methods should be  unique from the haplotypes methods.


#### Paths table:

The paths table holds data created when a haplotype graph and a set of read mappings are used to infer the most likely path through the graph.  When that path is determined, it is stored to the PHG paths table.  The method name associated with this path should provide indication of the read mappings used to create this path.  Often the path method name is defined as &lt;read_mapping_method>_PATH.  The path table methods should be  unique from the haplotypes and read_mapping methods.

**NOTE**: As mentioned above, paths are automatically created for each taxon that is loaded to the haplotypes table.  The method name created for these paths is &lt;hapltoype_method>_PATH, where &lt;haplotype_method> is the method name provided by the user to store the initial haplotypes.  This method is a path method specific to the set of haplotype ids that comprise a specific genome’s path through the graph.  It should not be reused when creating paths based on read_mappings.


#### Ref_range_ref_range_method:

This table groups reference ranges.  It is initially populated when the database is created and reference ranges are stored from a user supplied bed file.  The 4th column of the reference range bed file provides the initial type (ie method) of each reference range, and this is stored as a method in the methods table.  Multiple ranges will have the same method, e.g. “genic” or “inter-genic”.  

Users may create additional groupings of reference ranges when running the AddRefRangeGroupPlugin.  This is useful when there is a desire to group specific ranges based on additional user criteria.  For example:  initial reference range methods may have been “genic” and “non-genic”.  Additional methods could group could put ranges into specific genic regions.  Reference ranges may belong to multiple groups.  All of these groups become methods in the methods table.


### Recommendations for Method Usage

Users should provide a unique method when loading haplotypes, read_mapping and paths.  The same haplotypes method may be used to load haplotypes from multiple genomes.  However a previously defined haplotypes method should not be used as a method to a read_mapping or paths table.  Likewise, a previously defined read_mapping method should not be used as a method to the paths table.