!!! warning "Legacy Documentation - PHG Version 1"

    This section contains documentation for **PHG version 1**, which is
    no longer actively developed. It is preserved here for archival and
    historical reference only. If you are looking to use the Practical
    Haplotype Graph, please refer to the [PHG v2 documentation](../../index.md),
    which reflects the current version of the software.

This Plugin takes a reference genome fasta, a genome data file and an anchors file.  The anchors file contains columns indicating chrom, start position, end position.  The plugin code grabs sequence from the reference genome fasta based on the coordinates from the anchors file and loads the genome interval sequence to the specified PHG database.  The genome data file contains genome specific details used when adding to the genotypes, gametes and method tables.
 
When finished loading the anchor (conserved) regions, inter-anchor regions are identified and loaded.  For the reference genome, the inter-anchor regions are simply all sequence in regions not included in an anchor region.

Tables populated via this method are:

* genotypes
* gametes
* gamete_groups
* methods
* genome_interval_versions
* genome_intervals
* haplotypes
* gamete_haplotypes

The parameters to this plugin are:

* -ref <Reference Genome File>  Reference Genome File for aligning against. (Default is null) (REQUIRED)
* -anchors <Anchors File> CSV file containing chrom, anchor start positions, anchor end position, gene start, gene end. The positions are physical positions (1-based, inclusive/inclusive).  (Default is null) (REQUIRED)
* -genomeData <Genome Data File> A tab-delimited file containing genome specific data with columns: Genotype Hapnumber Dataline Ploidy Reference GenePhased ChromPhased Confidence Method MethodDetails RefVersion.  (Default is null) (REQUIRED)

The genome data file contents are described below:

* Genotype:  the name of the line as you want it to appear in the db genotypes table "name" column, e.g. "B104" or "B104_haplotype_caller"
* Hapnumber:  The 0-based chromosome number.  For inbreds, this is always 0.
* Dataline: Text description defining the line. This data will be stored in the "description" field of the genotypes table.
* Ploidy:  Number of chromosomes for this genome, the value to be stored in the "ploidy" field of the genome_intervals table.  Should be 1 for the reference.
* Reference:  a boolean field indicating if this genome is the reference (should be true for this plugin)
* GenePhased: a boolean field indicating if the genes are phased.  Should be false for the reference.
* ChromPhased:  a boolean field indicating if the chromosomes are phased.  Should be false for the reference.
* Confidence:  a float field indicating confidence in phasing.  Generally 1 when no phasing.
* Method:  a mathod name by which this version of the reference data can be identified.
* MethodDetails:  Text description defining the method used to create the reference (or indicating the AGP version)
* RefVersion: a version name used to identify the version of this reference data.

This plugin expects a data base connection as input.  This can be obtained by chaining the GetDBConnectionPlugin with the LoadGenomeIntervalsToPHGdbPlugin when called from the command line with the TASSEL run_pipeline.pl script.

An example of chaining these plugins is below:

```
#!python

/tassel-5-standalone/run_pipeline.pl -Xmx10G -debug -GetDBConnectionPlugin -config dbConfigFile -create true  -endPlugin -LoadGenomeIntervalsToPHGdbPlugin -ref reference -anchors rangesFile -genomeData genomeDataFile  -outputDir <pathToOutPut> -endPlugin
```