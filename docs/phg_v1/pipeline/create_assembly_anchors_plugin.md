!!! warning "Legacy Documentation - PHG Version 1"

    This section contains documentation for **PHG version 1**, which is
    no longer actively developed. It is preserved here for archival and
    historical reference only. If you are looking to use the Practical
    Haplotype Graph, please refer to the [PHG v2 documentation](../../index.md),
    which reflects the current version of the software.

This plugin takes an assembly genome fasta file, splits the sequence into contigs based on the presence of N's.  Each contig begins at a defined allele (base of A,C,G or T) and ends when an N is encountered.  N's are skipped, a new contig begins at the next non-N allele.  The contig fasta is aligned against the reference genome using minimap2.  Samtools mpileup command is used to generate a pileup of read bases aligning to the reference, limiting the pileup to the conserved region intervals.  Bcftools are then used to create a gvcf file.  Fastas created from the gvcf are loaded to the db via the LoadHapSequencesToDBPlugin.

Assembly inter-anchor regions are not identified at this time.

Tables populated via this method are:

* genotypes
* gametes
* gamete_groups
* methods
* haplotypes
* gamete_haplotypes

The parameters to this plugin are:

* -ref <Reference Genome File>  Reference Genome File for aligning against. (Default is null) (REQUIRED)
* -genomeFile <Assembly Genome> Path to assembly genome fasta from which to pull sequence for aligning.
* -genomeData <Genome Data File> A tab-delimited file containing genome specific data with columns: Genotype Hapnumber Dataline Ploidy Reference GenePhased ChromPhased Confidence Method MethodDetails RefVersion.  
* -outputDir <Output directory>  Path to output directory including trailing / that will be used when creating temporary output files.
* -intervals <Intervals File> Path to BED formatted file containing the intervals used when creating the reference genome intervals. The positions are 0-based, inclusive/exclusive. This is used in mpileup. (Default is null) (REQUIRED)
* -configFile <DB Config File> Path to config file containing DB parameters host, user, password, DB, type.  Used for making the database connection.  Type must be wither "sqlite" or "postgres" to identify db type for connection.

An example config file looks as below:


```
#!java

host=localHost
user=sqlite
password=sqlite
DB=/tempFileDir/outputDir/phgSmallSeq.db
DBtype=sqlite
```


The genome data file contents are described below:

* Genotype:  the nmae of the line as you want it to appear in the db genotypes table "name" column, e.g. "B104" or "B104_haplotype_caller"
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
