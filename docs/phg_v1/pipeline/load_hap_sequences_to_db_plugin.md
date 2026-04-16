!!! warning "Legacy Documentation - PHG Version 1"

    This section contains documentation for **PHG version 1**, which is
    no longer actively developed. It is preserved here for archival and
    historical reference only. If you are looking to use the Practical
    Haplotype Graph, please refer to the [PHG v2 documentation](../../index.md),
    which reflects the current version of the software.

This plugin takes fasta files created from FilterGVCFPlugin or assembly processing and loads PHG data base.  The fasta id lines contain information identifying the reference genome interval of which the sequence represents, along with the name of the corresponding gvcf file.  The gvcf file is used to extract and stored variant contexts for the sequence.  The fasta id line must be formatted as below:

```
#!java
>refChrom:refstart:refEnd gvcfFileName

```

When running the PHG pipeline from the Docker image with the supplied scripts, this plugin is called directly from the FilterGVCFPlugin (for haplotypes) or the Minimap2Plugin (for assembly sequences).  The input files are created by these plugins.

The parameters to the LoadHapSequenceToDBPlugin  method are:

* -fasta <Fasta> Fasta file containing haplotype sequences.  (Default = null) (REQUIRED)
* -genomeData <Genome Data> Path to tab-delimited file containing genome speciic data with header line:\nGenotype Hapnumber Dataline Ploidy Reference GenePhased ChromPhased Confidence Method MethodDetails RefVersion. (Default = null) (REQUIRED)
* -gvcf <Gvcf> Directry containing GVCF file used to create the haplotype fasta file.  Directory path only including trailing /. This directory name is pre-pended to the gvcf file name that appears in the fasta id line for each sequence.

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



