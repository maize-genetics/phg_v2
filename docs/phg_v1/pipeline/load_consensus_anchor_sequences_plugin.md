!!! warning "Legacy Documentation - PHG Version 1"

    This section contains documentation for **PHG version 1**, which is
    no longer actively developed. It is preserved here for archival and
    historical reference only. If you are looking to use the Practical
    Haplotype Graph, please refer to the [PHG v2 documentation](../../index.md),
    which reflects the current version of the software.

The LoadConsensusAnchorSequencesPlugin takes an input file in fasta format and adds the sequences to the PHG data base.  The fasta file's idline must be of the format:

```
#!java
> refChrom:refStartPos:refEndPos;TaxaName_hapNumber:TaxaName2_hapNumber:... gvcfFile
```
A semi-colon separates the ref coordinates from the list of taxa.  The ref coordinates are used to identify the anchor to which the sequence belongs.  The taxa names indicate the taxa (haplotypes) whose sequences at the indicated anchor collapse to this consensus.  A space separates the taxaNames from the gvcf file, which is a file name that can be found in the directory indicated by the gvcf file dir parameter.

NOTE:  This will run very slow unless the consensus files are concatenated into 1 per chromosome and put into a separate directory containing on the concatenated files for processing.

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

* -inputFile <Input File> Input fasta file, or directory containing fasta files, with consensus sequences. (Default=null) (REQUIRED)
* -version <Version> Version name for this set of anchors as stored in genome_interval_versions table in db.  This is necessary to pull the genomge_interval ID for storing with each haplotype sequence in the haplotypes table. (Default=null) (REQUIRED)
* -vcfDir <Vcf Dir> Directory that contains the gvcf files for the consensus, including trailing /. (Default=null) (REQUIRED)
* -cmethod <Collapse Method> Name of method used to collapse the anchors.  If it doesn't exist, will be added to DB methods table. (Default=null) (REQUIRED)
* -method_details <Collapse Method Description> Description of methods used to collapse the anchor sequences. (Default=null) (OPTIONAL)

