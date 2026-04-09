!!! warning "Legacy Documentation - PHG Version 1"

    This section contains documentation for **PHG version 1**, which is
    no longer actively developed. It is preserved here for archival and
    historical reference only. If you are looking to use the Practical
    Haplotype Graph, please refer to the [PHG v2 documentation](../../index.md),
    which reflects the current version of the software.

## SCRIPT PURPOSE ##
Generate reference intervals for the Practical Haplotype Graph (PHG)

## NOTES ##
* This script assumes it is run inside the PHG Docker container with predefined I/O paths
* The .gff file is assumed to be in JGI format: gene models have the "gene" name, and an "ID=..." field is present in annotation

## RUNNING THE SCRIPT ##

```
#!bash

docker run --rm                                                                    \
-v /your_data_folder/:/tempFileDir/data                                            \
maizegenetics/phg                                                                  \
/CreateReferenceIntervals.sh -f your_reference.fasta -a your_reference.gene.gff3 [ ... optional parameters]
```


## REQUIRED PARAMETERS ##

```
#!bash

   -f <file name>  
      name of fasta file containing the reference sequence  
   -a <file name>  
      name of genome annotation file in .gff format containing gene model annotation 

```

## OPTIONAL PARAMETERS ##

```
#!bash

  -k <integer>  
     Length of kmer used for determining repetitive regions  
     Default: 11
  -e <integer>  
     Number of bases by which to expand gene models for initial reference interval selection  
     Default: 1000
  -m <integer>  
     Distance (in bp) between genes below which gene models are merged  
     Default: 100
  -p <double>  
     Proportion of kmers to be considered repetitive.  
     This determines the high kmer count tail which is considered repetitive (e.g. the top 0.05 most frequent)  
     Default: 0.1
  -n <integer>  
     Number of kmer copies (genome-wide) above which a kmer is considered repetitive. Overrides -p  
     Default: none, -p is used by default
  -l <integer>  
     The number of bases to consider when evaluating if a location in the genome is repetitive  
     Default: 100
  -s <integer>  
     The step size (in bp) by which to proceed outward from a gene model when evaluating flanking regions  
     Default: 10
```


## SCRIPT RESULTS ##
### Output location (subfolder in the input data folder) ###
* genomic_intervals_*unique-timestamp*

### Relevant output contents ###
*  reference_intervals_run.log -- a log file summarizing parameters for the run
*  *your_fasta*.gene.expand.trimmed.summary_report.tsv -- a summary of seed gene model expansion
*  *your_fasta*.*k*mer_count.tsv -- a complete list of kmer counts for kmers with count > 1
*  *your_fasta*.gene.expand.trimmed.bed -- the final reference intervals, in BED format