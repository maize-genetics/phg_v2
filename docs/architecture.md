# PHG Architecture
The PHGs are generally used by two type of scientists - one who leads in the building and curating the PHG, and then 
a much larger group of the users that include breeders and genetistists using the data through R or a JupyterHub.

## Data Structures
Three permanent data structures are used to represent the pangenomes.  We build off community standard VCF format
and its storage through TileDB. The key difference from standard VCF, which focuses on single base variants, we are
making extensive of the VCF's support for symbolic alternate alleles in haplotype VCF or [h.vcf](hvcf_specifications.md).
We have a second TileDB (variant DB) that does store single base variants for each of high quality genomes against either 
the reference genome or the closest alternate haplotype.

Persistent data stores:
* Haplotype TileDB - stores all the haplotypes, their origin, and the genotypic paths to construct the composite genomes
* Variants TileDB - stores the variants called either against the reference genome, or for rare variants called against the closest alternate haplotype
* AGC - [Assembled Genomes Compressor](https://github.com/refresh-bio/agc) that stores all assembly genomes in a compressed archive

Temporary data structures:
* Kmer Haplotype Index - an index of kmers across sets of haplotypes used for short read mapping and imputation
* h.vcf - haplotype VCF to record genotypes through a path of haplotypes
* g.vcf - genomic VCF to record variants against reference genome or composite genome
* fasta - stores either genomes or haplotypes at the sequence level
* maf - multiple alignment file for storing the whole genome alignment of the reference genomes to the alternate genomes

## Data Flow
There are three main curation tasks:
* **Build and Load** with the reference and assembly genomes
* **Imputation** of haplotypes with low coverage sequence and genotyping data
* **Rare Allele Discovery** with deep short read sequence data against the known assemblies

These curation tasks are done with the PHG command line tool.

The the haplotype graph are then used through either a server implementing the BrAPI or directly through R accessing the TileDBs.



