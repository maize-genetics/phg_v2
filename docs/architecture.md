# PHG Architecture

PHGs are utilized by two types of scientists: those who lead in building and curating the PHG, and a larger group of users comprising breeders and geneticists who access the data via R or JupyterHub.

## Data Structures

Three primary data structures represent the pangenomes. We build on the community standard VCF format and its storage in TileDB. The key distinction from the standard VCF, which primarily focuses on single base variants, is our extensive use of the VCF's support for symbolic alternate alleles in haplotype VCF or [h.vcf](hvcf_specifications.md). Additionally, we maintain a second TileDB (variant DB) that stores single base variants for each of the high-quality genomes against either the reference genome or the nearest alternate haplotype.

**Persistent Data Stores**:
* **Haplotype TileDB** - Contains all the haplotypes, their origins, and the genotypic paths to construct the composite genomes.
* **Variants TileDB** - Contains the variants called either against the reference genome or, for rare variants, against the nearest alternate haplotype.
* **AGC** - [Assembled Genomes Compressor](https://github.com/refresh-bio/agc), which stores all assembled genomes in a compressed format.

**Temporary Data Structures**:
* **Kmer Haplotype Index** - An index of kmers across sets of haplotypes, used for short read mapping and imputation.
* **h.vcf** - Haplotype VCF used to record genotypes through a path of haplotypes.
* **g.vcf** - Genomic VCF to record variants against the reference genome or composite genome.
* **fasta** - Stores genomes or haplotypes at the sequence level.
* **maf** - Multiple alignment file used for storing the whole genome alignment of the reference genomes to the alternate genomes.

## Data Flow

Three main curation tasks exist:
* **Build and Load** - Working with reference and assembly genomes.
* **Imputation** - Haplotype imputation using low coverage sequence and genotyping data.
* **Rare Allele Discovery** - Identifying rare alleles with deep short read sequence data against the known assemblies.

These curation tasks are executed using the PHG command line tool.

The haplotype graph can then be accessed either through a server implementing the BrAPI or directly via R accessing the TileDBs.

## API Endpoints
BrAPI

## Deployment
Containers.

## Technology Stack
Mostly in JVM languages (Java and Kotlin) calling BioKotlin, HTSJDK, GATK, 

## Known Limitations
