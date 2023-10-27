## PHG Architecture
The PHGs are generally used by two type of scientists - one who leads in the building and curating the PHG, and then 
a much larger group of the users that include breeders and genetistists using the data through R or a JupyterHub.

Three permanent data structures are used to represent the pangenomes.  We build off community standard VCF format
and its storage through TileDB. The key difference from standard VCF, which focuses on single base variants, we are
making extensive of the VCF's support for symbolic alternate alleles []



There are three main curation tasks:
* **Build and Load** with the reference and assembly genomes
* **Imputation** of haplotypes with low coverage sequence and genotyping data
* **Rare Allele Discovery** with deep short read sequence data against the known assemblies

These curation tasks are all supported by PHG command line tool, and the results are stored in one folder with several key files and 
data sources:
* Haplotype TileDB stores all the haplotypes, their origin, and the genotypic paths to construct the composite genomes
* Variants TileDB stores the variants called either against the reference genome, or for rare variants called against the closest alternate haplotype
* AGC - Assembled Genomes Compressor
