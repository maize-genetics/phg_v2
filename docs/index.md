---
hide:
  - navigation
  - toc
---

# PHGv2

The **Practical Haplotype Graph (PHG)** is a powerful tool for 
representing pangenomes. The PHG is optimized for plant breeding and 
genetics, where genomic diversity can be high, phased haplotypes are 
common (e.g., inbred lines), and imputation with low-density markers is 
essential for breeding efficiency. This complements other imputation 
tools (e.g., BEAGLE) designed explicitly for handling samples from 
unphased species characterized by low genetic diversity and high-density 
genotyping.


The PHG is a graph-based trellis representation of consecutive genic 
and intergenic regions (called reference ranges), representing diversity 
across and between samples. It can be used to:

* Create custom genomes for alignment
* Call rare alleles
* Impute genotypes
* Efficiently store genomic data from many samples 
  (i.e., reference, assemblies, and other lines)

The PHG also works well with community standards, including the Breeding 
API (BrAPI) and efficient tools for R, such as [rPHG2](https://rphg2.maizegenetics.net/) for pangenome 
extraction and [rTASSEL](https://rtassel.maizegenetics.net/) for connecting genotype to phenotype.
