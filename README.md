# PHG version 2
The Practical Haplotype Graph (PHG) is a powerful for representing pangenomes.  The PHG is optimized for the plant breeding and genetics context, where genomic diversity can be high and where imputation with low density markers is essential for breeding.

The PHG is a trellis graph based representation of genic and intergenic regions (called reference ranges or reference intervals) which represent diversity across and between taxa. It can be used to: create custom genomes for alignment, call rare alleles, impute genotypes, and efficiently store genomic data from many lines (i.e. reference, assemblies, and other lines). The PHG also works well with community standards including the Breeding API [BrAPI](https://brapi.org) and powerful tools for R such as [rPHG](https://github.com/maize-genetics/rPHG) for pangenome extraction and [rTASSEL](https://github.com/maize-genetics/rTASSEL) for connecting genotype to phenotype.

[PHGv1](https://bitbucket.org/bucklerlab/practicalhaplotypegraph/wiki/Home) was [published in 2022](https://doi.org/10.1093/bioinformatics/btac410) and solved many of the issues in aligning diverse genomes, efficient storage, and imputing across a pangenome.  However, it relied on a custom relational database that required custom formats and it did not scale well with large numbers of taxa and rare alleles. Additionally, after developing PHGs for six species, we realized there were substantial opportunities to streamline the platform for curation.

# PHGv2 design


# PHG terminoloy


# Example usage
```
## Initialize DBs
./phg initdb /path/to/dbs

## Build VCF data
./phg create-ranges --gff my.gff --boundary gene --pad 500 -o /path/to/bed/file.bed
./anchorwave (ref.fasta, asm.fasta) -o /path/to/maf/files/
./phg build-ref-vcf --bed /my/bed/file.bed --reference /my/ref.fasta -o /path/to/vcfs
./phg build-maf-vcf --maf /my/maf/files -o /path/to/vcfs

## Load data into DBs
./phg load-vcf --vcf /my/vcf/dir --dbpath /my/db/uri
```

```
## Index
./phg index-kmers --ancestor founder.h.vcf -o kmer_index.map // we need this

## Map
./phg map-kmers \
    --kmer-index kmer_index.map \
    --reads my_reads.fastq \ // possibly thousands of samples being inputted
    --output read_count_out.map \ // could we pipe this into impute method? // thousands of outputs
    // consider batch interface here ^^

## Impute
./phg impute \
    --hap-counts read_count_out.map \ // will users understand the di
    --diploid false \
    --ancestor founder.h.vcf \
    --max-anc-hap-num 20 \
    --max-anc-hap-prop 0.95 \
    --output-parent best_parents.txt \
    -o my_impute.h.vcf

## Load
./phg load-vcf --vcf my_impute.vcf --dbpath /my/db/uri
```
