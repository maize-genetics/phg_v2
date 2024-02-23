# PHG version 2
[![PHGv2 CI](https://github.com/maize-genetics/phg_v2/actions/workflows/phgv2_ci.yml/badge.svg)](https://github.com/maize-genetics/phg_v2/actions/workflows/phgv2_ci.yml) [![codecov](https://codecov.io/gh/maize-genetics/phg_v2/graph/badge.svg?token=4BVD2QXQ1A)](https://codecov.io/gh/maize-genetics/phg_v2) [![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)

The Practical Haplotype Graph (PHG) is a powerful tool for 
representing pangenomes. The PHG is optimized for plant breeding 
and genetics, where genomic diversity can be high, phased haplotypes 
are common (e.g. inbred lines), and imputation with low density 
markers is essential for breeding efficiency. This complements 
other imputation tools (e.g. [BEAGLE](https://faculty.washington.edu/browning/beagle/beagle.html)) 
designed explicitly for handling samples from unphased species 
characterized by low genetic diversity and high-density genotyping.

The PHG is a trellis graph based representation of consecutive genic 
and intergenic regions (called reference ranges) which represent 
diversity across and between samples. It can be used to:

* Create custom genomes for alignment
* Call rare alleles
* Impute genotypes 
* Efficiently store genomic data from many samples (i.e. reference, 
  assemblies, and other lines)

The PHG also works well with community 
standards including the Breeding API ([BrAPI](https://brapi.org)) and efficient 
tools for R such as [rPHG](https://github.com/maize-genetics/rPHG) for pangenome extraction and 
[rTASSEL](https://github.com/maize-genetics/rTASSEL) for connecting genotype to phenotype.


## Quick start

### Installation

Using a Linux distribution, download the latest release
[here](https://github.com/maize-genetics/phg_v2/releases/latest) or
use the command line:

```shell
curl -s https://api.github.com/repos/maize-genetics/phg_v2/releases/latest \
| awk -F': ' '/browser_download_url/ && /\.tar/ {gsub(/"/, "", $(NF)); system("curl -LO " $(NF))}'
```

Untar and add the wrapper script to your `PATH` variable. Detailed
information about these steps can be found [here](docs/installation.md).
### Build and load data

_Long-form documentation for this section can be found [here](docs/build_and_load.md)_

```shell
## Setup conda environment
./phg setup-environment

## Initialize TileDB DataSets
./phg initdb --db-path /path/to/dbs

## Preprocessing data
./phg annotate-fastas --keyfile /path/to/keyfile --output-dir /path/to/annotated/fastas --threads numberThreadstoRun

## Build VCF data
./phg create-ranges --reference-file Ref.fa --gff my.gff --boundary gene --pad 500 -o /path/to/bed/file.bed
./phg align-assemblies --gff anchors.gff --reference-file Ref.fa -a assembliesList.txt --total-threads 20 --in-parallel 4 -o /path/for/generatedFiles
./phg agc-compress --db-path /path/to/dbs --reference-file /my/ref.fasta --fasta-list /my/assemblyFastaList.txt 
./phg create-ref-vcf --bed /my/bed/file.bed --reference-file /my/ref.fasta --reference-url https://url-for-ref --reference-name B73 --db-path /path/to/tiled/dataset folder
./phg create-maf-vcf --db-path /path/to/dbs --bed /my/bed/file.bed --reference-file /my/ref.fasta --maf-dir /my/maf/files -o /path/to/vcfs
## Convert GVCF to HVCF: use this if you have GVCF files created by PHG, but do not have MAF or h.vcf files
./phg gvcf2hvcf --bed /my/bin/file.bed --gvcf-dir /my/gvcf/dir --reference-file /my/ref.fasta --reference-name B73 
## Load data into DBs
./phg load-vcf --vcf /my/vcf/dir --dbpath /path/to/dbs
```

### Imputation

> [!NOTE]
> This section is currently in progress and command input may be
> subject to change. The following pseudocode is a possible
> representation of the imputation workflow:

```shell
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

### Data retrieval

> [!NOTE]
> This section is currently in progress and command input may be
> subject to change. The following pseudocode is a possible
> representation of the retrieval workflow:

```shell
## Export from Tiledb
./phg export-vcf --db-path /my/db/uri --dataset-type hvcf --sample-Names LineA,LineB --output-dir /my/output/dir
```




## Design and history

[PHGv1](https://bitbucket.org/bucklerlab/practicalhaplotypegraph/wiki/Home) was [published in 2022](https://doi.org/10.1093/bioinformatics/btac410). It addressed many
challenges related to aligning diverse genomes, efficient storage,
and imputation across a pangenome. However, it depended on a custom
relational database that necessitated unique formats, and database
queries did not scale effectively with a large number of taxa and
rare alleles. Moreover, after developing PHGs for six species, we
identified significant opportunities to refine and streamline the
platform for curation.

The redesign leverages the performant TileDB-VCF database, which is
widely used in human genetics for extensive medical applications and
is highly proficient for rapid querying and storage of rare variants.
The PHG is now backed by two TileDB-VCF databases: one for tracking
haplotypes across all samples (`.h.vcf`), and another for tracking
variants relative to either the reference genomes or the closest
haplotype (`.g.vcf`). Our implementation of haplotype encoding in VCF
heavily relies on the VCF ALT haplotype specification defined in
[v4.2](http://samtools.github.io/hts-specs/VCFv4.2.pdf).

Other important things to note:
* **High-quality phased genome assemblies** (or similar) are available to
  initialize the PHG.
* Ancestral haplotypes are aligned to the reference genome for the
  identification of haplotypes.
* All PHG tools rely on public file standards - FASTA, VCF, BCF, BED,
  and MAF.
* We rely on public tools like TileDB, Minimap2, GATK, AnchorWave,
  BioKotlin, and HTSJDK.
* Genotyping with low-density markers is now done using a memory- and
  speed-efficient k-mer approach, followed by pathfinding (imputation)
  with [hidden Markov model](https://en.wikipedia.org/wiki/Hidden_Markov_model) methods. 
* Rare allele discovery with short reads is based on the above path,
  involving short read alignment to the inferred haplotype path
  genome and the GATK haplotype caller.


## Terminology

When describing components used in the PHG, certain terms are used to 
efficiently communicate more complicated ideas. Some common terms you 
may find are:

| Term             | Definition                                                |
|------------------|-----------------------------------------------------------|
| haplotype        | The sequence of part of an individual chromosome.         |
| path             | The phased set of haplotypes that represent a chromosome. |
| reference genome | A genome used for initial alignment and base coordinates. |
| reference range  | A segment of the reference genome.                        |

More commonly used terms can be found [here](docs/terminology.md).


## Long-form documentation

### PHG workflows
1. [Installation](docs/installation.md)
2. [Building and loading](docs/build_and_load.md)
3. [Imputation](docs/imputation.md) (_**WIP**_) 🚧
4. [Export data](docs/export_data.md) (_**WIP**_) 🚧

### Reference
* [hVCF format specifications](docs/hvcf_specifications.md)
* [haplotype region handling](docs/hvcf_region_handling.md)
* [PHGv2 terminology](docs/terminology.md)
* [PHGv2 architecture](docs/img/architecture/phg_v2_architecture-Diagram.drawio.svg)