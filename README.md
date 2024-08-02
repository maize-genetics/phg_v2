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

_Long-form documentation for this section can be found 
[here](docs/build_and_load.md). Additional information about **QC 
metrics** can be found [here](docs/qc_metrics.md)._

> [!NOTE]
> As of version 2.4.X, the phg utilizes a new version of anchorwave(1.2.3).
> This changes how ASM coordinates are handled. 
> If you are using old MAF files generated either from anchorwave 1.2.2 or from PHGv2 version 2.3 or eariler, 
> please use the --legacy-maf-file flag for create-maf-vcf.  More information can be found [here](docs/build_and_load.md).

```shell
## Setup conda environment
./phg setup-environment

## Initialize TileDB DataSets
./phg initdb --db-path /path/to/dbs

## Preprocessing data
./phg prepare-assemblies --keyfile /path/to/keyfile --output-dir data/updated_assemblies --threads numberThreadstoRun

## Build VCF data
./phg create-ranges --reference-file data/updated_assemblies/Ref.fa --gff my.gff --boundary gene --pad 500 --range-min-size 500 -o /path/to/bed/file.bed
./phg align-assemblies --gff anchors.gff --reference-file data/updated_assemblies/Ref.fa --assembly-file-list assembliesList.txt --total-threads 20 --in-parallel 4 -o /path/for/generatedFiles
./phg agc-compress --db-path /path/to/dbs --reference-file data/updated_assemblies/Ref.fa --fasta-list /my/assemblyFastaList.txt 
./phg create-ref-vcf --bed /my/bed/file.bed --reference-file data/updated_assemblies/Ref.fa --reference-url https://url-for-ref --reference-name B73 --db-path /path/to/tiled/dataset folder
./phg create-maf-vcf --db-path /path/to/dbs --bed /my/bed/file.bed --reference-file data/updated_assemblies/Ref.fa --maf-dir /my/maf/files -o /path/to/vcfs

## OPTIONAL: Convert GVCF to HVCF: use this instead of create-maf-vcf if you have GVCF files created by PHG, but do not have MAF or h.vcf files
./phg gvcf2hvcf --bed /my/bin/file.bed --gvcf-dir /my/gvcf/dir --reference-file data/updated_assemblies/Ref.fa --db-path /path/to/dbs
 
## Load data into DBs
./phg load-vcf --vcf /my/vcf/dir --dbpath /path/to/dbs
```

### Imputation

_Long-form documentation for this section can be found [here](docs/imputation.md)_


```shell
## Export
./phg export-vcf --db-path /my/db/uri --dataset-type hvcf --sample-names LineA,LineB --output-dir /my/hvcf/dir

## Index
./phg build-kmer-index --db-path /my/db/uri --hvcf-dir /my/hvcf/dir

## Map
./phg map-kmers --hvcf-dir /my/hvcf/dir --kmer-index /my/hvcf/dir/kmerIndex.txt --key-file /my/path/keyfile --output-dir /my/mapping/dir

## Find paths (impute)
./phg find-paths --path-keyfile /my/path/keyfile --hvcf-dir /my/hvcf/dir --reference-genome /my/ref/genome --path-type haploid --output-dir /my/imputed/hvcfs

## Load in DB
./phg load-vcf --vcf /my/imputed/hvcfs --dbpath /my/db/uri
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
3. [Imputation](docs/imputation.md)
4. [Export data](docs/export_data.md)

### Reference
* [haplotype region handling](docs/hvcf_region_handling.md)
* [hVCF format specifications](docs/hvcf_specifications.md)
* [Ktor specifications](...)
* [PHGv2 terminology](docs/terminology.md)
* [PHGv2 architecture](docs/img/architecture/phg_v2_architecture_20240411.svg)
* [QC metrics](docs/qc_metrics.md)
* [SLURM Usage with `align-assemblies`](docs/slurm_usage.md)