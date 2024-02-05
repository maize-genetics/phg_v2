# PHGv2 - Exporting Data

In this document, we will discuss general strategies for exporting
data from a PHG database.

> [!NOTE]
> This section will assume you have a pre-existing database loaded
> with haplotype data. See the ["Building and Loading"](build_and_load.md)
> documentation for further information.

## Quickstart

* Export hVCF files from database
```shell
phg export-vcf \
  --db-path /path/to/dbs \
  --dataset-type hvcf \ # can also be 'gvcf'
  --sample-names LineA,LineB \
  -o /path/to/output/directory
```

* Create FASTA files from hVCF data or database
```shell
phg create-fasta-from-hvcf \
  --hvcf-file my_sample.h.vcf \
  --fasta-type composite \ # can also be 'haplotype'
  -o my_sequence.fa
```

* Data retrieval using BrAPI endpoints and [rPHG2](https://maize-genetics.github.io/rPHG2/)

```shell
phg start-server \
  --db-that /path/to/dbs \
  --port 8080
```

``` r
library(rPHG2)

PHGServerCon("localhost", 8080) |>
    readHaplotypeData()
```

## Export hVCF data

## Create FASTA data

## Data retrieval using BrAPI endpoints and rPHG