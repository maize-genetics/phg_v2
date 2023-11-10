# PHGv2 - Building and Loading

In this document, we will discuss the steps needed to:
1. Set up a working Conda environment containing the required
   external software
2. Initialize a [TileDB-VCF](https://docs.tiledb.com/main/integrations-and-extensions/genomics/population-genomics) instance
3. Build VCF data
4. Loading data into TileDB-VCF instances

## Quick start
* Set up the PHGv2 Conda environment:
   ```shell
   ./phg setup-environment
   ```
* Initialize TileDB instances:
    ```shell
    ./phg initdb --db-path /path/to/dbs
    ```
* Create BED file from GFF for reference range coordinates:
    ```shell
    ./phg initdb \
        --gff my.gff \
        --boundary gene \
        --pad 500 \
        -o /path/to/bed/file.bed
    ```
* Align assemblies:
    ```shell
    ./phg align-assemblies \
        --gff anchors.gff \
        --reference-file Ref.fa \
        --assemblies assemblies_list.txt \
        -o /path/for/generated_files
    ```
* Compress FASTA files
    ```shell
    ./phg agc-compress
    ```
* Create VCF files
    ```shell
    ./phg create-ref-vcf
    ./phg create-maf-vcf
    ```
* Load data into DBs
    ```shell
    ./phg load-vcf
    ```

## Setting up Conda
Once you have downloaded and set up all the necessary requirements,
you will first need to set up a Conda environment to accommodate the
necessary pieces of software for alignment and storage:

| Software   | Purpose                                                    |
|------------|------------------------------------------------------------|
| agc        | Performant FASTA genome compression                        |
| AnchorWave | Sensitive aligner for genomes with high sequence diversity |
| bcftools   | Utilities for indexing VCF data                            |
| samtools   | bgzip compression for VCF data                             |
| TileDB     | Performant storage core for array data                     |
| TileDB-VCF | API for storing and querying VCF data                      |



## Setup

## Steps

## Miscellaneous