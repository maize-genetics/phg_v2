# Convenience commands

In addition to the primary commands for the build, imputation, and
resequencing pipelines, PHGv2 also provides a suite of "convenience
commands" for miscellaneous "quality of life (QoL)" improvements. In 
this document, we will discuss the currently available external
commands for performing highly used tasks.


## Convert gVCF files to hVCF files

> Create hVCF files from existing gVCF files created by the PHG

**Command** - `gvcf2hvcf`

**Example**

``` shell
phg gvcf2hvcf \
    --bed my/bed/file.bed \
    --reference-file my/updated/ref/fasta.fa \
    --gvcf-dir gvcf/directory \
    --db-path my/phg/db
```

**Parameters**

| Parameter name       | Description                                                                                                         | Default value                      | Required?        |
|----------------------|---------------------------------------------------------------------------------------------------------------------|------------------------------------|------------------|
| `--bed`              | BED file with entries that define the haplotype boundaries.                                                         | `""`                               | :material-check: |
| `--gvcf-dir`         | Directory containing bgzipped and CSI indexed gVCF files.                                                           | `""`                               | :material-check: |
| `--reference-file`   | Path to local Reference FASTA file.                                                                                 | `""`                               | :material-check: |
| `--conda-env-prefix` | Prefix for the Conda environment to use. If provided, this should be the full path to the Conda environment.        | _Current active Conda environment_ |                  |
| `--db-path`          | Folder name where TileDB datasets and AGC record is stored. If not provided, the current working directory is used. | _Current working dir_              |                  |


## Convert hVCF files to gVCF files

> Create gVCF files from existing hVCF files created by the PHG

**Command** - `hvcf2gvcf`

**Example**

``` shell
phg hvcf2gvcf \
    --reference-file my/updated/ref/fasta.fa \
    --hvcf-dir hvcf/directory \
    --db-path my/phg/db
    --output-dir output/directory/for/gvcfs
```

**Parameters**

| Parameter name       | Description                                                                                                         | Default value                      | Required?        |
|----------------------|---------------------------------------------------------------------------------------------------------------------|------------------------------------|------------------|
| `--hvcf-dir`         | Path to directory holding hVCF files. Data will be pulled directly from these files instead of querying TileDB.     | `""`                               | :material-check: |
| `--reference-file`   | Path to local Reference FASTA file.                                                                                 | `""`                               | :material-check: |
| `--conda-env-prefix` | Prefix for the Conda environment to use. If provided, this should be the full path to the Conda environment.        | _Current active Conda environment_ |                  |
| `--db-path`          | Folder name where TileDB datasets and AGC record is stored. If not provided, the current working directory is used. | _Current working dir_              |                  |
| `--output-dir`       | Output directory for the gVCF files. If not provided, the current working directory is used.                        | _Current working dir_              |                  |


## Merge gVCF files

> Merge multiple gVCF files into a single gVCF file

**Command** - `merge-gvcfs`

**Example**

```shell
phg merge-gvcfs \
    --input-dir my/gvcf/directory \
    --output-file output/merged_gvcfs.g.vcf \
```

**Parameters**

| Parameter name | Description                                | Default value | Required?        |
|----------------|--------------------------------------------|---------------|------------------|
| `--input-dir`  | Path to input gVCF file directory.         | `""`          | :material-check: |
| `--output-dir` | Path and/or filename for merged gVCF file. | `""`          | :material-check: |


## Merge hVCF files

> Merge multiple hVCF files into a single hVCF file

**Command** - `merge-hvcfs`

**Example**

```shell
phg merge-hvcfs \
    --input-dir my/hvcf/directory \
    --output-file output/merged_gvcfs.g.vcf \
    --id-format CHECKSUM \
    --reference-file \
    --range-bed-file 
```

**Parameters**

| Parameter name     | Description                                                                                                   | Default value         | Required?        |
|--------------------|---------------------------------------------------------------------------------------------------------------|-----------------------|------------------|
| `--input-dir`      | Path to input gVCF file directory.                                                                            | `""`                  | :material-check: |
| `--output-dir`     | Path and/or filename for merged gVCF file.                                                                    | `""`                  | :material-check: |
| `--id-format`      | ID format for hVCF files. Options are: `CHECKSUM` or `RANGE_SAMPLE_GAMETE` (_see notes for further details_). | `RANGE_SAMPLE_GAMETE` |                  |
| `--reference-file` | Path to reference FASTA file.                                                                                 | `""`                  |                  |
| `--range-bed-file` | Path to [reference range BED file](build_and_load.md#create-reference-ranges).                                | `""`                  |                  |

!!! note "Note - `id-fomat`"
    If you select `CHECKSUM` for the `--id-format` parameter, the `ID`
    values will be MD5 checksums in the `##ALT` header:

    ```
    ##ALT=<ID=06ae4e937668d301e325d43725a38c3f, ...>
    ```

    If you select `RANGE_SAMPLE_GAMETE`, the `ID` values will change
    to a `reference range/sample/gamete` ID format:

    ```
    ##ALT=<ID=R000001_LineA_G01, ...>
    ```
    