# Convenience Commands

In addition to the primary commands for the build, imputation, and
resequencing pipelines, PHGv2 also provides a suite of "convenience
commands" for miscellaneous "quality of life (QoL)" improvements. In 
this document, we will discuss the currently available external
commands for performing highly used tasks.

<div class="grid cards" markdown>
-   :material-arrow-left-right:{ .lg .middle } __Conversion__

    ---

    Commands for converting one file type to another.

    [:octicons-arrow-right-24: Go to section](#conversion)


-   :material-table-merge-cells:{ .lg .middle } __Merging__

    ---

    Commands for merging various PHG file types.

    [:octicons-arrow-right-24: Go to section](#merging)


-   :material-table:{ .lg .middle } __Statistics__

    ---

    Commands for generating summary information in tabular format.

    [:octicons-arrow-right-24: Go to section](#statistics)


-   :material-flask:{ .lg .middle } __Experimental__

    ---

    Prototype commands that are under construction and may have
    significant changes in future updates.

    [:octicons-arrow-right-24: Go to section](#experimental)
</div>



## :material-arrow-left-right: Conversion

### Convert gVCF files to hVCF files

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

<br>
<hr/>

### Convert hVCF files to gVCF files

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
| `--batch-size`       | Number of sample vcf files to export in a single batch from tiledb                                                  | `5`                                |                  |

<br>
<hr/>

### Create a GFF file from an imputed hVCF file

> Create a _path-specific_ [GFF file](https://en.wikipedia.org/wiki/General_feature_format) 
> from an imputed hVCF file and existing sample GFF files. Useful for
> creating a GFF file that contains only annotations for paths present
> in the imputed hVCF file.


**Command** - `paths-to-gff`

**Example**

```shell
phg paths-to-gff \
    --hvcf-file my/hvcf/file.h.vcf \
    --key-file my/samples/keyfile.txt \
    --output-file output/path_specific.gff
```

**Parameters**

| Parameter name  | Description                                                                                                                                  | Default value | Required?        |
|-----------------|----------------------------------------------------------------------------------------------------------------------------------------------|---------------|------------------|
| `--hvcf-file`   | Path to hVCF file for which the GFF will be created                                                                                          | `""`          | :material-check: |
| `--key-file`    | Path to key file containing 2 tab-delimited columns: <ol><li>Sample ID</li><li> The **full path** to the GFF3 file for that sample</li></ol> | `""`          | :material-check: |
| `--output-file` | Full path to the file where the new GFF3 file will be written                                                                                | `""`          | :material-check: |                                |                  |

!!! note "Advanced API use"
    For advanced users who would like to leverage this GFF-based data
    structure in memory for downstream Kotlin pipelines in contrast
    to handling output files, you may use the following example code:

    ``` kotlin
    import net.maizegenetics.phgv2.utils.loadGFFsToGff3Feature
    import net.maizegenetics.phgv2.utils.makeGffFromHvcf
    
    // Same as CLI parameter inputs
    val keyFile  = "my/samples/keyfile.txt"
    val hvcfFile = "my/hvcf/file.h.vcf"
    
    // Create GFF 'TreeMap' object
    val resTreeMap = loadGFFsToGff3Feature(keyFile)
    
    // Create HTSJDK 'Gff3Feature' set object
    val taxonPathGFF = makeGffFromHvcf(hvcfFile, resTreeMap)
    ```
    
    In the above example, `taxonPathGFF` is an in-memory [HTSJDK](https://github.com/samtools/htsjdk)
    [`Gff3Feature`](https://javadoc.io/doc/com.github.samtools/htsjdk/2.24.1/htsjdk/tribble/gff/package-summary.html) 
    object that can be used for downstream purposes. See the [`PathsToGff`](https://github.com/maize-genetics/phg_v2/blob/main/src/main/kotlin/net/maizegenetics/phgv2/cli/PathsToGff.kt)
    class source code for further details.

<br>
<hr/>

### Create a PS4G file from read mapping data

> Convert [read mapping](imputation.md#read-mapping) data into a
> [PS4G (positional support for gamete) file](ps4g_specifications.md).

**Command** - `convert-rm2ps4g`

**Example**

```shell
phg convert-rm2ps4g \
    --read-mapping-file /path/to/readmapping.txt \
    --output-dir /dir/for/ps4g/output/ \
    --hvcf-dir /path/to/hvcf/files/
```

**Parameters**

| Parameter name        | Description                                                             | Default value | Required?        |
|-----------------------|-------------------------------------------------------------------------|---------------|------------------|
| `--read-mapping-file` | Path to [read mapping](imputation.md#read-mapping) file.                | `""`          | :material-check: |
| `--output-dir`        | Output directory for the generated [PS4G](ps4g_specifications.md) file. | `""`          | :material-check: |
| `--hvcf-dir`          | Directory containing hVCF files.                                        | `""`          | :material-check: |

<br>
<hr/>

### Create a PS4G file from ropebwt3 BED data

> Convert a [ropebwt3 BED](https://github.com/lh3/ropebwt3?tab=readme-ov-file#finding-maximal-exact-matches) 
> file into a [PS4G (positional support for gamete) file](ps4g_specifications.md).

!!! note
    This command will only work with ropebwt3 files where the reads
    are aligned to the whole assembly chromosome using the 
    [`mem`](https://github.com/lh3/ropebwt3?tab=readme-ov-file#finding-maximal-exact-matches)
    command. MEMs (**M**aximal **E**xact **M**atche**s**) are used
    to determine what the optimal mappping is. One downside to this
    approach is that if a SNP is in the middle of the read, the
    mappings will be ignored. We may integrate running this in
    conjunction with ropebwt3's [Burrows-Wheeler Aligner's Smith-Waterman Alignment (BWA-SW)](https://doi.org/10.1093/bioinformatics/btp698)
    approach (i.e., the [`sw`](https://github.com/lh3/ropebwt3?tab=readme-ov-file#local-alignment) 
    command) in a future update.

**Command** - `convert-ropebwt2ps4g`

**Example**

```shell
phg convert-ropebwt2ps4g \
    --ropebwt-bed /path/to/readmapping.txt \
    --output-dir /dir/for/ps4g/output/ \
    --hvcf-dir /path/to/hvcf/files/
```

**Parameters**

| Parameter name     | Description                                                                                                                                             | Default value | Required?        |
|--------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------|---------------|------------------|
| `--ropebwt-bed`    | Path to [ropebwt3 BED](https://github.com/lh3/ropebwt3?tab=readme-ov-file#finding-maximal-exact-matches) file.                                          | `""`          | :material-check: |
| `--output-dir`     | Output directory for the generated [PS4G](ps4g_specifications.md) file.                                                                                 | `""`          | :material-check: |
| `--hvcf-dir`       | Directory containing hVCF files.                                                                                                                        | `""`          | :material-check: |
| `--min-mem-length` | Minimum length of a possible match to be considered a match. Default value is the average length of a short read (150 bp) - 2 bp for possible variance. | `148`         |                  |
| `--max-num-hits`   | Maximum number of hits to report.                                                                                                                       | `50`          |                  |

!!! note
    ropebwt3 can hit more than the value provided in the 
    `--max-num-hits` parameter but any alignment hitting more
    haplotypes than this will be ignored.



## :material-table-merge-cells: Merging

### Merge gVCF files

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

<br>
<hr/>

### Merge hVCF files

> Merge multiple hVCF files into a single hVCF file

**Command** - `merge-hvcfs`

**Example**

```shell
phg merge-hvcfs \
    --input-dir my/hvcf/directory \
    --output-file output/merged_hvcfs.h.vcf \
    --id-format CHECKSUM \
    --reference-file \
    --range-bedfile 
```

**Parameters**

| Parameter name     | Description                                                                                                   | Default value         | Required?        |
|--------------------|---------------------------------------------------------------------------------------------------------------|-----------------------|------------------|
| `--input-dir`      | Path to input hVCF file directory.                                                                            | `""`                  | :material-check: |
| `--output-dir`     | Path and/or filename for merged hVCF file.                                                                    | `""`                  | :material-check: |
| `--id-format`      | ID format for hVCF files. Options are: `CHECKSUM` or `RANGE_SAMPLE_GAMETE` (_see notes for further details_). | `RANGE_SAMPLE_GAMETE` |                  |
| `--reference-file` | Path to reference FASTA file.                                                                                 | `""`                  |                  |
| `--range-bedfile`  | Path to [reference range BED file](build_and_load.md#create-reference-ranges).                                | `""`                  |                  |

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



## :material-table: Statistics

### List Sample names from datasets

> List the sample names from the AGC compressed file, the TileDB gVCF 
> dataset, and/or the TileDB hVCF dataset

**Command** - `list-samples`

**Example**

```shell
phg list-samples \
    --db-path my/phg/db \
    --output-file output/hvcf_samples.txt \
    --data-set hvcf 
```

**Parameters**

| Parameter name       | Description                                                                                                         | Default value                      | Required?        |
|----------------------|---------------------------------------------------------------------------------------------------------------------|------------------------------------|------------------|
| `--db-path`          | Folder name where TileDB datasets and AGC record is stored. If not provided, the current working directory is used. | `""`                               | :material-check: |
| `--output-file`      | Path and/or filename for the samples list file.                                                                     | `""`                               | :material-check: |
| `--data-set`         | Storage from which to pull sample names. Must be one of `all`, `agc`, `gvcf`, `hvcf`                                | `hvcf`                             |                  |
| `--conda-env-prefix` | Prefix for the Conda environment to use. If provided, this should be the full path to the Conda environment.        | _Current active Conda environment_ |

**Example output**

```
B73
Ki11
Mo18W
Ki3
```

<br>
<hr/>

### Create a table of haplotype IDs by reference range

> Creates a tab-delimited table of haplotype IDs by reference range
> coordinates and sample IDs

**Command** - `sample-hapid-by-range`

**Example**

```shell
phg sample-hapid-by-range \
    --input-dir my/hvcf/directory \
    --output-file output/hapid_table
```

**Parameters**

| Parameter name | Description                                  | Default value | Required?        |
|----------------|----------------------------------------------|---------------|------------------|
| `--input-dir`  | Path to directory holding hVCF files.        | `""`          | :material-check: |
| `--output-dir` | Path and/or filename for haplotype ID table. | `""`          | :material-check: |

**Example output**

```
#CHROM START   END   B73                                SEEDGWAS1                          SEEDGWAS10
chr1   1       5000  <c9ecfe3967a71282f3ad7c41d48e0bbf> <b19364bc9a4c07a80986b1ee181446c2> <5c8e72b2e9f11ecc652d5b8e8d0e5bf3>
chr1   5001    6000  <f162e742c4d30f151ae6276fbebe762c> <fdfdaa361c39cf5b6f13fad195d0e519> <283a8261c193212fd5cf43d208673322>
chr1   6001    9000  <471d4abbf0545dede647e65915345648> <d6dd5ecea7fb4e6f77f9e630f601b7a8> <13e0ac1a8d12e1aedd6a5302d1e221fd>
```


### Create a table of haplotype IDs to sample

> Creates a tab-delimited table of haplotype IDs to sample gamete. Can
> be one or multiple samples mapping to each haplotype ID.

**Command** - `hapid-sample-table`

**Example**

```shell
phg hapid-sample-table \
    --hvcf-dir my/hvcf/directory \
    --output-file output/hapid_table
```

**Parameters**

| Parameter name | Description                                  | Default value | Required?        |
|----------------|----------------------------------------------|---------------|------------------|
| `--hvcf-dir`   | Path to directory holding hVCF files.        | `""`          | :material-check: |
| `--output-dir` | Path and/or filename for haplotype ID table. | `""`          | :material-check: |


!!! note "Note - `--hvcf-dir`"
    This is intended for use with hVCF files created from aligning assemblies. While
    this will work with hVCF files from imputation, **all** the sample names will be
    the imputed file name and not the sample names associated with the hapids when
    they were created.

**Example output**

```
"a81a7df7340ae0f14a6dccce0d9632db"  Ki3
"c45452b07db68928da6f4e14d50ba1e3"  Mo18W
"1b6d29dbb7b380e67b15c5a0f0142cf0"  Ms71,R2D2
"a935ee46a1a1118df309fc34bdb9e5a5"  B73,Ky21,Ki11
"b878dec3587e24c4714fec5131d4dbbb"  C3PO
```



## :material-flask: Experimental

### Run ropebwt3 indexing by full length chromosomes

> Creates a ropeBWT3 index for a set of assemblies using the full 
> length indexing method. Each FASTA is taken one at a time and is 
> processed and indexed into the ropeBWT3 index. Once the initial 
> index is finished, the `.fmr` index is converted to `.fmd` and the 
> suffix array is built.

**Command** - `rope-bwt-chr-index`

**Example**

```shell
phg rope-bwt-chr-index \
    --keyfile keyfile.txt \
    --output-dir /path/to/output/bed/files/ \
    --index-file-prefix phgIndex \
    --threads 20 
```

**Parameters**

| Parameter name        | Description                                                                                                                                               | Default value | Required?        |
|-----------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------|---------------|------------------|
| `--keyfile`           | Tab-delimited file containing 2 columns: `Fasta` and `SampleName`. `Fasta` is the full path to the FASTA file, and `SampleName` is the name for assembly. | _None_        | :material-check: |
| `--output-dir`        | Output directory where the index files will be written.                                                                                                   | _None_        | :material-check: |
| `--index-file-prefix` | Prefix for the ropebwt3 index file. This prefix will be added to the output directory and used for generated index files.                                 | _None_        | :material-check: |
| `--threads`           | Number of threads to use for index creation.                                                                                                              | `3`           |                  |
| `--delete-fmr-index`  | Delete the `.fmr` index file after converting to `.fmd`.                                                                                                  | `true`        |                  |
| `--conda-env-prefix`  | Prefix for the Conda environment to use. If provided, this should be the full path to the Conda environment.                                              | `""`          |                  |

!!! note
    * `--keyfile` is a tab-delimited file with 2 columns with names 
      `Fasta` and `SampleName`. `Fasta` needs the **full path** for each 
      assembly FASTA file and `SampleName` needs to be the name you want 
      included in the contig name. The first step of the indexer is to 
      open up each FASTA file and rename the contigs to include the 
      provided sample name separated by an '_' (e.g., `lineA_chr1`).
    * `--index-file-prefix` is the prefix for all the output index files. 
      This tool will make a number of files (some temporary) while it is 
      running each with this prefix. **There should not be an extension 
      here as this will be added as need be**.


<br>
<hr/>

### Initialize custom TileDB instance for hVCFs

> Creates a TileDB array instance to house hVCF header data

**Command** - `init-hvcf-array`

**Example**

```shell
phg init-hvcf-array \
    --db-path /path/for/tiledb/instance/
```

**Parameters**

| Parameter name | Description                                                                                                    | Default value               | Required? |
|----------------|----------------------------------------------------------------------------------------------------------------|-----------------------------|-----------|
| `--db-path`    | Directory name under which TileDB datasets will be created. If this folder does not exist, it will be created. | _Current working directory_ |           |

<br>
<hr/>

### Load hVCF data into custom TileDB instance

> Loads hVCF data into TileDB array instance that was created with
> the command `init-hvcf-array`

**Command** - `load-hvcf`

**Example**

```shell
phg load-hvcf \
    --db-path /path/for/tiledb/instance/ \
    --hvcf-dir /path/to/hvcf/files/
```

**Parameters**

| Parameter name | Description                                                                                                    | Default value               | Required?        |
|----------------|----------------------------------------------------------------------------------------------------------------|-----------------------------|------------------|
| `--db-path`    | Directory name under which TileDB datasets will be created. If this folder does not exist, it will be created. | _Current working directory_ |                  |
| `--hvcf-dir`   | Full path to an hVCF file directory                                                                            | ""                          | :material-check: |

<br>
<hr/>

### Query hVCF arrays

> Query TileDB arrays created from `load-hvcf`

**Command** - `query-hvcf-arrays`

**Example**
```shell
phg query-hvcf-arrays \
    --db-path /path/for/tiledb/instance/ \
    --query-type distinctSamples \
    --array-type variants \
    --ref-range-file /path/to/refrange/bed/file.bed \
    --output-file query_results.txt
```

**Parameters**

| Parameter name     | Description                                                                                                                                                                       | Default value               | Required?        |
|--------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|-----------------------------|------------------|
| `--db-path`        | Directory name under which TileDB datasets will be created. If this folder does not exist, it will be created.                                                                    | _Current working directory_ |                  |
| `--query-type`     | Type of query to perform. Options to choose from are: <ul><li>`distinctSamples`</li><li>`distinctRanges`</li></ul>                                                                | ""                          | :material-check: |
| `--array-type`     | Type of array to query. Options to choose from are: <ul><li>`altHeader`</li><li>`variants`</li></ul>                                                                              | `variants`                  |                  |
| `--ref-range-file` | Full path to a [BED](https://en.wikipedia.org/wiki/BED_(file_format)) file-formatted list of reference ranges to query. If no file is provided, all reference ranges are queried. | _All reference ranges_      |                  |
| `--output-file`    | Name of file for query results.                                                                                                                                                   | ""                          | :material-check: |


