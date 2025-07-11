# Imputation using Machine Learning


!!! danger "Disclaimer"
This Imputation using Machine Learning (ML) section is a work in progress. It is not yet complete and may contain inaccuracies or incomplete information. 


In this document, we will discuss the steps needed to build PS4G files using the full assemblies and then to perform
Machine Learning Based imputation using the PHG:

1. Run ropebwt3 indexing by full length chromosomes
2. Build Splines for the assemblies
3. Align Reads against the Index
4. Create a PS4G file from ropebwt3 BED data
5. Run ML Model - *Future Work*

!!! note
The steps detailed in this document build on the materials from
the "[Building and Loading](build_and_load.md)" documentation.
Please review this if you have not worked with the PHG before!


## Quick start

* Run ropebwt3 indexing by full length chromosomes:
  ```shell
  phg rope-bwt-chr-index \
      --keyfile /my/keyfile.txt \
      --output-dir /my/phgIndex/ \
      --index-file-prefix phg_index \
      --threads 20
  ```
  
* Build Splines for the assemblies:
  ```shell
    phg build-spline-knots \
        --vcf-dir /my/vcfs \
        --vcf-type gvcf \
        --output-dir /my/spline_output_dir \
        --min-indel-length 10 \
        --num-bps-per-knot 50000 \
        --contig-list chr1,chr2,chr3
    ```
  
* Align Reads against the Index:
  ```shell
    phg align-reads \
        --db-path /my/db/uri \
        --data-set hvcf \
        --sample-names /my/sample_names.txt \
        --output-dir /my/alignment_output_dir
  ```
  
* Create a PS4G file from ropebwt3 BED data:
  ```shell
    phg convert-bed-to-ps4g \
        --db-path /my/db/uri \
        --data-set hvcf \
        --sample-names /my/sample_names.txt \
        --output-dir /my/ps4g_output_dir
  ```

## Detailed walkthrough

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


!!!note
Be sure to include your reference fasta file in the `--keyfile` as well.

!!!note
The SampleName should not have any underscores in it.  We rename the contigs in the assembly fasta file temporarily so
that RopeBWT3 can differentiate between contigs of the same name but coming from different assemblies(like chr1, chr2, ... etc.).


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

