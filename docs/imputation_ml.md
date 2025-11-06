# Imputation using Machine Learning


!!! danger "Disclaimer"
    This Imputation using Machine Learning (ML) section is a work in progress. 
    It is not yet complete and may contain inaccuracies or incomplete 
    information. 


In this document, we will discuss the steps needed to build PS4G files 
using the full assemblies and then to perform Machine Learning Based 
imputation using the PHG:

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

  !!! note "Getting gVCF files without reference ranges"
      When not using reference ranges (skipping `phg create-ranges`), you cannot
      run the standard `phg create-ref-vcf` and `phg create-maf-vcf` commands.
      Instead, you can convert MAF alignment files directly to gVCF format using
      the `maf-to-gvcf-converter` command from [biokotlin-tools](https://github.com/maize-genetics/biokotlin-tools?tab=readme-ov-file#1-maf-to-gvcf-converter).

* Find maximal exact matches between ropebwt3 index and reads:
  ```shell
  phg align-reads \
      --index ropebwt.fmd \
      --query query.fastq \
      --min-smem-len 31 \
      --threads 4 \
      --output example.bed
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


!!! note
    Be sure to include your reference fasta file in the `--keyfile` as well.

!!! note
    The SampleName should not have any underscores in it. We rename the 
    contigs in the assembly fasta file temporarily so that RopeBWT3 can 
    differentiate between contigs of the same name but coming from 
    different assemblies (e.g., `chr1`, `chr2`, ... etc.).


<br>
<hr/>


### Build spline knots

> Build spline knot points from gVCF or hVCF data

**Command** - `build-spline-knots`

**Example**

```bash
build-spline-knots \
  --vcf-dir /data/hvcf_files \
  --output-dir /results/spline_knots \
  --vcf-type hvcf \
  --num-bps-per-knot 75000 \
  --random-seed 42
```

**Parameters**

| Parameter name       | Description                                                                                                                                                                               | Default value | Required?        |
|----------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|---------------|------------------|
| `--vcf-dir`          | Directory containing the hVCF or gVCF files.                                                                                                                                              | `""`          | :material-check: |
| `--vcf-type`         | Type of VCFs to build the splines from. Accepts `"hvcf"` or `"gvcf"`.                                                                                                                     | `"hvcf"`      |                  |
| `--output-dir`       | Output directory to write the spline knots to.                                                                                                     | `""`          | :material-check: |
| `--min-indel-length` | Minimum length of an indel to break up the running block for spline creation of gVCFs. Ignored if `--vcf-type` is `hvcf`.                                                                 | `10`          |                  |
| `--num-bps-per-knot` | Maximum number of base pairs per knot for each contig’s spline. The actual number may be lower if a contig has fewer bases.                                                               | `50000`       |                  |
| `--contig-list`      | Comma-separated list of chromosomes to include in spline generation. If not provided, all chromosomes will be included.                                                                    | `""`          |                  |
| `--random-seed`      | Random seed used for downsampling the number of points per chromosome. Ensures reproducibility.                                                                                           | `12345`       |                  |

<br>
<hr/>

### Align reads to ropebwt3 index

> This command serves as a high-level wrapper around the
> [`ropebwt3 mem`](https://github.com/lh3/ropebwt3?tab=readme-ov-file#finding-maximal-exact-matches) 
> algorithm, providing a streamlined interface for aligning short 
> sequencing reads to a pre-built **FM-index** reference.
> 
> This command identifies **Super-Maximal Exact Matches (SMEMs)** between 
> query reads and the indexed reference, efficiently mapping reads to genomic 
> coordinates. Results are exported in **BED** format for easy visualization 
> and downstream PS4G construction (_see next section_).

**Command** - `align-reads`

**Example**

```shell
phg align-reads \
    --index ropebwt.fmd \
    --query query.fastq \
    --min-smem-len 31 \
    --threads 4 \
    --output example.bed
```

**Parameters**

| Parameter name     | Description                                                                                                                                     | Default value | Required?        |
|--------------------|-------------------------------------------------------------------------------------------------------------------------------------------------|---------------|------------------|
| `--index`          | Path to the [ropebwt3](https://github.com/lh3/ropebwt3) **FM-index** (`.fmd`) file used as the reference for read alignment.                    | `""`          | :material-check: |
| `--query`          | Input FASTQ file containing reads to align. Supports both uncompressed and `.gz`-compressed files.                                              | `""`          | :material-check: |
| `--min-smem-len`   | Minimum **SMEM (Super-Maximal Exact Match)** length used as a seed for alignment. Larger values produce fewer but more specific matches.        | `31`          |                  |
| `--threads`        | Number of threads to use for parallel processing during alignment. Improves performance on multicore systems.                                   | `1`           |                  |
| `--output`         | Output file in **BED** format containing aligned read intervals, mapping quality, and strand information.                                       | `""`          | :material-check: |

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

**Command** - `convert-ropebwt2-ps4g-file`

**Example**

```shell
phg convert-ropebwt2-ps4g-file \
    --ropebwt-bed /data/alignments/sample.mem.bed \
    --output-dir /results/ps4g \
    --spline-knot-dir /refs/spline_knots \
    --min-mem-length 148 \
    --max-num-hits 50 \
    --sort-positions
```

**Parameters**

| Parameter name        | Description                                                                                                                                                                                                 | Default value | Required?        |
|-----------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|---------------|------------------|
| `--ropebwt-bed`       | Path to the RopeBWT3 **BED** file (MEM hits per read) to convert into PS4G.                                                                                                                                 | `""`          | :material-check: |
| `--output-dir`        | Output directory where the generated **PS4G** file will be written.                                                                                                                                         | `""`          | :material-check: |
| `--spline-knot-dir`   | Directory containing **Spline Knot** lookup files (per contig / sample) used to transform MEM positions into consensus PS4G positions.                                                                      | `""`          | :material-check: |
| `--min-mem-length`    | Minimum length of a match (MEM) to be considered for consensus. Shorter MEMs are ignored.                                                                                                                   | `148`         |                  |
| `--max-num-hits`      | Maximum total number of hits allowed (sum across best MEMs). If a read hits more haplotypes than this, it is ignored.                                                                                       | `50`          |                  |
| `--sort-positions`    | Sort positions in the resulting PS4G file. Use `--no-sort-positions` to disable.                                                                                                                            | `true`        |                  |

!!! note
    ropebwt3 can hit more than the value provided in the
    `--max-num-hits` parameter but any alignment hitting more
    haplotypes than this will be ignored.

