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
    ./phg setup-environment --env-file phg_environment.yml
    ```
* Initialize TileDB instances:
    ```shell
    ./phg initdb --db-path /path/to/dbs
    ```
* Create BED file from GFF for reference range coordinates:
    ```shell
    ./phg create-ranges \
        --reference-file /my/ref.fasta \
        --gff my.gff \
        --boundary gene \
        --pad 500 \
        -o /path/to/bed_file.bed
    ```
* Align assemblies:
    ```shell
    ./phg align-assemblies \
        --gff anchors.gff \
        --reference-file /my/ref.fasta \
        --assemblies assemblies_list.txt \
        --total-threads 20 \
        --in-parallel 4 \
        -o /path/for/generated_files
    ```
* Compress FASTA files
    ```shell
    ./phg agc-compress \
        --db-path /path/to/dbs \
        --reference-file \
        --fasta-list /my/assembly_fasta_list.txt
    ```
* Create VCF files
    ```shell
    # Reference VCF
    ./phg create-ref-vcf \
        --bed /path/to/bed_file.bed \
        --reference-file /my/ref.fasta \
        --reference-name B73 \
        -o /path/to/ref_vcf.vcf
  
    # MAF alignments VCF
    ./phg create-maf-vcf \
        --db-path /path/to/dbs \
        --bed /path/to/bed_file.bed \
        --reference-file /my/ref.fasta \
        --maf-dir /my/maf/files \
        -o /path/to/vcfs
    ```
* Load data into DBs
    ```shell
    ./phg load-vcf --vcf /my/vcf/dir --db-path /path/to/dbs
    ```

## Detailed walkthrough

### Preamble
For the following steps, we will first make an example directory to
house our example input data. The overall structure of the directory
looks like the following

```shell
$ tree phgv2_example
...
```

We will also assume that the PHGv2 application is placed in your 
system's `PATH` variable (_see installation documentation for further 
details_).

### Set up Conda environment
Once you have downloaded the latest release of PHGv2 and have Conda 
installed (_see installation documentation_), you will first need to set up a 
Conda environment to accommodate the necessary pieces of software for 
alignment and storage. Here is a list of essential pieces of
software that will be installed in your Conda environment:

| Software   | Purpose                                                    |
|------------|------------------------------------------------------------|
| agc        | Performant FASTA genome compression                        |
| AnchorWave | Sensitive aligner for genomes with high sequence diversity |
| bcftools   | Utilities for indexing VCF data                            |
| samtools   | bgzip compression for VCF data                             |
| TileDB     | Performant storage core for array data                     |
| TileDB-VCF | API for storing and querying VCF data                      |

Instead of setting this up manually, we can use the 
`setup-environment` command which will automate this process. To
run, use the following command:

```shell
./phg setup-environment --env-file phg_environment.yml
```

This command takes one parameter, `--env-file`. This is a path to the 
Conda environment file. You can create the file 
from scratch (`phg_environment.yml`)  and copy over the contents in 
the following block:

```
name: phgv2-conda
channels:
  - conda-forge
  - bioconda
  - tiledb
dependencies:
  - python=3.8
  - tiledb-py=0.22
  - tiledbvcf-py=0.25
  - anchorwave
  - bcftools
  - samtools
  - agc
```

Another option is to pull in the environment file directly from the 
[GitHub repository](https://github.com/maize-genetics/phg_v2):

```shell
curl https://raw.githubusercontent.com/maize-genetics/phg_v2/main/src/main/resources/phg_environment.yml -o phg_environment.yml
```

After the setup step is complete, we can activate our environment
using the following conda command:

```shell
conda activate phgv2-conda
```

If we look in our example project directory, you will also see two
new logging (`.log`) files which will record all the logging and
error output from the PHGv2 command steps:

* `condaCreate_error.log`
* `condaCreate_output.log`


### Initialize TileDB instances
In PHGv2, we leverage [TileDB](https://tiledb.com/) and 
[TileDB-VCF](https://docs.tiledb.com/main/integrations-and-extensions/genomics/population-genomics) 
for efficient storage and querying of VCF data. For downstream 
pipelines, we will need to set up 2 databases for storing different
VCF information: (1) haplotype and (2) genomic variant information.

Similar to setting
up our Conda environment from the prior section, we can automate this
process using the `initdb` command:

```shell
./phg initdb --db-path vcf_dbs
```

This command takes one parameter, `--db-path` which is the path or
subdirectory to where we want to place our databases. In this
example project, I will set up the databases in a subdirectory called
`vcf_dbs`.

After initialization is complete, we will have two empty TileDB-VCF
instances and a `temp` directory in our `vcf_dbs` subdirectory:

| Directory      | Purpose                                 |
|----------------|-----------------------------------------|
| `gvcf_dataset` | Genomic variant storage                 |
| `hvcf_dataset` | Haplotype variant storage               |
| `temp`         | Creation output and error logging files |

### Create reference ranges
Next, we must define ordered sequences of genic and inter-genic 
ranges across the reference genome. These ordered ranges, which we 
will call "reference ranges", are used to specify haplotype sequence 
coordinate information for the haplotype and genomic variant call 
format data. These genomic positions are structured using the 
[BED](https://en.wikipedia.org/wiki/BED_(file_format)) format.

Generally, this data is not readily available as BED files. PHGv2
can generate this file directly from 
[GFF](https://en.wikipedia.org/wiki/General_feature_format)-formatted 
data using the `create-ranges` command:

```shell
./phg create-ranges \
    --gff data/anchors.gff \
    --reference-file data/Ref.fa \
    --boundary gene \
    --pad 500 \
    -o output/ref_ranges.bed
```

This command uses several parameters:

* `--gff` - Our reference GFF file of interest. Since 
  reference ranges are usually established by genic/intergenic 
  regions, we can leverage coordinate and feature information from
  GFF files.
* `--reference-file` - Our reference genome in 
  [FASTA](https://en.wikipedia.org/wiki/FASTA_format) format. This
  is needed for determining the intergenic regions near the ends of 
  chromosomes.
* `--boundary` - How do you want to define the boundaries of the 
  reference ranges. Currently, these can be defined as either
  `gene` or `cds` regions:

  <img src="img/build_and_load/create_ranges_01.svg" width="500"/>
  
  + In the above figure, if `--boundary` is set to `gene`, the start
    and end positions are at the UTR regions for each gene feature from
    the GFF file, while `cds` will begin and end at the open reading 
    frame. By default, the `cds` option will try to identify the
    transcript with the longest open reading frame if there are
    multiple transcript "children" for a gene "parent" in the GFF
    file.
* `--pad` - The number of basepairs you would like to flank regions:

  <img src="img/build_and_load/create_ranges_02.svg" width="500"/>

  + For example, if we were to set the `--pad` parameter to `500`,
    we would extend the region 500 base pairs upstream and downstream
    of the defined boundary (in this case, `gene`).

> [!NOTE]
> There is a possibility that overlaps of regions will occur. If
> this does happen, `create-ranges` will identify any overlapping
> regions and merge regions together:
> ```
> [---Overlap 1----]
>            [-------Overlap 2-------]
> ====================================
> [---------Merged overlap-----------]
>  ```
  
* `-o` - Name for the output BED file.

In the above example, I am using test data from the 
[PHGv2 GitHub repository](https://github.com/maize-genetics/phg_v2/tree/main/data/test/smallseq).
After the command is finished, we have produced a BED file (which
in my case, is located in a subdirectory labelled `output`). This BED 
file contains 4 columns of information:

| Column | Value               |
|--------|---------------------|
| `1`    | Sequence name       |
| `2`    | Start position (bp) |
| `3`    | End position (bp)   |
| `4`    | Range id            |


### Align assemblies
Next, we must align our collection of genome assemblies against a
single reference genome in order to have a common coordinate system
across all genomes within our PHG databases. While whole genome
alignment can be performed in a multitude of ways, PHGv2 provides an
opinionated wrapper to [AnchorWave](https://doi.org/10.1073/pnas.2113075119), 
which uses an **Anchor**ed **Wave**front alignment approach for 
genomes with high sequence diversity, extensive structural 
polymorphism, or whole-genome duplications. Since this software is
already set up during the Conda environment step, there is no need
to install this manually.

To run the aligner step, we can call the `align-assemblies` command:

```shell
./phg align-assemblies \
    --gff data/anchors.gff \
    --reference-file data/Ref.fa \
    --assemblies data/assemblies_list.txt \
    --total-threads 20 \
    --in-parallel 2 \
    -o output/alignment_files
```

This command uses several parameters:
* `--gff` - GFF file for the reference genome. This is used to
  identify full-length coding sequences to use as anchors
* `--reference-file` - The reference genome in 
  [FASTA](https://en.wikipedia.org/wiki/FASTA_format) format.
* `--assemblies` - A text file containing a list of assembly genomes.
  The contents of this file should be either full or relative paths
  to each assembly you would like to align. For example, since I am
  using the example data found on the 
  [PHGv2 GitHub repository](https://github.com/maize-genetics/phg_v2/tree/main/data/test/smallseq),
  I can create a text file called `assemblies_list.txt` and populate
  it with the following lines:

  ```
  data/LineA.fa
  data/LineB.fa
  ```
  Here, I am planning on aligning two genomes called `LineA` and 
  `LineB`. Since these are located in a subdirectory called `data`
  relative to my working directory, I will also add that to the path.

* `--total-threads` - How many threads would you like to allocate for
  the alignment step? More information about this step and the 
  `--in-parallel` step can be found in the following **Details - 
  threads and parallelization** section.
* `--in-parallel` - How many genomes would like to run in parallel?
  More information about this step and the `--total-threads` step can
  be found in the following **Details - threads and parallelization** 
  section.
* `-o` - The name of the directory for the alignment outputs. 

> [!WARNING]
> The directory that you specify in the output (`-o`) section must
> be a an existing directory.

Once alignment is finished, and you have navigated into the output 
directory (in my case, this would be `output/alignment_files`),
you will see several different file types. In addition to various 
logging files (`.log`) for the AnchorWave procedures, you will notice 
that each assembly will have a collection of different file types:

| File extension | Property                                                                                   |
|----------------|--------------------------------------------------------------------------------------------|
| `.sam`         | [sequence alignment map](https://en.wikipedia.org/wiki/SAM_(file_format)) (SAM) file       |
| `.maf`         | [multiple alignment format](https://genome.ucsc.edu/FAQ/FAQformat.html#format5) (MAF) file |
| `.anchorspro`  | alignment blocks between reference and assembly genomes (used for dotplot generation)      |

The MAF files from this output will be used in the VCF creation step.


#### Details - threads and parallelization

***WIP*** - revisit once we can agree on naming and parameter conventions

<img src="img/build_and_load/align_assemblies_01.svg" width="300"/>

***WIP*** - revisit once we can agree on naming and parameter conventions

<img src="img/build_and_load/align_assemblies_02.svg" width="300"/>


The table below shows the memory usage for a single assembly alignment
based on processor type:

| Processor | peak memory (Gb) | wall time |
|-----------|------------------|-----------|
| SSE2      | 20.1             | 26:47:17  | 
| SSE4.1    | 20.6             | 24:05:07  |
| AVX2      | 20.1             | 21:40:00  |
| AVX512    | 20.1             | 18:31:39  |
| ARM       | 34.2             | 18:08:57  |


### Compress FASTA files
For optimal storage of sequence information, we can convert our 
"plain-text" FASTA files into a more compact representation using 
compression. For PHGv2, we can use a command called `agc-compress`, 
which is a wrapper for the
[Assembled Genomes Compressor](https://github.com/refresh-bio/agc) 
(AGC). AGC provides performant and efficient compression ratios for 
our assembly genomes. Like AnchorWave, AGC is also installed during
the Conda environment setup phase, so there is no need to install 
this manually.

To run the compression step, we can call the `align-assemblies` 
command:

```shell
./phg agc-compress \
    --db-path vcf_dbs \
    --fasta-list data/assemblies_list.txt \
    --reference-file data/Ref.fa
```

This command takes in 3 parameters:
* `--db-path` - path to directory storing the TileDB instances. The
  AGC compressed genomes will be placed here on completion.

> [!NOTE]
> The directory specified here should be the same directory used to 
> initialize the TileDB instances in the database initialization 
> (`initdb`) step.

* `--fasta-list` - List of assembly FASTA genomes to compress.

> [!TIP]
> The list specified in `--fasta-list` can be the same list used
> in the alignment (`align-assemblies`) step.

* `--reference-file` - Reference FASTA genome.

After compression is finished, we can navigate to the directory
containing the TileDB instances. In my case, this would be the
subdirectory, `vcf_dbs`. Here, you will see a new file created:
`assemblies.agc`. This is the compressed AGC file containing our 
assemblies. This file will be used later to query for haplotype
sequence regions and composite genome creation.


### Create VCF files
Now that we have (1) created alignments of our assemblies against a
single reference genome and (2) created compressed representations
of our assembly genomes, we can now create VCF information to
populate our TileDB instances. This process is performed using two
commands:

1. Create hVCF data from reference genome:

```shell
./phg create-ref-vcf \
    --bed output/ref_ranges.bed \
    --reference-file data/Ref.fa \
    --reference-name B73 \
    -o output/vcf_files
```

2. Create hVCF and gVCF data from assembly alignments against reference
   genome:

```shell
./phg create-maf-vcf \
    --db-path vcf_dbs \
    --bed output/ref_ranges.bed \
    --reference-file data/Ref.fa \
    --maf-dir output/alignment_files \
    -o output/vcf_files
```

> [!TIP]
> For more information about the haplotype VCF (hVCF) specification,
> please refer the hVCF specification documentation.






### Load VCF data into DBs