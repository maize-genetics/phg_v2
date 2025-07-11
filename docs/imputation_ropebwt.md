# Imputation using RopeBWT3

!!! note
    This is the recommended imputation method. Please 
    [let us know](https://github.com/maize-genetics/phg_v2/issues/new/choose) if 
    you have any questions or issues.

In this document, we will discuss the steps needed to perform
imputation using the [ropebwt3](https://github.com/lh3/ropebwt3) tool within the PHG. The steps:

1. Export hVCF data
2. Index the pangenome from exported hVCF data
3. Map short reads
4. Find paths

!!! note
    The steps detailed in this document build on the materials from
    the "[Building and Loading](build_and_load.md)" documentation.
    Please review this if you have not worked with the PHG before!

## Quick start

* Get list of samples to export:
  ```shell
  phg list-samples \
      --db-path /my/db/uri \
      --data-set hvcf \
      --output-file /my/sample_names.txt
  ```

* Export hVCF data:
  ```shell
  phg export-vcf \
      --db-path /my/db/uri \
      --dataset-type hvcf \
      --sample-file /my/sample_names.txt \
      --output-dir /my/hvcf/dir
  ```

* Index pangenome:
  ```shell
  phg rope-bwt-index \
      --db-path /my/db/uri \
      --hvcf-dir /my/hvcf/dir \
      --output-dir /my/index/dir \
      --index-file-prefix myindex
  ```

* Map short reads
  ```shell
  phg map-reads \
      --hvcf-dir /my/hvcf/dir \
      --index /my/index/dir/myindex.fmd \
      --key-file /my/path/keyfile \
      --output-dir /my/mapping/dir
  ```

* Find paths
  ```shell
  phg find-paths \
      --path-keyfile /my/path/keyfile \
      --hvcf-dir /my/hvcf/dir \
      --reference-genome /my/ref/genome \
      --path-type haploid \
      --output-dir /my/imputed/hvcfs
  ```

* OPTIONAL: Get SNPs from imputed hVCFs
  ```shell
  phg hvcf2gvcf \
      --hvcf-dir /my/imputed/hvcfs \
      --db-path /my/db/uri \
      --reference-file /my/ref/genome \
      --output-dir /my/gvcf/dir
  ```

## Detailed walkthrough

### hVCF export

!!! note
    This step is currently needed, but will be removed in the future as
    we can do direct connections to the hVCF TileDB instance.

Where we last left off in the "[Build and Load](build_and_load.md)"
steps, we had an example directory that looks similar to the
following:

```
phg_v2_example/
├── data
│   ├── anchors.gff
│   ├── Ref-v5.fa
│   ├── LineA-final-01.fa
│   └── LineB-final-04.fa
├── output
│   ├── alignment_files/
│   ├── ref_ranges.bed
│   ├── updated_assemblies
│   │   ├── Ref.fa
│   │   ├── LineA.fa
│   │   └── LineB.fa
│   └── vcf_files/
└── vcf_dbs
    ├── assemblies.agc
    ├── gvcf_dataset/
    ├── hvcf_dataset/
    ├── hvcf_files/
    ├── reference/
    └── temp/
```

!!! note
    The following steps in the rest of this document will assume we are
    in the top level of our `phg_v2_example` working directory.


For this example imputation pipeline, we will be using haplotypes
generated from two samples in our database:

* `LineA`
* `LineB`

If you do not have [hVCF files](hvcf_specifications.md) for samples
you wish to impute against already generated, we will first need to
export this haplotype data in the form of hVCF data from the TileDB
instance. This is done using the `export-vcf` command:

```shell
./phg export-vcf \
    --db-path vcf_dbs \
    --dataset-type hvcf \
    --sample-file /my/sample_names.txt \
    --output-dir output/vcf_files
```
This command takes in 4 parameters:

* `--db-path` - path to directory storing the TileDB instances. The
  AGC compressed genomes will be placed here on completion.
* `--dataset-type` - the type of dataset to export. In this case, we
  are exporting the hVCF data.
* `--sample-file` - text file with list of sample names to export, 
  one per line.
* `--output-dir` - the directory to place the exported hVCF files.

In our above example, we use the `--sample-file` parameter, but if
your sample list is small, you can use the `--sample-names` parameter.
An example of a sample file is as follows:

```
$ cat sample_names.txt

LineA
LineB
```

Example command using --sample-names parameter:

```shell
./phg export-vcf \
    --db-path vcf_dbs \
    --dataset-type hvcf \
    --sample-names LineA,LineB \
    --output-dir output/vcf_files
```

After the export command, our directory structure will have new
files added to the `vcf_files` directory:

```
phg_v2_example/
├── data
│   ├── anchors.gff
│   ├── Ref-v5.fa
│   ├── LineA-final-01.fa
│   └── LineB-final-04.fa
├── output
│   ├── alignment_files/
│   ├── ref_ranges.bed
│   ├── updated_assemblies
│   │   ├── Ref.fa
│   │   ├── LineA.fa
│   │   └── LineB.fa
│   └── vcf_files
│       ├── LineA.h.vcf *
│       └── LineB.h.vcf *
└── vcf_dbs
    ├── assemblies.agc
    ├── gvcf_dataset/
    ├── hvcf_dataset/
    ├── hvcf_files/
    ├── reference/
    └── temp/
```



### Pangenome indexing

!!! info "Upcoming v2.5 Release Notice"
    We will soon update how indexing is performed for the imputation 
    pipeline by replacing the `rope-bwt-index` command with 
    `rope-bwt-chr-index`. If you wish to try it out now, the indexer 
    is currently usable as a **standalone command**. More information 
    about how to use `rope-bwt-chr-index` can be found **[here](imputation_ml.md#run-ropebwt3-indexing-by-full-length-chromosomes)**.

In order to run the later [ropebwt3 read mapping steps](#read-mapping),
we will first need to **build** a ropebwt3 index. 
This is performed using the `rope-bwt-index` command:

  ```shell
  phg rope-bwt-index \
      --db-path vcf_dbs \
      --hvcf-dir output/vcf_files \
      --output-dir output/index_files \
      --index-file-prefix myindex
  ```

This command has **4** required parameters:

* `--db-path` - path to directory storing the TileDB instances and
  `assemblies.agc` file made in the
  "[Compress FASTA files](build_and_load.md#compress-fasta-files)"
  section of the "[Build and Load](build_and_load.md)" documentation.
* `--hvcf-dir` - the directory containing the hVCF files. This is the
  output directory from the `export-vcf` command.  Right now this is
  required, but will be optional in the future.
* `--output-dir` - the directory to place the index files.
* `--index-file-prefix`  - the prefix for the index file
  name. A number of files will be created all of which start with 
  `--index-file-prefix`.

Running rope-bwt-index creates an index file with the extension `.fmd` 
and some additional files for a sequence suffix array(`.ssa`) and a 
contig length file (`.len.gz`). Both of these are needed to map reads 
to this index. Since we have used defaults, the index files will show 
up in the following `output/index_files` directory of our example:

```
phg_v2_example/
├── data
│   ├── anchors.gff
│   ├── Ref-v5.fa
│   ├── LineA-final-01.fa
│   └── LineB-final-04.fa
├── output
│   ├── alignment_files/
│   ├── ref_ranges.bed
│   ├── updated_assemblies
│   │   ├── Ref.fa
│   │   ├── LineA.fa
│   │   └── LineB.fa
│   └── vcf_files
│   │   ├── LineA.h.vcf
│   │   └── LineB.h.vcf
│   └── index_files
│       ├── myindex.fmd *
│       ├── myindex.fmd.len.gz *
│       └── myindex.fmd.ssa *
└── vcf_dbs
    ├── assemblies.agc
    ├── gvcf_dataset/
    ├── hvcf_dataset/
    ├── hvcf_files/
    ├── reference/
    └── temp/
```


#### Optional parameters
In addition to `--index-file`, this command can take other optional
parameters:

| Parameter name             | Description                                                                                                                                                                 | Default value                             |
|----------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------|-------------------------------------------|
| `--threads`                | Number of threads to use during indexing.                                                                                                                                   | `3`                                       |
| `--delete-fmr-index`       | RopeBWT3 will originally output a .fmr index file.  This is converted to fmd for efficiency.  If this flag is set to false it will keep that original file.                 | `true`                                     |
| `--no-diagnostics` or `-n` | A flag that eliminates the diagnostic report                                                                                                                                | Disabled (report is written)              |


### Read mapping

Now that we have a generated ropeBWT3 index file, we can efficiently
align short reads against the PHG using the `map-reads` command. To
demonstrate this, I will retrieve some example data from our
[PHGv2 GitHub test directory](https://github.com/maize-genetics/phg_v2/tree/main/data/test/kmerReadMapping/simulatedReads).
The files I will be using from this directory for this walkthrough
will be the following:

* `LineA_LineB_1.fq`
* `LineA_LineB_2.fq`

These files are simulated paired-end (e.g., `_1.fq` and `_2.fq`)
short reads with a length of 150 bps in
[FASTQ](https://en.wikipedia.org/wiki/FASTQ_format) format. I will
place these files under the `data` directory in sub folder called
`short_reads`:

```
phg_v2_example/
├── data
│   └── short_reads          *
│   │   ├── LineA_LineB_1.fq *
│   │   └── LineA_LineB_2.fq * 
│   ├── anchors.gff
│   ├── Ref-v5.fa
│   ├── LineA-final-01.fa
│   └── LineB-final-04.fa
├── output
│   ├── alignment_files/
│   ├── ref_ranges.bed
│   ├── updated_assemblies
│   │   ├── Ref.fa
│   │   ├── LineA.fa
│   │   └── LineB.fa
│   └── vcf_files
│   │   ├── LineA.h.vcf
│   │   └── LineB.h.vcf
│   └── index_files
│       ├── myindex.fmd
│       ├── myindex.fmd.len.gz
│       └── myindex.fmd.ssa
└── vcf_dbs
    ├── assemblies.agc
    ├── gvcf_dataset/
    ├── hvcf_dataset/
    ├── hvcf_files/
    ├── reference/
    └── temp/
```

Now that we have both short read data and our k-mer index file, we
can pass these to the `map-reads` command:

```shell  
./phg map-reads \
    --hvcf-dir output/vcf_files \
    --index output/index_files/myindex.fmd \
    --key-file data/key_files/read_mapping_data.txt \
    --output-dir output/read_mappings
```

This command has the following parameters:

* `--hvcf-dir` - the directory containing the hVCF files.
* `--key-file` - a tab-delimited list of FASTQ files for a collection
  of samples.
    + In the above example, I have made a keyfile and placed it in a
      subdirectory under the `data` folder called `key_files`.
    + My example keyfile would look like the following:
      ```
      sampleName  filename  filename2
      LineA_B data/short_reads/LineA_LineB_1.fq  data/short_reads/LineA_LineB_2.fq
      ```
    + If you have more than one sample, you would place additional
      lines at the bottom of the keyfile. For example:
      ```
      sampleName  filename  filename2
      LineA_B data/short_reads/LineA_LineB_1.fq  data/short_reads/LineA_LineB_2.fq
      CrossC data/short_reads/cross_c_1.fq  data/short_reads/cross_c_2.fq
      CrossD data/short_reads/cross_d_1.fq  data/short_reads/cross_d_2.fq
      ```
    + > ℹ️ **Note**  
      The keyfile for this parameter needs **column names**. If you
      are using single-reads, the column names would be:
      `sampleName` and `filename`. If you have paired-end reads (like
      the example above), the column names would be `sampleName`,
      `filename`, and `filename2`.
    + > ℹ️ **Note**  
      File names must be of type "FASTQ". In other words, files must
      end in the permitted extensions: `.fq`, `.fq.gz`, `.fastq`, and
      `.fastq.gz`.
* `--output-dir` - the directory to place the read-mapping files.
* `--index` - the ropeBWT3 index file created by the
  `ropebwt-index` command.



!!! tip
If you do not want to create a keyfile and only have one sample,
you can replace the `--key-file` parameter with the `--read-files`
parameter. This parameter will take the paths to the FASTQ reads as
a comma separated list.

    For example, if we were to modify the prior example, the command
    structure would look like the following:
    
    ```shell
    ./phg map-reads \
       --hvcf-dir output/vcf_files \
       --index output/index_files/myindex.fmd \
       --read-files data/short_reads/LineA_LineB_1.fq,data/short_reads/LineA_LineB_2.fq \
       --output-dir output/read_mappings 
    ```


#### Optional parameters
The Map Reads  command can take the following optional parameters:

| Parameter name     | Description                                                                                                                               | Default value |
|--------------------|-------------------------------------------------------------------------------------------------------------------------------------------|---------------|
| `--threads`        | Number of threads used for mapping                                                                                                        | `5`           |
| `--min-mem-length` | Minimum length of the Maximal Exact Match(MEM) to retain.  Any alignments below this value are ignored                                    | `148`         |
| `--max-num-hits`  | Maximum number of alignments to consider for a given read.  If there are more than `--max-num-hits` alignments, the read will be ignored. | `50`          |




Now that we have performed read mapping, we will have a new
read-mapping file in our example working directory that will be
used for the path finding algorithm in the next step:

```
phg_v2_example/
├── data
│   └── short_reads
│   │   ├── LineA_LineB_1.fq
│   │   └── LineA_LineB_2.fq
│   ├── anchors.gff
│   ├── Ref-v5.fa
│   ├── LineA-final-01.fa
│   └── LineB-final-04.fa
├── output
│   ├── alignment_files/
│   ├── ref_ranges.bed
│   ├── updated_assemblies
│   │   ├── Ref.fa
│   │   ├── LineA.fa
│   │   └── LineB.fa
│   ├── vcf_files
│   │   ├── LineA.h.vcf
│   │   └── LineB.h.vcf
│   └── read_mappings                     *
│       └── LineA_LineB_1_readMapping.txt *
└── vcf_dbs
    ├── assemblies.agc
    ├── gvcf_dataset/
    ├── hvcf_dataset/
    ├── hvcf_files/
    ├── reference/
    └── temp/
```

!!! note
    The naming scheme for the read-mapping files follows the first
    sample in the keyfile or comma-delimited list and:
    1. strips the path and FASTQ extension text
    2. adds the `readMapping.txt` suffix to the stripped sample name

    Additionally, $n$ amount of files will be generated for $n$ amount 
    of short read samples found in the keyfile or comma-delimited list.
    For example, if you have 3 samples you are trying to impute, 3 read 
    mapping files will be generated in the given output. 



### Find Paths
`find-paths` is the command used to impute paths through a
[`HaplotypeGraph`](https://github.com/maize-genetics/phg_v2/blob/main/src/main/kotlin/net/maizegenetics/phgv2/api/HaplotypeGraph.kt)
object based on a set of read-mappings. The method uses read-mapping
files generated by the `map-reads` command. Running `find-paths` uses
the [Viterbi algorithm](https://en.wikipedia.org/wiki/Viterbi_algorithm)
to solve a [hidden Markov model (HMM)](https://en.wikipedia.org/wiki/Hidden_Markov_model)
to identify the set of haplotypes most likely to have generated the
observed read-mappings. Two path types are implemented:

* `haploid` - A haploid path is a **single path** through the graph.
* `diploid` - A diploid path is **two paths** through the graph that
  jointly explain the read-mappings.

!!! note
    `find-paths` can run the old kmer read mapping pipeline.  
    Eventually this will be transitioned over to using the ropebwt3 index 
    and alignment process.  For now run the index and mapping step manually.

#### Input files
`find-paths` can take either read-mapping files generated by the command
`map-kmers` or FASTQ files of sequence reads. If FASTQ files are used,
`map-kmers` is run to generate read-mappings but the intermediate
read-mapping files are not saved. **To save a copy of the read-mapping
files, `map-kmers` and `find-paths` must be run separately.** For
more information, please refer to the prior "[Read Mapping](#read-mapping)"
section for additional input details.

!!! note
    If you using read-mapping files, the files **must** end with the
    `_readMapping.txt` suffix or errors will occur.


#### Example usage
Using our prior working example, we can set up the path finding
commands in the following code block. For reference, I have made
a new subdirectory in the `output` folder called `vcf_files_imputed`
and will place imputed hVCF files there:

```shell
./phg find-paths \
    --path-keyfile data/key_files/path_finding_data.txt \
    --hvcf-dir output/vcf_files \
    --reference-genome output/updated_assemblies/Ref.fa \
    --path-type haploid \
    --output-dir output/vcf_files_imputed
```

This command has the following required parameters:

* `--path-keyfile` **or** `--read-file` **but not both**:
    + `--path-keyfile` - name and path of a tab-delimited keyfile that
      contains a list of the read-mapping files or read files (FASTQ).
        + In the working example, I have made a keyfile and placed it in
          the `key_files` folder under the `data` directory. The keyfile
          looks something like this:
          ```
          sampleName  filename
          LineA_B output/read_mappings/LineA_LineB_1_readMapping.txt
          ```
        + If you wish to skip the prior manual `map-reads` step, the data
          in the keyfile can also be paths to FASTQ data (_see the
          "[Read Mapping](#read-mapping)" section for further details_)
        + > ℹ️ **Note**  
          The keyfile must have two columns labeled `sampleName` and
          `filename` if you are specifying paths to read-mapping files.
          If you are using paths to FASTQ data, you may add another
          column to your keyfile, labelled `filename2`, only if your
          FASTQ data is paired-end.
        + > ℹ️ **Note**  
          The samples in the `sampleName` must be unique and must match
          prior keyfile data from the "[Read Mapping](#read-mapping)"
          steps.
        + > ℹ️ **Note**  
          If an output hVCF for a sample name already exists in the
          output directory, the sample will be skipped.
    + `--read-file` - FASTQ samples or read-mapping files as a
      comma-separated list.
        + If bypassing the manual `map-reads` step, paths to FASTQ
          files can be used. Either 1 (for single-end) or 2 (for
          paired-end) comma-separated files can be input this way.
        + If specifying read-mapping files, comma-separated paths to
          files must have the `_readMapping.txt` suffix to work.
* `--hvcf-dir` - The directory containing the hvcf used to build the
  haplotype graph used for imputation.
* `--reference-genome` -  The name and path to the reference genome
  fasta or fastq file.
* `--output-dir` - The directory where the output hVCFs will be
  written. One file will be written for each sample, and the output
  file name will be `<sample name>.h.vcf`.
* `--path-type`: The type of path to be imputed. The value must be
  either `haploid` or `diploid`


#### Optional parameters
Additional optional parameters may also be set to optimize your
imputation pipeline. Please see the proceeding section
("[Likely Ancestors](#likely-ancestors)") for more information on
how to leverage a set of the following parameters.

| Parameter name             | Description                                                                                                                                                                                                                                                                                                  | Default value       |
|----------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|---------------------|
| `--out-parents-dir`        | The directory where the imputed parents (ancestors) will be written for each sample. File names will be <sampleName>_imputed_parents.txt. If no directory name is supplied, the files will not be written.                                                                                                   | `""`                |
| `--prob-correct`           | The probability that a mapped read was mapped correctly.                                                                                                                                                                                                                                                     | `0.99`              |
| `--prob-same-gamete`       | The probability of transitioning to the same gamete (sample) in the next reference range. This should be equal to 1 - (probability of a recombination). Probability of a recombination can be estimated as the total number of expected recombinations in a sample divided by the number of reference ranges | `0.99`              |
| `--min-gametes`            | The minimum number of gametes with a haplotype in a reference range. Reference ranges with fewer gametes will not be imputed.                                                                                                                                                                                | `1`                 |
| `--min-reads`              | The minimum number of reads per reference range. Reference ranges with fewer reads will not be imputed. If `--min-reads = 0`, all reference ranges will be imputed                                                                                                                                           | `0`                 |
| `--inbreeding-coefficient` | The estimated coefficient of inbreeding for the samples being evaluated. Only used for diploid path type. **The value must be between `0.0` and `1.0`**                                                                                                                                                      | `0.0`               |
| `--max-reads-per-kb`       | ReferenceRanges with more than max-reads-per-kb will not be imputed.                                                                                                                                                                                                                                         | `1000`              |
| `--use-likely-ancestors`   | The value must be `true` or `false`. This indicates whether the most likely ancestors of each sample will be used for path finding.                                                                                                                                                                          | `false`             |
| `--max-ancestors`          | If `--use-likely-ancestors = true`, use at most max-ancestors.                                                                                                                                                                                                                                               | `Integer.MAX_VALUE` |                                                                                                                                                                                                               
| `--min-coverage`           | If `--use-likely-ancestors = true`, use the fewest number of ancestors that together have this proportion of mappable reads. The values must be between `0.5` and `1.0`                                                                                                                                      | `1.0`               |
| `--likely-ancestor-file`   | If `--use-likely-ancestors = true`, a record of the ancestors used for each sample will be written to this file, if a file name is provided.                                                                                                                                                                 | `""`                |
| `--threads`                | The number of threads used to find paths.                                                                                                                                                                                                                                                                    | `3`                 |

#### Likely ancestors
The term "ancestors" refers to the taxa or samples (i.e., assemblies)
used to construct the haplotype graph used for imputation. The term
"ancestors" is used because the model considers the haplotypes in the
graph to represent the potential ancestral haplotypes from which the
samples to be imputed were derived.

In some scenarios, **limiting the number of ancestors** used to
impute a sample will be useful. Some examples include:

* if some samples were derived from an F1 cross and the individuals
  used to build the graph include those parents, using only those
  parents for imputation can improve accuracy.
* Also, because imputing diploid
  paths is more computationally intensive than haploid paths, limiting
  the number of ancestors used per sample may be necessary to control
  the amount of time required to impute each sample.

To restrict the number of ancestors used, set the
`--use-likely-ancestors` parameter to `true`, and provide a value for
either `--max-ancestors` or `--min-coverage`. For the case where the
samples have been derived from a limited number of ancestors setting
`--min-coverage` to `0.95` or `--max-ancestors` to the expected
number of ancestors is a useful strategy. In either case, providing a
name for the output file saves a record of the ancestors used for
each sample and should be checked to make sure the samples behaved as
expected.

When using this method, ancestors are chosen by first counting the
number of reads for a sample that map to each ancestor in the
haplotype graph. The ancestor with the most mapped reads is chosen.
Next the reads not mapping to first chosen ancestor are used to count
reads mapping to the remaining ancestors and the ancestor with the
most reads is selected. This is repeated until either the proportion
of reads mapping to the selected ancestors >= `--min-coverage` or the
number ancestors selected equals `--max-ancestors`. If a name is
provided for the `--likely-ancestor-file`, the chosen ancestors are
written to that file in the order they were chosen along with the
cumulative coverage, which is the proportion of reads mapping to each
ancestor and the previously chosen ones.


#### General output and next steps
Now that we have newly imputed hVCF data, our example working
directory looks like the following:

```
phg_v2_example/
├── data
│   └── short_reads
│   │   ├── LineA_LineB_1.fq
│   │   └── LineA_LineB_2.fq
│   ├── anchors.gff
│   ├── Ref-v5.fa
│   ├── LineA-final-01.fa
│   └── LineB-final-04.fa
├── output
│   ├── alignment_files/
│   ├── ref_ranges.bed
│   ├── updated_assemblies
│   │   ├── Ref.fa
│   │   ├── LineA.fa
│   │   └── LineB.fa
│   ├── vcf_files
│   │   ├── LineA.h.vcf
│   │   └── LineB.h.vcf
│   │   vcf_files_imputed *
│   │   └── LineA_B.h.vcf * 
│   └── read_mappings
│       └── LineA_LineB_1_readMapping.txt
└── vcf_dbs
    ├── assemblies.agc
    ├── gvcf_dataset/
    ├── hvcf_dataset/
    ├── hvcf_files/
    ├── reference/
    └── temp/
```

Similar to the "[Build and Load](build_and_load.md#load-vcf-data-into-dbs)"
steps, we can now load this data into our TileDB instance using the
`load-vcf` command:

```shell
./phg load-vcf \
    --vcf-dir output/vcf_files_imputed \
    --db-path vcf_dbs \
    --threads 10
```

This command takes three parameters:

* `--vcf-dir` - Directory containing hVCF data. Since I have imputed
  data in a specific subdirectory, I will use the path
  `output/vcf_files_imputed`
* `--db-path` - Path to the directory containing the TileDB instances.
* `--threads` - Number of threads for use by the TileDB loading
  procedure.

### Create g.vcf files (OPTIONAL)
Our imputed hVCF files provide data on a haplotype level. If desired we can take
the hVCF files and create gVCF files. This provides SNP level data and is done using
the `hvcf2gvcf` command:

```shell
./phg hvcf2gvcf \
    --hvcf-dir output/vcf_files_imputed \
    --db-path vcf_dbs \
    --reference-file output/updated_assemblies/Ref.fa \
    --output-dir output/gvcf_files
```

This command takes 4 parameters:

* `--hvcf-dir` - Directory containing hVCF data. Since I have imputed
  data in a specific subdirectory, I will use the path
  `output/vcf_files_imputed`
* `--db-path` - Path to the directory containing the TileDB instances.
* `--reference-file` - The reference genome fasta file.
* `--output-dir` - The directory to place the gVCF files.
