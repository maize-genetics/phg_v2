# Imputation


!!! danger "Imputation disclaimer"
    This Imputation method is deprecated. Please use the new "[RopeBWT3 Imputation](imputation_ropebwt.md)" method instead.

    With this imputation method, we are experiencing issues, where some PHGs have low mapping rates (~20% with WGS reads), 
    causing errors in the path-finding steps.
    
    PHGv2 is under active development and we expect to encounter bugs. We aim to report bugs as soon as they are discovered.



In this document, we will discuss the steps needed to perform
imputation using the PHG:

1. Export hVCF data
2. Index k-mers from exported hVCF data
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
  
* Index k-mers:
  ```shell
  phg build-kmer-index \
      --db-path /my/db/uri \
      --hvcf-dir /my/hvcf/dir
  ```

* Map short reads
  ```shell
  phg map-kmers \
      --hvcf-dir /my/hvcf/dir \
      --kmer-index /my/hvcf/dir/kmerIndex.txt \
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
* `--sample-file` - text file with list of sample names to export, one per line.
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



### K-mer indexing

In order to run the later [k-mer read mapping steps](#read-mapping) 
(_and run them quickly with a low memory footprint!_),
we will first need to **build** an index containing all 31 bp length 
[k-mers](https://en.wikipedia.org/wiki/K-mer) identified in
the haplotypes exported from the prior step and store this as a flat
file. This is performed using the `build-kmer-index` command:

```shell
./phg build-kmer-index \
    --db-path vcf_dbs \
    --hvcf-dir output/vcf_files
```


This command has **2** required parameters and **1** optional 
parameter that I will specify for example purposes:

* `--db-path` - path to directory storing the TileDB instances and
  `assemblies.agc` file made in the 
  "[Compress FASTA files](build_and_load.md#compress-fasta-files)"
  section of the "[Build and Load](build_and_load.md)" documentation.
* `--hvcf-dir` - the directory containing the hVCF files. This is the
  output directory from the `export-vcf` command.  Right now this is 
  required, but will be optional in the future.
* `--index-file` (_optional_) - The full path of the k-mer index
  file. If not specified, the default path will be:
  * `<--hvcf-dir input>/kmerIndex.txt`
  * In my case, this would be `output/vcf_files/kmerIndex.txt`

Running build-kmer-index creates an index file and a diagnostic file, kmerIndexStatistics.txt, 
that is written to the same directory as the index file.
Since we have used defaults, a new index file and a statistics file will show in the
following `output/vcf_files` directory of our example:

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
│       ├── kmerIndex.txt *
│       ├── kmerIndexStatistics.txt *
│       ├── LineA.h.vcf
│       └── LineB.h.vcf
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

| Parameter name            | Description                                                                                                                                                                                       | Default value              |
|---------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|----------------------------|
| `--max-hap-proportion`    | Note: This parameter is currently unused.  Will be removed in a future update. Only k-mers mapping to less than or equal to maxHapProportion of haplotypes in a reference range will be retained. | `0.75`                     |
| `--max-arg-length`        | The maximum argument length for a call to the [AGC program](https://github.com/acoleman2000/agc).                                                                                                 | `200000`                   |
| `--no-diagnostics` or `-n` | A flag that eliminates the diagnostic report                                                                                                                                                      | Disabled (report is written) |
| `--max-sample-proportion ` | This parameter will allow for kmers seen below --max-sample-proportion * numSamples to be retained.  Any kmers occurring more than that ratio will be filtered out.                               | `0.5`                      |
| `--use-big-discard-set`    | A flag to use a different data structure to hold the discard set at the expense of RAM and speed.  If the default option fails due to an error consider enabling this flag.                       | Disabled (default data structure is used) |


!!! tip
    If you get an error caused by a call to AGC being too long,
    try reducing the `--max-arg-length` value.

The following optional parameters affect how k-mers are pre-filtered 
to determine which are used for indexing. They would only need to be 
adjusted if the number of k-mers in the index is too low or too high:

| Parameter name  | Description                                                                                                                                                                         | Default value |
|-----------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|---------------|
| `--hash-mask`   | In conjunction with `--hash-filter`, used to mask k-mers for filtering. Default uses only the last k-mer nucleotide. **Only change this value unless you know what you are doing.** | `3`           |
| `--hash-filter` | Only hashes that pass the filter ((hashValue and hashMask) == hashFilter) will be considered. **Only change this value unless you know what you are doing.**                        | `1`           |

#### Kmer index diagnostics

Looking at kmer counts by reference range can help determine if diagnostic kmers provide adequate coverage for mapping 
or if there are reference ranges or genomic regions with inadequate coverage. One potential issue is that similar sequence could be
assigned to adjacent reference ranges in different assemblies instead of to the same range. This could result in the kmers
being discarded because they map to adjacent ranges. The column "adjacentCount" lists, for each reference range, the 
number of kmers discarded because they were detected in the previous reference range.

### Read mapping

Now that we have a generated k-mer index file, we can efficiently 
align short reads against the PHG using the `map-kmers` command. To
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
│       ├── kmerIndex.txt
│       ├── LineA.h.vcf
│       └── LineB.h.vcf
└── vcf_dbs
    ├── assemblies.agc
    ├── gvcf_dataset/
    ├── hvcf_dataset/
    ├── hvcf_files/
    ├── reference/
    └── temp/
```

Now that we have both short read data and our k-mer index file, we
can pass these to the `map-kmers` command:

```shell  
./phg map-kmers \
    --hvcf-dir output/vcf_files \
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
* `--kmer-index` (_optional_) - the k-mer index file created by the
  `build-kmer-index` command.
    + > ℹ️ **Note**  
      The `--kmer-index` parameter should only be used if you
      specify a given name and path for your k-mer index file
      generated in the [last step](#k-mer-indexing). If you have
      chosen the "default" option (i.e., leaving this blank), the 
      `map-kmer` command will automatically detect the 
      `kmerIndex.txt` file found in the directory specified by the 
      `--hvcf-dir` parameter.


!!! tip
    If you do not want to create a keyfile and only have one sample,
    you can replace the `--key-file` parameter with the `--read-files`
    parameter. This parameter will take the paths to the FASTQ reads as
    a comma separated list.
    
    For example, if we were to modify the prior example, the command
    structure would look like the following:
    
    ```shell
    ./phg map-kmers \
       --hvcf-dir output/vcf_files \
       --kmer-index output/kmer_index.txt \
       --read-files data/short_reads/LineA_LineB_1.fq,data/short_reads/LineA_LineB_2.fq \
       --output-dir output/read_mappings 
    ```

!!! note
    If no value is passed to the `--kmer-index` parameter, it will 
    default to `<--hvcf-dir value>/kmerIndex.txt`, the default value 
    used by the `build-kmer-index` [command](#k-mer-indexing). If a 
    non-default value was used for the `--index-file` parameter in 
    `build-kmer-index`, that same value needs to be set here for the 
    `--kmer-index` parameter.


#### Optional parameters
The Map Kmers  command can take the following optional parameters:

| Parameter name                          | Description                                                                                                         | Default value |
|-----------------------------------------|---------------------------------------------------------------------------------------------------------------------|---------------|
| `--threads`                             | Number of threads used for mapping                                                                                  | `5`           |
| `--min-proportion-of-max-count`         | Minimum proportion of the maximum kmer count for a read to be considered a match                                    | `1.0`         |
| `--min-proportion-same-reference-range` | Minimum proportion of the read that must align to the same reference range                                          | `0.9`         |


!!! note
    We have removed `--limit-single-ref-range` as it introduces a very high error rate.
    `map-kmers` will only use kmers from a single reference range for each read.



Now that we have performed k-mer mapping, we will have a new 
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
files generated by the `map-kmers` command. Running `find-paths` uses 
the [Viterbi algorithm](https://en.wikipedia.org/wiki/Viterbi_algorithm) 
to solve a [hidden Markov model (HMM)](https://en.wikipedia.org/wiki/Hidden_Markov_model) 
to identify the set of haplotypes most likely to have generated the 
observed read-mappings. Two path types are implemented: 

* `haploid` - A haploid path is a **single path** through the graph. 
* `diploid` - A diploid path is **two paths** through the graph that 
  jointly explain the read-mappings.


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
        + If you wish to skip the prior manual `map-kmers` step, the data
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
        + If bypassing the manual `map-kmers` step, paths to FASTQ
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
