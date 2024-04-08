# PHGv2 - Imputation

In this document, we will discuss the steps needed to perform
imputation using the PHG:

1. Indexing k-mers
2. Mapping short reads
3. Finding paths

> [!NOTE]
> The steps detailed in this document build on the materials from
> the "[Building and Loading](build_and_load.md)" documentation. 
> Please review this if you have not worked with the PHG before!

## hVCF export

> [!NOTE]
> This step is currently needed, but will be removed in the future as 
> we can do direct connections to the hVCF TileDB instance.

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
    ├── gvcf_dataset # gVCF db storage
    ├── hvcf_dataset # hVCF db storage
    ├── hvcf_files
    │   ├── Ref.h.vcf.gz
    │   └── Ref.h.vcf.gz.csi
    ├── reference
    │   ├── Ref.bed
    │   └── Ref.sam
    └── temp/
```

> [!NOTE]
> The following steps in the rest of this document will assume we are 
> in the top level of our `phg_v2_example` working directory.


For this example imputation pipeline, we will be using haplotypes 
generated from two samples in our database:

* `LineA`
* `LineB`

If you do not have hVCF files for samples you wish to impute against 
already generated, we will first need to export this haplotype data 
in the form of [hVCF](hvcf_specifications.md) data from the TileDB 
instance. This is done using the `export-vcf` command:

```shell
./phg export-vcf \
    --db-path vcf_dbs \
    --dataset-type hvcf \
    --sample-names LineA,LineB \
    --output-dir output/vcf_files
```
This command takes in 4 parameters:
* `--db-path` - path to directory storing the TileDB instances. The
  AGC compressed genomes will be placed here on completion.
* `--dataset-type` - the type of dataset to export. In this case, we
  are exporting the hVCF data.
* `--sample-names` - a comma-separated list of sample names to export.
* `--output-dir` - the directory to place the exported hVCF files.

In our above example, we can use the `--sample-names` parameter, but
for certain scenarios this may become inefficient (e.g., exporting
hundreds of samples). To remedy this, we can bypass the 
`--sample-names` parameter and replace it with the `--sample-file`
parameter. This can take in a file that contains a list of samples
you would like to export from the database. For example, if we would
want to export the sample samples from above but in file form, the
file would look like the following:

```
$ cat sample_names.txt

LineA
LineB
```

When we 

## K-mer Indexing

Once we have the hVCF data exported, we need to index the k-mers. 
This is done using the `build-kmer-index` command.

```shell
./phg build-kmer-index \
    --db-path /my/db/uri \
    --hvcf-dir /my/hvcf/dir 
```

This command has 2 required parameters:
* `--db-path` - path to directory storing the TileDB instances. The
  AGC compressed genomes are also stored here as well under 
  `db-path/assemblies.agc`.
* `--hvcf-dir` - the directory containing the hVCF files. This is the
  output directory from the `export-vcf` command.  Right now this is 
  required, but will be optional in the future.

This will store the k-mer index as `kmerIndex.txt` in the `--hvcf-dir` 
directory.

In addition, this command can take optional parameters:

* `--index-file` - The full path of the k-mer index file. 
  Default = <hvcf-dir>/`kmerIndex.txt`.  
* `--maxHapProportion` or `-p` - only k-mers mapping to less than or 
  equal to maxHapProportion of haplotypes in a reference range will 
  be retained. Default = `0.75`.  
* `--max-arg-length` - The maximum argument length for a call to AGC. 
  Default = `200000`. If you get an error caused by a call to agc 
  being too long try reducing this value."

The following optional parameters affect how k-mers are pre-filtered 
to determine which are used for indexing. They would only need to be 
adjusted if the number of k-mers in the index is too low or too high.

* `--hashMask` or `-m` - with hashFilter, used to mask k-mers for 
  filtering. Default uses only the last k-mer nucleotide. Only change 
  this if you know what you are doing. Default = `3`.
* `--hashFilter` or `-f` - Only hashes that pass the filter 
  ((hashValue and hashMask) == hashFilter) will be considered. Do not 
  change this value unless you know what you are doing. 
  Default = `1`.


## Read Mapping

Now that we have the index we can align short reads against the PHG 
using the `map-kmers` command.

```shell  
./phg map-kmers \
    --hvcf-dir /my/hvcf/dir \
    --kmer-index /my/hvcf/dir/kmerIndex.txt \
    --read-files /path/to/reads/LineA_R1.fq /path/to/reads/LineA_R2.fq \
    --output-dir /my/mapping/dir
```

This command has the following parameters:
* `--hvcf-dir` - the directory containing the hVCF files.
* `--kmer-index` - the k-mer index file created by the 
  `build-kmer-index` command.
* `--read-files` - a comma separated list of fastq files for a single 
  sample.  Either 1(for single end) or 2(for paired end) files can be 
  input at a time this way. Any more and an error will be thrown.
* `--output-dir` - the directory to place the read mapping files.

Instead of using `--read-files`, you can use `--key-file` to specify 
a key file that contains a list of the read mapping files. 
Columns for samplename and filename are required.  If using paired end 
FASTQs, a filename2 column can be included. File names must end in 
one of `.fq`, `.fq.ga`, `.fastq`, or `.fastq.gz`.

The value of `--kmer-index` defaults to `--hvcf-dir`/kmerIndex.txt, 
the default value used by BuildKmerIndex. If a non-default value was 
used for the `kmerIndex` file in BuildKmerIndex, that same value 
needs to be set here using this parameter.



## Find Paths
FindPaths is the command used to impute paths through a 
HaplotypeGraph based on a set of read mappings. The method uses read 
mapping files generated by the map-kmers command. Running FindPaths 
uses the Viterbi algorithm to solve a hidden Markov model (HMM) to 
identify the set of haplotypes most likely to have generated the 
observed read mappings. A haploid path is a single path through the 
graph. A diploid path is two paths through the graph that jointly 
explain the read mappings.


### Input Files
FindPaths can take either read-mapping files generated by the command 
MapKmers or fastq files of sequence reads. If fastq files are used, 
MapKmers is run to generate read mappings but the intermediate 
read-mapping files are not saved. To save a copy of the read-mapping 
files, MapKmers and FindPaths must be run separately. For either 
input type, the key-file must contain columns named sampleName and 
filename. Optionally, key files for fastqs can also contain a 
"filename2" column for paired-end reads. In either case, the key file 
used to specify the files to be processed can take single files or 
comma-separated lists for filename and filename2. Reads from lists of 
files will be combined to impute one path.

Using a key-file to specify inputs for multiple samples will be the 
most flexible. However, in some scenarios specifying a single input 
file on the command line is useful. In that case, the sample name is 
taken from the filename. For fastq files the sample name is filename 
with the .fastq or .fq extension removed. For read mapping files, the 
sample name will be the filename with the _readMapping.txt extension 
removed.

Fastq files must end in one of `.fq`, `.fq.gz`, `.fastq`, or 
`.fastq.gz`. Read mapping files must end in `_readMapping.txt`.


### Likely Ancestors
The term "ancestors" refers to the taxa/assemblies used to construct 
the haplotype graph used for imputation. The term "ancestors" is used 
because the model considers the haplotypes in the graph to represent 
the potential ancestral haplotypes from which the samples to be 
imputed were derived. In some scenarios, limiting the number of 
ancestors used to impute a sample will be useful. For example, if 
some samples were derived from an F1 cross and the individuals used 
to build the graph include those parents, using only those parents 
for imputation can improve accuracy. Also, because imputing diploid 
paths is more computationally intensive than haploid paths, limiting 
the number of ancestors used per sample may be necessary to control 
the amount of time required per sample.

To restrict the number of ancestors used, set use-likely-ancestors to 
true, and provide a value for either max-ancestors or min-coverage. 
For the case where the samples have been derived from a limited 
number of ancestors setting min-coverage to 0.95 or max-ancestors to 
the expected number of ancestors is a useful strategy. In either 
case, providing a name for the output file saves a record of the 
ancestors used for each sample and should be checked to make sure the 
samples behaved as expected.

When using this method, ancestors are chosen by first counting the 
number of reads for a sample that map to each ancestor in the 
haplotype graph. The ancestor with the most mapped reads is chosen. 
Next the reads not mapping to first chosen ancestor are used to count 
reads mapping to the remaining ancestors and the ancestor with the 
most reads is selected. This is repeated until either the proportion
of reads mapping to the selected ancestors >= min-coverage or the 
number ancestors selected equals max-ancestors. If a name is provided 
for the likely-ancestor-file, the chosen ancestors are written to 
that file in the order they were chosen along with the cumulative 
coverage, which is the proportion of reads mapping to each ancestor 
and the previously chosen ones.


### Required Parameters
* `--path-keyfile` **or** `--read-file` **but not both**: 
  + `--path-keyfile` is the name and path of a keyfile that contains 
    a list of the read mapping files or read files (FASTQ). The 
    keyfile must have two columns labeled sampleName and filename. 
    Any additional columns will be ignored. SampleName must be unique.
    Additionally, if an output hvcf for a sample name already exists 
    in the output directory, the sample will be skipped. See the 
    **Input Files** section above for more detail.

* `--hvcf-dir`: The directory containing the hvcf used to build the 
  haplotype graph used for imputation.

* `--reference-genome`: The name and path to the reference genome 
  fasta or fastq file.

* `--output-dir`: The directory where the output hVCFs will be 
  written. One file will be written for each sample, and the output 
  file name will be `(sampleName).h.vcf`.

* `--path-type`: The type of path to be imputed. The value must be 
  either "haploid" or "diploid" (without quotes).


### Optional Parameters
* `--prob-correct`: The probability that a mapped read was mapped 
  correctly. (Default = 0.99)

* `--prob-same-gamete`: The probability of transitioning to the same 
  gamete (sample) in the next reference range. This should be equal 
  to 1 - (probability of a recombination). Probability of a 
  recombination can be estimated as the total number of expected 
  recombinations in a sample divided by the number of reference 
  ranges (Default = `0.99`)

* `--min-gametes`: The minimum number of gametes with a haplotype in a reference range. Reference ranges with fewer gametes will not be imputed. (Default = 1)

* `--min-reads`: The minimum number of reads per ReferenceRange. Reference ranges with fewer reads will not be imputed.
If minReads = 0, all ReferenceRanges will be imputed. (Default = 0)

* `--inbreeding-coefficient`: The estimated coefficient of inbreeding for the samples being evaluated. Only used for diploid path type. The value must be between 0.0 and 1.0 (Default = 0.0)

* `--max-reads-per-kb`: ReferenceRanges with more than max-reads-per-kb will not be imputed. (Default = 1000)

* `--use-likely-ancestors`: The value must be "true" or "false" (no parentheses). This indicates whether the most
likely ancestors of each sample will be used for path finding. (Default = false)

* `--max-ancestors`: If use-likely-ancestors = true, use at most max-ancestors. (Default = Integer.MAX_VALUE)

* `--min-coverage`: If use-likely-ancestors = true, use the fewest number of ancestors that together have this proportion
of mappable reads. The values must be between 0.5 and 1.0 (Default = 1.0)

* `--likely-ancestor-file`: If useLikelyAncestors is true, a record of the ancestors used for each sample will be written
to this file, if a filename is provided. (Default = "")

* `--threads`: number of threads used to find paths. (Default = 3)


The command
````
phg find-paths --help
````  
will print the following:

```shell
Usage: find-paths [<options>]

Impute best path(s) using read mappings.

Options:  
--path-keyfile=\<value\>           Name of tab-delimited key file. Columns for samplename and filename are required. Files must be either read mapping files (ending in_readMapping.txt) or fastq files. If using paired end fastqs, a filename2 column can be included. A value must be entered for either --key-file or --read-files.  
--read-files=\<value\>             Comma separated list of fastq files for a single sample. Either 1(for single end)or 2(for paired end) files can be input at a time this way. Any more and an error will be thrown. If listing files from read mapping, they must end in _readMapping.txt.  
--hvcf-dir=\<text\>                The directory containing the hvcf files used to build a HaplotypeGraph for path finding. Required parameter.  
--reference-genome=\<text\>        path to reference genome (fasta or fastq). Required parameter.    
--output-dir=\<text\>              The directory where the output hvcfs will be written. The output file names will be<sampleName>.h.vcf. Required parameter.    
--path-type=(haploid|diploid)      The type of path to find. Must be lower case'haploid' or 'diploid' (without quotes).'haploid' infers a single path through the graph. 'diploid' infers a pair of paths. Required parameter.  
--kmer-index=\<text\>              The name and path to the kmerIndex. Default =<hvcfDir>/kmerIndex.txt  
--prob-correct=\<floa\>            The probability that a mapped read was mapped correctly. Default = 0.99  
--prob-same-gamete=\<float\>       The probability of transitioning to the same gamete (sample) in the next reference range. Default = 0.99  
--min-gametes=\<int\>              The minimum number of gametes with a haplotype in a reference range. Reference ranges with fewer gametes will not be imputed. Default = 1  
--min-reads=\<int\>                The minimum number of reads per ReferenceRange. Reference ranges with fewer reads will not be imputed. If minReads = 0, all ReferenceRanges will be imputed. Default = 0  
--inbreeding-coefficient=\<float\> The estimated coefficient of inbreeding for the samples being evaluated. Only used for diploid path type. Default = 0.0  
--max-reads-per-kb=\<int\>         ReferenceRanges with more than max-reads-per-kb will not be imputed. Default = 1000.  
--use-likely-ancestors=true|false  Use only the most likely ancestors of each sample for path finding. Default = false  
--max-ancestors=\<int\>            If use-likely-ancestors = true, use at most max-ancestors. Default = Int.MAX_VALUE.  
--min-coverage=\<float\>           If use-likely-ancestors = true, use the fewest number of ancestors that together have this proportion of mappable reads. Default = 1.0  
--likely-ancestor-file=\<text\>    If useLikelyAncestors is true, a record of the ancestors used for each sample will be written to this file. Default = ''  
--threads=\<int\>                  number of threads used to find paths. Default = 3.  
-h, --help                         Show this message and exit    
```