# PHGv2 - Exporting Data

In this document, we will discuss general strategies for exporting
data from a PHG database.

> [!NOTE]
> This section will assume you have a pre-existing database loaded
> with haplotype data. See the ["Building and Loading"](build_and_load.md)
> documentation for further information.

## Quickstart

* Export hVCF files from database
* NOTE: --sample-names can be replaced with --sample-file <fileName.txt> where fileName.txt is a file containing the sample names, one per line

```shell
phg export-vcf \
  --db-path /path/to/dbs \
  --dataset-type hvcf \ # can also be 'gvcf'
  --sample-names LineA,LineB \ # comma separated list of sample names
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

## Detailed walkthrough

### Export VCF data
In PHGv2, we leverage [TileDB](https://tiledb.com/) and
[TileDB-VCF](https://docs.tiledb.com/main/integrations-and-extensions/genomics/population-genomics)
for efficient storage and querying of VCF data. For this example,
let's assume I have a pre-existing PHG database (located in my
`vcf_dbs` directory) that contains several 
samples:

* LineA
* LineB
* LineC

If I want to export [hVCF](hvcf_specifications.md) files for a given set of samples 
(`LineA` and `LineC`) to an output directory (in my case 
`output/hvcf_files`), I can use the `export-vcf` command:

```shell
phg export-vcf \
  --db-path vcf_dbs \
  --dataset-type hvcf \
  --sample-names LineA,LineC \
  -o output/hvcf_files
```

This command uses several parameters:

* `--db-path` - path to directory storing the TileDB instances.
* `--dataset-type` - what type of data do you want to extract? This
  can either be:
  + hVCF (`hvcf`) data (**_default parameter_**)
  + gVCF (`gvcf`) data
* `--sample-names` - a comma (`,`) separated list of sample IDs.
* `-o` - output directory of VCF data.

> [!NOTE]
> Make sure there is no whitespace between sample IDs. For example:
> * `SampleA,SampleB` ✅
> * `SampleA , SampleB` ❌

Users may instead use the `--sample-file` parameter to specify a file that contains the sample names, one per line.

```shell
phg export-vcf \
  --db-path vcf_dbs \
  --dataset-type hvcf \
  --sample-file sample_names.txt \
  -o output/hvcf_files
```

### Create FASTA data
While haplotype sequences are abstracted to MD5 hashes in hVCF
files, sequence information can be recapitulated from these hash 
values using the `create-fasta-from-hvcf` command:

```shell
phg create-fasta-from-hvcf \
  --hvcf-file my_sample.h.vcf \
  --fasta-type composite \
  -o my_sequence.fa
```

As the name of this command implies, we are creating FASTA files
of nucleotide sequence data from a single hVCF file or a collection
of hVCF files by specifying a directory:

* `--hvcf-file` - path to an hVCF file. **Can be substituted with
  `--hvcf-dir`**.
* `--hvcf-dir` - path to a directory containing hVCF files. **Can be
  substituted with `--hvcf-file`**.
* `--fasta-type` - what type of FASTA format do you want to use?
  + `composite` - generate a FASTA file that contains all haplotypes 
    concatenated together by consecutive reference ranges. This 
    composite or "pseudo" genome can be **used for rare allele 
    discovery**.
  + `haplotype` - generate a FASTA file where each haplotype is a
    seperate FASTA entry. **Useful for read mapping, imputation
    or simple haplotype sequence retrieval**.
* `-o` - output path for generated FASTA file.


### Data retrieval using BrAPI endpoints and rPHG

While the above commands allow for individual-oriented access to PHG
data, another option is to start a "[REST](https://en.wikipedia.org/wiki/REST)ful"
web service. This service can provide access to a centralized PHG 
database, allowing multiple individuals in a team to simultaneously
retrieve PHG-relevant information. The following web service 
leverages the Breeding API ([BrAPI](https://brapi.org/]\)) which
provides a standard, community-driven collection of web retrieval
calls relevant to plant breeding.

To create a web service for serving PHG data, we can use the
`start-server` command:

```shell
phg start-server \
  --db-path vcf_dbs \
  --port 8080
```

This command takes only two arguments:

* `--db-path` - path to directory storing the TileDB instances.
* `--port` - [web server port](https://en.wikipedia.org/wiki/Port_(computer_networking)) 
  for the network connection. Defaults to `8080`.

Once this command is run, a web service to `localhost` will start
and data can be retrieved:
* manually using 
[BrAPI](https://brapi.org/specification) endpoints and 
[cURL](https://en.wikipedia.org/wiki/CURL):

  ```shell
  # An example pointing to the 'samples' BrAPI endpoint
  $ curl http://localhost:8080/brapi/v2/samples
  ```

* using the R package, `rPHG2`:

  ```r
  # Retrieving same data as prior cURL example
  library(rPHG2)
  
  PHGServerCon("localhost", 8080) |> readSamples()
  ```
  + Since `rPHG2` is its own library, more information about using
    it can be found [here](https://maize-genetics.github.io/rPHG2/).
