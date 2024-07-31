# Exporting Data

In this document, we will discuss general strategies for exporting
data from a PHG database.

!!! note
    This section will assume you have a pre-existing database loaded
    with haplotype data. See the ["Building and Loading"](build_and_load.md)
    documentation for further information.

## Quickstart

* Export hVCF files from database
```shell
phg export-vcf \
    --db-path /path/to/dbs \
    --dataset-type hvcf \ # can also be 'gvcf'
    --sample-names LineA,LineB \ # comma separated list of sample names
    -o /path/to/output/directory
```
!!! note
    `--sample-names` can be replaced with `--sample-file <file_name.txt>` 
    where `file_name.txt` is a file containing the sample names, one 
    per line.

* Create FASTA files from hVCF data or database
```shell
phg create-fasta-from-hvcf \
  --hvcf-dir my/hvcf_dir \ # can also be an individual file ('--hvcf-file')
  --fasta-type composite \ # can also be 'haplotype'
  -o /path/to/output_folder
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
* `--regions-file` - a file of positions to be exported. Can be a 
   BED file or a VCF file.

!!! note
    Make sure there is no whitespace between sample IDs. For example:

    * `LineA,LineB` ✅
    * `LineA , LineB` ❌

Users may instead use the `--sample-file` parameter to specify a file 
that contains the sample names, one per line. For example, if I have
a text file called `sample_names.txt`, the contents of the file would
look like the following:

```
LineA
LineB
```

...and would be passed to the `export-vcf` command:

```shell
phg export-vcf \
    --db-path vcf_dbs \
    --dataset-type hvcf \
    --sample-file sample_names.txt \
    -o output/hvcf_files
```
If input is specified for the `--regions-file` parameter, only 
variants overlapping those positions will be exported. 
The regions-file must be either a 
[BED file](https://en.wikipedia.org/wiki/BED_(file_format)) or a 
[VCF file](https://en.wikipedia.org/wiki/Variant_Call_Format), and 
must have either a `.bed` or a `.vcf` extension.

For example, if I want the regions from `1` to `5000` base pairs (bp) on
chromosome 3 (in my case the ID would be `chr03`), I could make a BED 
file:

```
chr03 0 5000
```

!!! note
    BED files are 0-based, so plan accordingly!

...or this could be a VCF file that contains a data line for `chr03`
region information for the `CHROM`, `POS`, and
`INFO` columns with the `INFO` column containing a `END` field. For
example:

```
#CHROM  POS ID  REF ALT QUAL  FILTER  INFO  ...
chr03   1   r1  A   T   50    PASS    END=5000
```


### Create FASTA data
While haplotype sequences are abstracted to MD5 hashes in hVCF
files, sequence information can be recapitulated from these hash 
values using the `create-fasta-from-hvcf` command:

```shell
phg create-fasta-from-hvcf \
  --hvcf-file my_sample.h.vcf \
  --fasta-type composite \
  -o /path/to/outputFolder
```

As the name of this command implies, we are creating FASTA files
of nucleotide sequence data from a single hVCF file or a collection
of hVCF files by specifying a directory. The output FASTA files will
be written (one FASTA per hVCF file) to the specified output directory
(`-o`). The format of the file names will be `sample_name_type.fa` 
where `sample_name` is the name of the sample from the hVCF file name 
and `type` is the type of fasta file created (`composite` or 
`haplotype`). The following parameters may be used:

* input type (**you can only select one**):
    + `--hvcf-file` - path to an hVCF file. **Can be substituted with
       `--hvcf-dir`**.
    + `--hvcf-dir` - path to a directory containing hVCF files. **Can
      be substituted with `--hvcf-file`**.
* `--fasta-type` - what type of FASTA format do you want to use?
    + `composite` - generate a FASTA file that contains all haplotypes 
      concatenated together by consecutive reference ranges. This 
      composite or "pseudo" genome can be **used for the resequencing 
      pipeline**.
    + `haplotype` - generate a FASTA file where each haplotype is a
      seperate FASTA entry. **Useful for read mapping, imputation
      or simple haplotype sequence retrieval**.
* `-o` - output path to directory for the created fasta files.


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
  ```shell
  # An example pointing to a composite hVCF file
  $ curl http://localhost:8080/brapi/v2/variantsets
  ```

* using the R package, `rPHG2`. Since this is a separate library,
  more information about the library and retrieval methods can be
  found [here](https://rphg2.maizegenetics.net/articles/rPHG2.html).
