# hVCF - Haplotype Variant Call Format Specification

* **_Specification version:_** `v2.2`
* **_Date:_** 2023-12-07

## Overview
hVCF stands for **h**aplotype **V**ariant **C**all **F**ormat. This 
format is used to store and encode haplotype information across 
samples from the PHG ([Practical Haplotype Graph](https://github.com/maize-genetics/phg_v2)).
An hVCF file is based on the standards of a VCF file, specifically, 
[VCF v4.2](http://samtools.github.io/hts-specs/VCFv4.2.pdf). This 
format leverages VCF's symbolic allele information from the `ALT` 
field.

## The hVCF specification
hVCF files can be broken into 3 main components:
* Meta-information lines
* Header line
* Data lines containing information for each reference range
  * Fixed fields
  * Haplotype fields

### An example
The following code block illustrates a formatted example hVCF file:

```
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##ALT=<ID=06ae4e937668d301e325d43725a38c3f,Description="haplotype data for line: Ref",Source="data/test/smallseq/Ref.fa",Regions=1:45001-49500,Checksum=Md5,RefRange=06ae4e937668d301e325d43725a38c3f>
##ALT=<ID=073286a82fe47d6a370e8a7a3803f1d3,Description="haplotype data for line: Ref",Source="data/test/smallseq/Ref.fa",Regions=1:39501-44000,Checksum=Md5,RefRange=073286a82fe47d6a370e8a7a3803f1d3>
##ALT=<ID=105c85412229b45439db1f03c3f064f4,Description="haplotype data for line: Ref",Source="data/test/smallseq/Ref.fa",Regions=1:27501-28500,Checksum=Md5,RefRange=105c85412229b45439db1f03c3f064f4>
##ALT=<ID=105e63346a01d88e8339eddf9131c435,Description="haplotype data for line: Ref",Source="data/test/smallseq/Ref.fa",Regions=2:50501-55000,Checksum=Md5,RefRange=105e63346a01d88e8339eddf9131c435>
##ALT=<ID=2c4b8564bbbdf70c6560fdefdbe3ef6a,Description="haplotype data for line: Ref",Source="data/test/smallseq/Ref.fa",Regions=2:34001-38500,Checksum=Md5,RefRange=2c4b8564bbbdf70c6560fdefdbe3ef6a>
##ALT=<ID=347f0478b1a553ef107243cb60a9ba7d,Description="haplotype data for line: Ref",Source="data/test/smallseq/Ref.fa",Regions=2:11001-12000,Checksum=Md5,RefRange=347f0478b1a553ef107243cb60a9ba7d>
##ALT=<ID=39f96726321b329964435865b3694fd2,Description="haplotype data for line: Ref",Source="data/test/smallseq/Ref.fa",Regions=2:49501-50500,Checksum=Md5,RefRange=39f96726321b329964435865b3694fd2>
##ALT=<ID=43687e13112bbe841f811b0a9de82a94,Description="haplotype data for line: Ref",Source="data/test/smallseq/Ref.fa",Regions=2:22001-23000,Checksum=Md5,RefRange=43687e13112bbe841f811b0a9de82a94>
##ALT=<ID=546d1839623a5b0ea98bbff9a8a320e2,Description="haplotype data for line: Ref",Source="data/test/smallseq/Ref.fa",Regions=1:1-1000,Checksum=Md5,RefRange=546d1839623a5b0ea98bbff9a8a320e2>
##ALT=<ID=57705b1e2541c7634ea59a48fc52026f,Description="haplotype data for line: Ref",Source="data/test/smallseq/Ref.fa",Regions=1:1001-5500,Checksum=Md5,RefRange=57705b1e2541c7634ea59a48fc52026f>
##ALT=<ID=1bda8c63ae8e2f3678b85bac0ee7b8b9,Description="haplotype data for line: B97",Source="data/test/smallseq/B97.fa",Regions=1:1250-6750,Checksum=Md5,RefRange=57705b1e2541c7634ea59a48fc52026f>
##ALT=<ID=5fedf293a1a5443cc896d59f12d1b92f,Description="haplotype data for line: CML231"Source="data/test/smallseq/CML231.fa",Regions=2:22001-23000,Checksum=Md5,RefRange=43687e13112bbe841f811b0a9de82a94>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##INFO=<ID=END,Number=1,Type=Integer,Description="Stop position of the interval">
##contig=<ID=1,length=55000>
##contig=<ID=2,length=55000>
##reference=https://s3.amazonaws.com/maizegenetics/phg/phgV2Test/Ref.fa
#CHROM POS   ID REF ALT                                                                   QUAL FILTER INFO      FORMAT Ref B97 CML231
1      1     .  G   <546d1839623a5b0ea98bbff9a8a320e2>                                    .    .      END=1000  GT     1   1   1
1      1001  .  A   <57705b1e2541c7634ea59a48fc52026f>,<1bda8c63ae8e2f3678b85bac0ee7b8b9> .    .      END=5500  GT     1   2   1
1      27501 .  G   <105c85412229b45439db1f03c3f064f4>                                    .    .      END=28500 GT     1   1   1
1      39501 .  G   <073286a82fe47d6a370e8a7a3803f1d3>                                    .    .      END=44000 GT     1   1   1
1      45001 .  G   <06ae4e937668d301e325d43725a38c3f>                                    .    .      END=49500 GT     1   1   1
2      11001 .  A   <347f0478b1a553ef107243cb60a9ba7d>                                    .    .      END=12000 GT     1   1   1
2      22001 .  A   <43687e13112bbe841f811b0a9de82a94>,<5fedf293a1a5443cc896d59f12d1b92f> .    .      END=23000 GT     1   1   2
2      34001 .  A   <2c4b8564bbbdf70c6560fdefdbe3ef6a>                                    .    .      END=38500 GT     1   1   1
2      49501 .  G   <39f96726321b329964435865b3694fd2>                                    .    .      END=50500 GT     1   1   1
2      50501 .  G   <105e63346a01d88e8339eddf9131c435>                                    .    .      END=55000 GT     1   1   1
```
!!! note
    In the prior example, the hVCF output columns below the header line
    (e.g. below the line starting with `#CHROM`) are formatted for 
    visual clarity. In a real example, delimiters are tab (`\t`) based.

### Meta-information lines
The header portion of an hVCF file contain rows of "meta-information"
, which are lines that start with `##` and must appear first in the 
file. Like a VCF file, hVCF files can contain both **unstructured** 
and **structured** meta-information.

**Unstructured** meta-information is characterized by a 
straightforward pairing of `key=value` logic. For instance, in the 
previous illustration, the `##fileformat=VCF4.4` represents 
unstructured meta-information, with `##=fileformat` serving as the 
_key_ and `VCF4.4` as the corresponding _value_.

**Structured** meta-information also consists of a key-value pair,
but in this case, the value is a collection of additional key-value
pairs separated by a comma (`,`) and enclosed with angle brackets
(`<` and `>`). In our prior example, sequence information fields
(e.g. `##contig=<ID=chr7,length=461>`) represent structured 
meta-information where `##contig` is the primary key and 
`<ID=chr7,length=461>` is the value containing nested key-value pairs:

* `ID=chr7`
* `length=461`

#### File format (`##fileformat`) field
Similar to VCF, a single line containing file format information
(e.g. `fileformat`) _must_ be the _first line in the file_. In the
case of hVCF files, the version must be version 4.4 of the VCF
specification (`VCFv4.4`).


#### Alternative allele (`##ALT`) field
The primary driver of the hVCF specification is information stored
within the structured alternative allele field. At its core, the 
alternative allele field contains two primary key-values pairs, the 
`ID` and `Description` which describe symbolic alternate alleles in 
the `ALT` column of VCF records. While this field is usually used 
to describe possible structural variants and 
[IUPAC ambiguity codes](https://www.bioinformatics.org/sms/iupac.html),
here it is used to represent a haplotype sequence and ID for a 
given reference range. This is achieved by defining the `ID` value
with an [MD5](https://en.wikipedia.org/wiki/MD5) checksum of the
given haplotype sequence and defining the `Description` value with
information about the origin of the haplotype sequence. 

Since this haplotype sequence is (I) derived from a particular 
sample, (II) related to reference range information, 
and (III) has its own positional information, we can populate 
the alternative allele field with additional key-value information. 
Take the following example:

```
##ALT=<ID=06ae4e937668d301e325d43725a38c3f,Description="haplotype data for line: Ref",Source="data/test/smallseq/Ref.fa",Regions=1:45001-49500,Checksum=Md5,RefRange=06ae4e937668d301e325d43725a38c3f>
```

Here, we have the following information:

| Key           | Value                              | Description                                                                                                                                   |
|---------------|------------------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------|
| `ID`          | `06ae4e937668d301e325d43725a38c3f` | MD5 checksum identifier for haplotype sequence                                                                                                |
| `Description` | `"haplotype data for line: Ref"`   | Information about the origin of the haplotype sequence                                                                                        |
| `Source`      | `"data/test/smallseq/Ref.fa"`      | Fasta file ID and path containing haplotype sequence                                                                                          |
| `Regions`     | `1:45001-49500`                    | List of regions which make up the haplotype.  They are in the format contig:startPos-endPos.  All the sub-regions will be separated by commas |
| `Checksum`    | `Md5`                              | Identifier describing the digest hashing algorithm for the `ID` key                                                                           |
| `RefRange`    | `20b13a1d5b466ff438174aa897c985f7` | MD5 checksum identifier for reference sequence                                                                                                |

#### Individual format (`##FORMAT`) field
The meta-information contained in the individual format field closely 
adheres to the VCF specification. This structured field provides a 
description of the IDs found within the `FORMAT` column of the data 
rows. The necessary keys are as follows:

| Key           | Description                                                                          |
|---------------|--------------------------------------------------------------------------------------|
| `ID`          | Identifier for `FORMAT` entry                                                        |
| `Number`      | Number (integer) of values representing `ID`                                         |
| `Type`        | [Data type](http://samtools.github.io/hts-specs/VCFv4.4.pdf#subsection.1.3) for `ID` |
| `Description` | Descriptive information about `ID`                                                   |


#### Information (`##INFO`) field 
Much like the `##FORMAT` field, the `##INFO` field is a structured 
meta-information field that provides details pertaining to each 
reference range and the corresponding haplotype data contained within 
those reference ranges. Similar to the `##FORMAT` field, the 
necessary keys are as follows:

| Key           | Description                                                                          |
|---------------|--------------------------------------------------------------------------------------|
| `ID`          | Identifier for `FORMAT` entry. Defaults to `END`                                     |
| `Number`      | Number (integer) of values representing `ID`                                         |
| `Type`        | [Data type](http://samtools.github.io/hts-specs/VCFv4.2.pdf#subsection.1.3) for `ID` |
| `Description` | Descriptive information about `ID`                                                   |

In order to properly represent the information regarding reference
ranges, the following values are required:

| Value        | Description                          |
|--------------|--------------------------------------|
| `End`        | End position of reference range (bp) |

By combining the values identified within the `POS` column and the
`END` value, we can specify the total length of the reference range
along with assembly information.


#### Sequence information (`##contig`) field
The contig field is used to detail additional attributes for
each sequence represented within the haplotype data. For now,
this is a structured field requiring the identifier (`ID`) for the 
sequence and the length of the mentioned sequence (`length`)


### Header
Like the VCF specifications, the 8 mandatory tab-delimited (`\t`) 
column headers are required:

* `#CHROM`
* `POS`
* `ID`
* `REF`
* `ALT`
* `QUAL`
* `FILTER`

Since genotype (i.e. haplotype information) data is also required,
the `FORMAT` column is also required with the `GT` identifier.
More information about each of these fields is discussed in the
next section.

!!! note
    The end of the line must have no tab characters (`\t`).


### Data lines

#### Fixed fields
There are 8 fixed fields for each reference range record:

| Field    | Description                                                                                       |
|----------|---------------------------------------------------------------------------------------------------|
| `CHROM`  | sequence identifier                                                                               |
| `POS`    | start position of reference range                                                                 |
| `ID`     | an optional identifier for a reference range record                                               | 
| `REF`    | The allele at the start position for reference haplotype sequence                                 |
| `ALT`    | MD5 hash sums for each possible haplotype sequence detected at a given reference range record     |
| `QUAL`   | quality score (_needed to satisfy VCF specifications_)                                            |
| `FILTER` | filter status (_needed to satisfy VCF specifications_)                                            |
| `INFO`   | information field used to represent the end position (`END`) value for the reference range record |

#### Haplotype fields
Haplotype information must be specified by first creating a format
field (`FORMAT`) field along with a genotype (`GT`) identifier.

The following fields proceeding the `FORMAT` field are specified with 
the sample (e.g. taxa) identifier for each given sample referenced
in the hVCF file. For example, let's take a look at a given record
in the prior example (_with added header for additional clarity_):

```
#CHROM POS   ID REF ALT                                                                   QUAL FILTER INFO      FORMAT Ref B97 CML231
1      1001  .  A   <57705b1e2541c7634ea59a48fc52026f>,<1bda8c63ae8e2f3678b85bac0ee7b8b9> .    .      END=5500  GT     1   2   1
```

One thing you will notice is that there are no calls to the 
"reference" allele field; only calls to the alternate field
since these allele values represent the haplotype sequence in
MD5 hash form. Allele values, **if using haploid path finding**,
are represented using 
[singular values](https://samtools.github.io/hts-specs/VCFv4.2.pdf#subsubsection.1.4.2) 
(e.g. `1`, `2`) which represent the indexed order of haplotype 
sequences in the `ALT` field. In other terms, if a sample has an 
allele value of `1`, this would refer to the _first_ symbolic allele 
in the `ALT` field for the haploid value.

Using this information with the prior example, we can infer the 
following haplotype sequence information for the given reference 
range record (`1:1001-5500`):

| Sample ID | Allele values | MD5 symbolic allele                |
|-----------|---------------|------------------------------------|
| `Ref`     | `1`           | `57705b1e2541c7634ea59a48fc52026f` |
| `B97`     | `2`           | `1bda8c63ae8e2f3678b85bac0ee7b8b9` |
| `CML231`  | `1`           | `57705b1e2541c7634ea59a48fc52026f` |

Alternatively, allele values in hVCF files can be generated **using 
diploid path finding** during the PHGv2 imputation process. Here is
an example entry of this:

```
#CHROM POS   ID REF ALT                                                                   QUAL FILTER INFO      FORMAT Ref B97 CML231
1      1001  .  A   <57705b1e2541c7634ea59a48fc52026f>,<1bda8c63ae8e2f3678b85bac0ee7b8b9> .    .      END=5500  GT     1|1 2|1 1|1
```

Allele values are separated with a "phased" indicator (`|`) and
never with an "unphased" indicator (`/`). Similar to haploid path
finding, allele values represent the indexed order of haplotype
sequences in the `ALT` field. In other terms, if a sample has an 
allele value of `2|1`, this would refer to the _second_ symbolic 
allele in the `ALT` field for the first gamete and the _first_ 
symbolic allele for the second gamete.

Using this information with the prior example, we can infer the
following haploid sequence information for the given reference
range record (`1:1001-5500`) using diploid values:

| Sample ID | Allele values | MD5 symbolic allele (gamete 1)     | MD5 symbolic allele (gamete 2)     |
|-----------|---------------|------------------------------------|------------------------------------|
| `Ref`     | `1\|1`        | `57705b1e2541c7634ea59a48fc52026f` | `57705b1e2541c7634ea59a48fc52026f` |
| `B97`     | `2\|1`        | `1bda8c63ae8e2f3678b85bac0ee7b8b9` | `57705b1e2541c7634ea59a48fc52026f` |
| `CML231`  | `1\|1`        | `57705b1e2541c7634ea59a48fc52026f` | `57705b1e2541c7634ea59a48fc52026f` |
