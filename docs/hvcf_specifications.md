# hVCF - Haplotype Variant Call Format Specification

* **_Specification version:_** `v2.1`
* **_Date:_** 2023-10-31

## Overview
hVCF stands **h**aplotype **V**ariant **C**all **F**ormat. This 
format is used to store and encode haplotype information across 
samples from the PHG ([Practical Haplotype Graph](https://github.com/maize-genetics/phg_v2)).
An hVCF file is based on the standards of a VCF file, specifically, 
[VCF v4.4](http://samtools.github.io/hts-specs/VCFv4.4.pdf). This 
format leverages VCF's symbolic allele information from the `ALT` 
field.

## The hVCF specification
hVCF files can be broken into 4 main components:
* Meta-information lines
* Header line
* Data rows containing information for each reference range
* Haplotype information on samples for each reference range

### An example
The following code block illustrates an example hVCF file:

```
##fileformat=VCFv4.4
##ALT=<ID=4c81af3c2f102652b58a3af21355bc25,Description="haplotype data for line: B97",Number=9,Source="data/test/buildMAFVCF/B97_ASM_Test.fa",Contig=chr7,Start=451,End=456,Asm_Contig=chr4,Asm_Start=5247,Asm_End=5252,Checksum=Md5,RefRange=20b13a1d5b466ff438174aa897c985f7>
##ALT=<ID=b935ff8898a7942d09bfc11e73e8a555,Description="haplotype data for line: B97",Number=9,Source="data/test/buildMAFVCF/B97_ASM_Test.fa",Contig=chr1,Start=1,End=40,Asm_Contig=chr6,Asm_Start=98,Asm_End=142,Checksum=Md5,RefRange=e2cb401a9bd8a1bcca158761b755a6a8>
##ALT=<ID=c864de7dcaa52600fe611032f79d20e3,Description="haplotype data for line: B97",Number=9,Source="data/test/buildMAFVCF/B97_ASM_Test.fa",Contig=chr7,Start=15,End=48,Asm_Contig=chr4,Asm_Start=4245,Asm_End=4281,Checksum=Md5,RefRange=2af5a69357a7eb4a537fa3d08dea18de>
##FORMAT=<ID=AD,Number=3,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth (only filtered reads used for calling)">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">
##INFO=<ID=AF,Number=3,Type=Integer,Description="Allele Frequency">
##INFO=<ID=ASM_Chr,Number=1,Type=String,Description="Assembly chromosome">
##INFO=<ID=ASM_End,Number=1,Type=Integer,Description="Assembly end position">
##INFO=<ID=ASM_Start,Number=1,Type=Integer,Description="Assembly start position">
##INFO=<ID=ASM_Strand,Number=1,Type=String,Description="Assembly strand">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=END,Number=1,Type=Integer,Description="Stop position of the interval">
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
##contig=<ID=chr7,length=461>
##contig=<ID=chr1,length=40>
##contig=<ID=chr10,length=40>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	B97
chr1	1	.	G	<b935ff8898a7942d09bfc11e73e8a555>	.	.	END=40	GT:AD:DP	1:0,2,0:2
chr10	1	.	G	<b935ff8898a7942d09bfc11e73e8a555>	.	.	END=40	GT:AD:DP	1:0,2,0:2
chr7	15	.	A	<c864de7dcaa52600fe611032f79d20e3>	.	.	END=48	GT:AD:DP	1:0,2,0:2
chr7	451	.	T	<4c81af3c2f102652b58a3af21355bc25>	.	.	END=456	GT:AD:DP	1:0,2,0:2
```

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
sample, (II) related to related to reference range information, 
and (III) has its own positional information, we can populate 
the alternative allele field with additional key-value information. 
Take the following example:

```
##ALT=<ID=4c81af3c2f102652b58a3af21355bc25,Description="haplotype data for line: B97",Number=9,Source="data/test/buildMAFVCF/B97_ASM_Test.fa",Contig=chr7,Start=451,End=456,Asm_Contig=chr4,Asm_Start=5247,Asm_End=5252,Checksum=Md5,RefRange=20b13a1d5b466ff438174aa897c985f7>
```

Here, we have the following information:

| Key           | Value                                     | Description                                                         |
|---------------|-------------------------------------------|---------------------------------------------------------------------|
| `ID`          | `4c81af3c2f102652b58a3af21355bc25`        | MD5 checksum identifier for haplotype sequence                      |
| `Description` | `"haplotype data for line: B97"`          | Information about the origin of the haplotype sequence              |
| `Number`      | `9`                                       | Number of fields following the `Description` key                    |
| `Source`      | `"data/test/buildMAFVCF/B97_ASM_Test.fa"` | Fasta file ID and path containing haplotype sequence                |
| `Contig`      | `chr7`                                    | Sequence identifier for reference range                             |
| `Start`       | `451`                                     | Start position (bp) for reference range                             |
| `End`         | `456`                                     | End position (bp) of reference range                                |
| `ASM_Contig`  | `chr4`                                    | Sequence identifier for assembly data containing haplotype sequence |
| `ASM_Start`   | `5247`                                    | Start position (bp) for assembly data containing haplotype sequence |
| `ASM_End`     | `5252`                                    | End position (bp) for assembly data containing haplotype sequence   |
| `Checksum`    | `Md5`                                     | Identifier describing the digest hashing algorithm for the `ID` key |
| `RefRange`    | `20b13a1d5b466ff438174aa897c985f7`        | MD5 checksum identifier for reference sequence                      |

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
those reference ranges. The necessary keys are as follows:

| Key           | Description                                                                          |
|---------------|--------------------------------------------------------------------------------------|
| `ID`          | Identifier for `FORMAT` entry                                                        |
| `Number`      | Number (integer) of values representing `ID`                                         |
| `Type`        | [Data type](http://samtools.github.io/hts-specs/VCFv4.4.pdf#subsection.1.3) for `ID` |
| `Description` | Descriptive information about `ID`                                                   |

In order to properly represent the information regarding reference
ranges, the following values are required:

| Value        | Description                          |
|--------------|--------------------------------------|
| `AF`         | Allele frequency                     |
| `ASM_Chr`    | Assembly chromosome                  |
| `ASM_Start`  | Assembly start position (bp)         |
| `ASM_End`    | Assembly end position (bp)           |
| `ASM_Strand` | Assembly strand                      |
| `DP`         | Total depth                          |
| `End`        | End position of reference range (bp) |
| `NS`         | Number of samples with data          |

By combining the values identified within the `POS` column and the
`End` value, we can specify the total length of the reference range
along with assembly information.


#### Sequence information (`##contig`) field
The contig field is used to detail additional attributes for
each sequence represented within the haplotype data. For now,
this is a structured field requiring the identifier (`ID`) for the 
sequence and the length of the mentioned sequence (`length`)


### Header
Like the VCF specifications, the 8 mandatory tab-delimited column 
headers are required:

* `#CHROM`
* `POS`
* `ID`
* `REF`
* `ALT`
* `QUAL`
* `FILTER`

Since genotype (i.e. haplotype information) data is also required,
the `FORMAT` column is also required.

> [!NOTE]
> The end of the line must have no tab characters (`\t`).


### Data rows
_WIP_

### Haplotype information
_WIP_

