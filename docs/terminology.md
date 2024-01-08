# PHG v2 Terminology

## General terms
| Term                | Definition                                                                                          |
|---------------------|-----------------------------------------------------------------------------------------------------|
| Reference genome    | genome used for initial alignment and base coordinates                                              |
| Reference range     | segment of the reference genome                                                                     |
| Haplotype           | sequence of part of an individual chromosome with its start and stop defined by the reference range |
| Reference Haplotype | haplotype from the reference genome                                                                 |
| Alternate genome    | high quality genomes used to identify alternate haplotypes                                          |
| Alternate haplotype | haplotype derived from a genome assembly                                                            |
| Composite genome    | inferred genome based on its composite set of alternate and reference haplotypes                    |
| Haplotype ID        | MD5 checksum for the haplotype sequence                                                             |
| Sample              | genotype (haploid or diploid or higher), taxon, individual                                          |
| Path                | phased set of haplotype ids through the pangenome graph                                             |

## File types
| File Type | Acronym definition                       | Usage                                                                                                             |
|-----------|------------------------------------------|-------------------------------------------------------------------------------------------------------------------|
| `.agc`    | **A**ssembled **G**enomes **C**ompressor | Efficient genome sequence compression.                                                                            |
| `.bam`    | **B**inary **A**lignment **M**ap         | Binary representation of a SAM file; useful for efficient processing.                                             |
| `.bed`    | **B**rowser **E**xtensible **D**ata      | Genomic feature coordinate (e.g. reference ranges) storage.                                                       |
| `.bcf`    | **B**inary **C**all **F**ormat           | Binary representation of a VCF file; useful for efficient processing.                                             |
| `.fasta`  | **FAST**-**A**ll                         | Sequence representation and storage.                                                                              |
| `.g.VCF`  | **g**enomic **VCF** file                 | Variant and non-variant genomic storage.                                                                          |
| `.h.VCF`  | **h**aplotype **VCF** file               | Haplotype information representation and storage. More information can be found [here](./hvcf_specifications.md). |
| `.maf`    | **M**ultiple **A**lignment **F**ormat    | Multiple alignment storage; basis for gVCF and hVCF creation.                                                     |
| `.sam`    | **S**equence **A**lignment **M**ap       | Sequence alignment to a reference sequence.                                                                       |
| `.vcf`    | **V**ariant **C**all **F**ormat          | Genetic variant representation and storage.                                                                       |

## Software
| Software                                                | Purpose                                                    |
|---------------------------------------------------------|------------------------------------------------------------|
| [agc](https://github.com/refresh-bio/agc)               | Performant FASTA genome compression                        |
| [AnchorWave](https://github.com/baoxingsong/AnchorWave) | Sensitive aligner for genomes with high sequence diversity |
| [bcftools](https://github.com/samtools/bcftools)        | Utilities for indexing VCF data                            |
| [samtools](https://github.com/samtools/samtools)        | bgzip compression for VCF data                             |
| [TileDB](https://github.com/TileDB-Inc/TileDB)          | Performant storage core for array data                     |
| [TileDB-VCF](https://github.com/TileDB-Inc/TileDB-VCF)  | API for storing and querying VCF data                      |
