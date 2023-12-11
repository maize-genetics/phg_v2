# PHG version 2
> [!TIP]
> The quickest way to install PHG2 is from Conda: `conda install -c maize-genetics phg2`

[![PHGv2 CI](https://github.com/maize-genetics/phg_v2/actions/workflows/phgv2_ci.yml/badge.svg)](https://github.com/maize-genetics/phg_v2/actions/workflows/phgv2_ci.yml) [![codecov](https://codecov.io/gh/maize-genetics/phg_v2/graph/badge.svg?token=4BVD2QXQ1A)](https://codecov.io/gh/maize-genetics/phg_v2) [![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)

The Practical Haplotype Graph (PHG) is a powerful tool for representing pangenomes.  The PHG is optimized for the plant breeding and genetics, where genomic diversity can be high, phased haplotypes are common (e.g. inbred lines), and imputation with low density markers is essential for breeding efficiency. This is powerful complement to the excellent tools such as [BEAGLE](https://faculty.washington.edu/browning/beagle/beagle.html) that is used extensively in low diversity, unphased species with high density genotyping.

The PHG is a trellis graph based representation of genic and intergenic regions (called reference ranges) which represent diversity across and between samples. It can be used to: create custom genomes for alignment, call rare alleles, impute genotypes, and efficiently store genomic data from many samples (i.e. reference, assemblies, and other lines). The PHG also works well with community standards including the Breeding API [BrAPI](https://brapi.org) and powerful tools for R such as [rPHG](https://github.com/maize-genetics/rPHG) for pangenome extraction and [rTASSEL](https://github.com/maize-genetics/rTASSEL) for connecting genotype to phenotype.

[PHGv1](https://bitbucket.org/bucklerlab/practicalhaplotypegraph/wiki/Home) was [published in 2022](https://doi.org/10.1093/bioinformatics/btac410). It addressed many challenges related to aligning diverse genomes, efficient storage, and imputation across a pangenome. However, it depended on a custom relational database that necessitated unique formats, and database queries did not scale effectively with a large number of taxa and rare alleles. Moreover, after developing PHGs for six species, we identified significant opportunities to refine and streamline the platform for curation.

## Quick Start
> [!IMPORTANT]
> PHG2 is still under active development, please regularly check back for updates!

You can download PHG2 from [the releases](https://github.com/maize-genetics/phg_v2/releases) page, but we recommend installing through Conda:

```
conda install --channel maize-genetics phg2
```

We recommend creating a dedicated Conda environment for each of your bioinformatics projects, as opposed to installing directly into the base environment; please refer to Conda's documentation [on managing environments](https://conda.io/projects/conda/en/latest/user-guide/getting-started.html#managing-environments). To create an environment and install PHG2 in a single command, do as follows:

```
conda create --name my_env --channel maize-genetics phg2
```

If you've already installed PHG2, you can update to the newest version via:

```
conda update phg2
```

Please refer to the [building and loading documentation](docs/build_and_load.md) for instructions on getting started with using PHG2.

## Design
The redesign leverages the powerful TileDB-VCF database, which is widely used in human genetics for extensive medical applications and is highly performant for rapid querying and storage of rare variants. The PHG is now backed by two TileDB-VCF databases: one for tracking haplotypes across all samples (h.vcf), and another for tracking variants relative to either the reference genomes or the closest haplotype (g.vcf). Our implementation of haplotype encoding in VCF heavily relies on the VCF ALT haplotype specification [v4.3](http://samtools.github.io/hts-specs/VCFv4.3.pdf).

* High-quality phased genome assemblies (or similar) are available to initialize the PHG.
* Ancestral haplotypes are aligned to the reference genome for the identification of haplotypes.
* All PHG tools rely on public file standards - fasta, vcf, bcf, bed, and maf.
* We rely on public tools like TileDB, Minimap2, GATK, AnchorWave, BioKotlin, and HTSJDK.
* Genotyping with low-density markers is now done using a memory- and speed-efficient kmer approach, followed by pathfinding (imputation) with HMM, BWT, or our ML model.
* Rare allele discovery with short reads is based on the above path, involving short read alignment to the inferred haplotype path genome and the GATK haplotype caller.

## Terminology

    Reference genome - the genome used for initial alignment and base coordinates
    Reference range - a segment of the reference genome
    Haplotype - the sequence of part of an individual chromosome.
    Founder Paths - 
    Path - the phased set of haplotypes that represent a chromosome (phased haplotype scaffold in BEAGLE)
    Composite Reference Haplotypes 

More information on terminology can be found [here](docs/terminology.md).
