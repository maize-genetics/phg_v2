# Imputing a high-density VCF from a lower-density VCF

!!! note
    Please [let us know](https://github.com/maize-genetics/phg_v2/issues/new/choose)
    if you have any questions or issues.

In this document, we will discuss how to raise the marker density of
genotyped samples using a reference panel of their founders. The
starting point is a VCF of genotype calls, such as a SNP array or a
low-coverage genotyping assay, rather than sequencing reads. The steps:

1. Infer each sample's **founder path**: which panel founders it
   descends from, position by position along each chromosome
2. Read the founders' alleles out of a **denser panel** to fill in the
   sample's genotype at every site that panel carries

A single command, `impute-vcf-from-vcf`, does both. The two halves are
also available as separate commands, `impute-path-from-vcf` and
`bed-to-vcf`, for when you want to inspect or reuse the founder paths.

!!! tip "Which imputation method should I use?"
    Use this method when your samples are already **genotyped** and you
    have a VCF of their calls. If you have **sequencing reads**, use
    [RopeBWT3 imputation](imputation_ropebwt.md) instead. That method
    works from the reads directly and does not require a genotyping
    step first.

## Quick start

* Impute the samples in one step:
  ```shell
  phg impute-vcf-from-vcf \
      --to-impute-vcf /my/samples.vcf \
      --panel-vcf /my/founder_panel.vcf \
      --high-density-panel-vcf /my/founder_panel_dense.vcf \
      --output-file /my/imputed.vcf
  ```

* OPTIONAL: Keep the founder paths as BED files by adding
  `--bed-dir /my/path/dir`.

* OPTIONAL: Run the two halves separately:
  ```shell
  phg impute-path-from-vcf \
      --to-impute-vcf /my/samples.vcf \
      --panel-vcf /my/founder_panel.vcf \
      --out-path-dir /my/path/dir

  phg bed-to-vcf \
      --bed-dir /my/path/dir \
      --reference-panel-vcf /my/founder_panel_dense.vcf \
      --output-file /my/imputed.vcf
  ```

## Detailed walkthrough

### Input files

Three VCF files are needed. They can be two, since the path panel and
the high-density panel may be the same file.

| Input | Contents | Option |
|---|---|---|
| Sample VCF | The genotyped samples to impute. | `--to-impute-vcf` |
| Path panel | One sample per founder, genotyped at (at least) the sample VCF's sites. Used to infer each sample's founder path. | `--panel-vcf` |
| High-density panel | One sample per founder, at every site you want in the output. Supplies the output's sites and alleles. | `--high-density-panel-vcf` (`--reference-panel-vcf` in `bed-to-vcf`) |

All three must:

* be **coordinate sorted**. The sample VCF and the path panel must also
  list their contigs in the same order.
* use the **same reference genome and the same contig names**. Sites
  are matched on contig, position and reference allele, so a site whose
  `REF` differs between files is treated as absent from the panel.

The panels should be **haploid, or diploid and mostly homozygous**, with
one sample per founder. This is the normal case for inbred founders. A
merged pangenome VCF, produced by running
[`merge-gvcfs`](convenience_commands.md) on the assembly gVCF files of a
PHG, is one source of a high-density panel.

The sample VCF is treated as **unphased**: genotypes are compared as
unordered pairs, so `0/1` and `1/0` are the same call. Missing
genotypes (`./.`) are allowed and are simply uninformative.

Every founder in the path panel must also be in the high-density panel,
under the same name. This is checked before any work starts. Extra
founders in the high-density panel are allowed; they are never used.

!!! note "Using one panel for both roles"
    You can pass the same file as `--panel-vcf` and
    `--high-density-panel-vcf`. Only the sites the path panel shares
    with the sample VCF are used to infer the path, so a dense panel
    works as a path panel. A separate path panel is useful when the
    founders were genotyped on the same platform as the samples, so
    that both are called the same way at the shared sites.

### Impute a high-density VCF

```shell
phg impute-vcf-from-vcf \
    --to-impute-vcf output/samples_array.vcf \
    --panel-vcf output/founders_array.vcf \
    --high-density-panel-vcf output/founders_dense.vcf \
    --output-file output/samples_imputed.vcf
```

For each sample, the command finds the most likely path through the
founders along each contig using a hidden Markov model (HMM). At each
shared site, the model scores every pair of founders by how well the
genotypes that pair could produce match the sample's observed genotype.
Between sites, it charges a penalty for switching founders, which
grows with the distance between the sites. The path is then used to
look up each founder's allele at every site in the high-density panel.

The founder paths are held in memory and nothing intermediate is
written, unless you add `--bed-dir`.

#### Parameters

| Parameter name | Description | Default value | Required? |
|---|---|---|---|
| `--to-impute-vcf` | VCF of the samples to impute. | | :material-check: |
| `--panel-vcf` | Founder panel used to infer the path. Its samples must be a subset of the high-density panel's. | | :material-check: |
| `--high-density-panel-vcf` | Founder panel supplying the output's sites and alleles. | | :material-check: |
| `--output-file` | The imputed VCF to write. | | :material-check: |
| `--bed-dir` | Directory for the founder paths, one `<sampleName>_imputed_path.bed` per sample. Nothing intermediate is written unless this is given. | | |
| `--path-type` | `diploid` infers a pair of founders at each position; `haploid` infers a single founder. See [Choosing a path type](#choosing-a-path-type). | `diploid` | |
| `--prob-correct` | The probability that a genotype call in the sample VCF is correct. | `0.98` | |
| `--prob-switch` | The probability of a path switch (a recombination) across `--prob-switch-distance` base pairs, scaled to each step's actual distance. | `1e-4` | |
| `--prob-switch-distance` | The distance, in base pairs, that `--prob-switch` is quoted over. | `1000000` | |
| `--extend-to-contig-ends` | Extend each contig's path out to both ends of the contig. See [Contig ends](#contig-ends). | off | |
| `--contigs-to-use` | Comma-separated contigs to impute, or a file with one contig per line. All contigs shared by the sample VCF and the path panel if omitted. | | |

#### Choosing a path type

* **`diploid`** (the default) infers a pair of founders at each
  position, so it can represent heterozygous samples such as F1
  hybrids and F2s.
* **`haploid`** infers a single founder, and is only suitable for
  fully inbred material such as doubled haploids and inbred lines.

!!! tip
    If you are unsure, use `diploid`. It imputes inbred lines as
    accurately as `haploid`, while `haploid` cannot represent a
    heterozygote and gives wrong genotypes for heterozygous samples.

#### Tuning the probabilities

The defaults work well in our tests, and accuracy changes very little
across a wide range of values (see
[How well does it work?](#how-well-does-it-work)).

* `--prob-switch` controls how readily the path changes founders.
  Although it is phrased as a recombination probability, it acts as a
  smoothing parameter. Lower values give fewer, longer founder
  segments; higher values let the path follow the data more closely.
* `--prob-correct` sets the cost of a site where the sample's genotype
  disagrees with the founder pair. Lower values make the path more
  tolerant of genotyping errors.

#### Contig ends

By default, each contig's path runs from the first to the last site the
sample shares with the path panel. High-density sites outside that span
are left uncalled (`./.`), because there is no evidence of ancestry
beyond the outermost markers.

`--extend-to-contig-ends` carries the first and last founder segments
out to the ends of the contig instead. Use it when your markers come
close to the contig ends, or when a complete call set matters more than
the risk of an undetected crossover near a chromosome end. The contig
lengths are read from the path panel's `##contig` header lines; if there
are none, the last segment is extended past any site the high-density
panel could hold.

### Output

The output VCF contains every site in the high-density panel and one
column per sample, sorted by sample name. Genotypes are written
unphased, with one allele from each founder in the sample's pair.

A founder contributes a no-call (`.`) for its haplotype where it has no
single allele:

* where the founder has no call in the high-density panel, or
* where the founder is **heterozygous** in the high-density panel, since
  it is not known which of its two alleles the sample inherited.

So a sample whose founders are one homozygous and one heterozygous at a
site is written as a half call, such as `1/.`.

!!! note
    Records are copied from the high-density panel, including its `INFO`
    fields. Fields such as `AF`, `DP` and `NS` describe the panel, not
    the imputed samples.

### Log diagnostics

The log reports how the sample VCF and the path panel matched up:

* **Sites used from the sample VCF**: the sites the paths were
  inferred from.
* **Sites in the sample VCF not present in the panel**: these are
  dropped. A large number, reported with a warning when more than half
  the sample's sites are missing, usually means the two files use
  different reference genomes, contig names or allele representations.
* **Sites in the panel not present in the sample VCF**: expected when
  the panel is denser than the samples.

### Founder paths (optional)

With `--bed-dir`, or from `impute-path-from-vcf`, each sample's path is
written to `<sampleName>_imputed_path.bed`. The coordinates are 0-based
and half-open, following the BED convention.

A diploid path has five columns, one founder per haplotype:

```
chrom	start	end	parent1	parent2
chr1	99999	41820554	LineA	LineB
chr1	41820554	216309981	LineA	LineA
```

A haploid path has four:

```
chrom	start	end	parent1
chr1	99999	216309981	LineA
```

### Running the two steps separately

`impute-path-from-vcf` infers and writes the founder paths:

```shell
phg impute-path-from-vcf \
    --to-impute-vcf output/samples_array.vcf \
    --panel-vcf output/founders_array.vcf \
    --out-path-dir output/founder_paths
```

It takes the same path-finding parameters as `impute-vcf-from-vcf`
(`--path-type`, `--prob-correct`, `--prob-switch`,
`--prob-switch-distance`, `--extend-to-contig-ends`,
`--contigs-to-use`), and writes the BED files to the required
`--out-path-dir`.

`bed-to-vcf` then composes the paths into a VCF:

```shell
phg bed-to-vcf \
    --bed-dir output/founder_paths \
    --reference-panel-vcf output/founders_dense.vcf \
    --output-file output/samples_imputed.vcf
```

| Parameter name | Description | Default value | Required? |
|---|---|---|---|
| `--bed-dir` | Directory of founder-path BED files. Recognized names are `<sample>_imputed_path.bed`, `<sample>_chr<contig>_imputed.bed` and `<sample>.bed`; several files for one sample are merged. | | :material-check: |
| `--reference-panel-vcf` | Founder panel supplying the output's sites and alleles. Its sample names must match the founder names in the BED files. | | :material-check: |
| `--output-file` | The VCF to write. | | :material-check: |

The chained command and the two separate commands produce the same
output. `bed-to-vcf` also accepts the BED files written by
[`impute-path-from-ps4g`](imputation_ropebwt.md), so it can compose a
VCF from read-based founder paths too.

## How well does it work?

We tested the method on simulated maize samples built as mosaics of
real founders with known crossovers. Chromosome 10 of a 25-founder panel
was reduced to 3,600 array-like SNPs, with 1% genotyping error and 5%
missing calls, and imputed back to 106,705 sites. At sites where a
sample's two parents differ:

| Population | Correct genotypes |
|---|---|
| F1 hybrids | 100.0% |
| F2 | 99.8% |
| Doubled haploids | 99.9% |
| Recombinant inbred lines | 99.8% |

Nearly all the remaining errors lie within 2 Mb of a crossover, where
its exact position falls between markers. Raising the genotyping error
to 5% and the missing calls to 15% cost less than 0.05 percentage
points. Imputing all 200 samples took about 10 seconds.

!!! warning "Founders missing from the panel"
    The method can only assign founders that are in the panel. Where a
    sample descends from a founder the panel lacks, its path patches
    together the most similar founders available. In our tests,
    accuracy at those positions fell to 60–80%, and the path switched
    founders far more often than the true crossovers. The imputed genotypes from those samples fit
    their own input genotypes noticeably worse than the others. If a
    few samples have many more path switches than the rest, or disagree
    with their own input genotypes at the shared sites much more often,
    check whether all of their parents are in the panel.
