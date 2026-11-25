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
one sample per founder. This is the normal case for inbred founders. For
heterozygous founders with a **phased** panel, give each haplotype its
own sample, with its own name: the method then treats each haplotype as
a haploid founder. A merged pangenome VCF, produced by running
[`merge-gvcfs`](convenience_commands.md) on the assembly gVCF files of a
PHG, is one source of a high-density panel.

A founder with **no call** at a site of the path panel is scored as if
it were a founder chosen at random from those called there: it carries
each allele at the frequency the called founders carry it. A missing
call is therefore neither evidence for the founder nor against it, so a
founder with many missing calls is not drawn into the path.

The sample VCF is treated as **unphased**: genotypes are compared as
unordered pairs, so `0/1` and `1/0` are the same call. Missing
genotypes (`./.`) are allowed and are simply uninformative, and so is a sample's
half call such as `0/.`: with one allele unknown, it is not read as `0/0`.

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
| `--output-file` | The imputed VCF to write. A name ending in `gz` (e.g. `imputed.vcf.gz`) is written block-compressed, ready for `tabix`; any other name as plain text. | | :material-check: |
| `--bed-dir` | Directory for the founder paths, one `<sampleName>_imputed_path.bed` per sample. Nothing intermediate is written unless this is given. | | |
| `--path-type` | `diploid` infers a pair of founders at each position; `haploid` infers a single founder. See [Choosing a path type](#choosing-a-path-type). | `diploid` | |
| `--prob-correct` | The probability that a genotype call in the sample VCF is correct. | `0.98` | |
| `--prob-switch` | The probability of a path switch (a recombination) across `--prob-switch-distance` base pairs, scaled to each step's actual distance. | `1e-4` | |
| `--prob-switch-distance` | The distance, in base pairs, that `--prob-switch` is quoted over. | `1000000` | |
| `--extend-to-contig-ends` | Extend each contig's path out to both ends of the contig. See [Contig ends](#contig-ends). | off | |
| `--contigs-to-use` | Comma-separated contigs to impute, or a file with one contig per line. All contigs shared by the sample VCF and the path panel if omitted. | | |
| `--keep-sample-calls` | Write each sample's own genotype call over the imputed one wherever it has a full call at a high-density site. See [Keeping the sample's own calls](#keeping-the-samples-own-calls). | off | |
| `--threads` | Number of samples to find paths for at the same time. See [Memory and threads](#memory-and-threads). | `1` | |

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

#### Keeping the sample's own calls

The path reproduces most of a sample's own calls, but not all of them.
A call the path cannot explain, such as a genotyping error, a new
mutation or a founder allele the panel has wrong, is outvoted by the
sites around it. With `--keep-sample-calls`, wherever the sample has a
**full call** at a site of the high-density panel, the output gets the
sample's own call instead of the imputed one. A call is used only when
the site's position and `REF` match and the panel record carries every
allele the sample calls. A half call such as `0/.` is not used, because
the missing allele could be either one.

Every genotype is then marked in a `SRC` FORMAT field: `O` for the
sample's own call (observed) and `I` for an imputed genotype.

A summary is written next to the output, named after it:
`imputed.vcf` gives `imputed_sample_calls.tsv`, with one row per sample:

| Column | Meaning |
|---|---|
| `sample_calls` | Sites where the sample's own call was written. |
| `agreed` | Of those, sites where the imputed genotype was the same. |
| `overwritten` | Of those, sites where the imputed genotype was a different full call. |
| `filled` | Of those, sites with no imputed genotype or only a half call. |
| `unmatched_alleles` | Full calls not used, because the panel record lacks one of their alleles. |
| `half_calls` | Half calls not used. |
| `discordance` | `overwritten / (agreed + overwritten)`: how often the path disagrees with the sample's calls. |

The discordance is an optimistic estimate of the imputation error, since
the same sites also informed the path. A sample whose discordance is much
higher than the rest is worth checking; see
[Founders missing from the panel](#how-well-does-it-work).

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
`--contigs-to-use`, `--threads`), and writes the BED files to the required
`--out-path-dir`. `--keep-sample-calls` is available only in the chained
command.

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
| `--output-file` | The VCF to write; block-compressed when the name ends in `gz`. | | :material-check: |

The chained command and the two separate commands produce the same
output. `bed-to-vcf` also accepts the BED files written by
[`impute-path-from-ps4g`](imputation_ropebwt.md), so it can compose a
VCF from read-based founder paths too.

### Memory and threads

Neither panel is held in memory: both are read one record at a time. The
path finder works one contig at a time, so the memory needed is set by the
**largest contig**, and by:

* **Shared sites**, `S`: the sites on that contig that the sample VCF has in
  common with the path panel. The sample VCF's density decides this, not
  the size of either panel.
* **Path states**: founders² for a diploid path (6,724 for 82 founders), or
  founders for a haploid one.
* **Founders and samples**: their alleles at the shared sites of the
  current contig, about `4 × (founders + samples)` bytes per site.

Finding a path keeps a traceback of 4 bytes per state per shared site,
for each sample being worked on. This is usually the largest term, and it
is **capped**. Each path may use up to a quarter of the Java heap, divided
among the `--threads`. A longer path keeps checkpoints instead and rebuilds
its traceback one block at a time. That gives exactly the same path, in
about twice the time. So any contig can be imputed: the heap decides only
whether it runs at full speed.

Roughly, the heap needed is

```
0.1 GB + S × 4 × (founders + samples) bytes                    at the least
     + threads × 16 × S × states bytes                           for full speed
```

Allow about 50% more than this when setting `-Xmx`
([Setting memory](installation.md#setting-memory)).

Some examples for a diploid path with 1,000 samples and one thread, where
the largest chromosome holds about 14% of the sites:

| Sample data | Founders | Shared sites, largest chromosome | Heap, at the least | Heap, full speed |
|---|---|---|---|---|
| 50K SNP array | 25 | about 7,000 | about 0.15 GB | about 0.2 GB |
| 600K SNP array | 25 | about 85,000 | about 0.5 GB | about 1.3 GB |
| 10 million variants (e.g. low-pass sequencing) | 25 | about 1.4 million | about 6 GB | about 20 GB |
| the same | 100 | about 1.4 million | about 6.5 GB | about 230 GB |
| 2.4 million variants | 82 haplotypes | about 340,000 | about 1.6 GB | about 38 GB |

`--threads` runs that many samples at once. Each needs its own
traceback, so more threads share the same quarter of the heap, and the
output does not depend on the thread count.

!!! tip "Reducing memory and time"
    * Use `--path-type haploid` for fully inbred material. This shrinks the
      largest term by a factor of the number of founders, with no loss of
      accuracy for inbred samples.
    * Thin a very dense sample VCF before imputing. The paths only need
      enough markers to place crossovers (array density imputed 99.8% of
      genotypes correctly in our tests), and the high-density panel still
      supplies every output site.
    * Splitting the samples into batches reduces only the smaller,
      per-sample term.

Running time grows with shared sites × states for each sample. In our
tests, 1,000 samples with 106,705 shared sites on a contig took about 8.5
minutes for a diploid path with 25 founders.

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
