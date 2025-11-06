# PS4G - **P**ositional **S**upport **for** **G**ametes Specification

* **_Specification version:_** `v2.0`
* **_Date:_** 2025-10-24

## Overview
A PS4G (**P**ositional **S**upport **for** **G**ametes) file is a standardized
tab-delimited format that tracks genomic support for reference panel gametes
across binned genomic positions. This format provides positional tracking for
genomic data from multiple sources (e.g., variants, read alignments), thereby
enhancing pathfinding support and enabling integration with machine learning-based
imputation engines.

PS4G files aggregate evidence from sequencing reads or genotype calls to determine
which reference panel haplotypes (gametes) are supported at different genomic
locations. This compressed representation enables efficient downstream imputation
and haplotype inference.

## File Format

### Structure
A PS4G file consists of:

1. **Header section** - Metadata lines prefixed with `#`
2. **Column header line** - Tab-delimited field names
3. **Data section** - Tab-delimited data lines

### Header Section
The header contains metadata about the file creation and reference gametes:

```
#PS4G
#version=2.0
#<metadata lines>
#Command: <CLI command used to generate this file>
#TotalUniqueCounts: <sum of all counts in the file>
#gamete	gameteIndex	count
#<SampleGamete>	<index>	<total_count>
...
```

**Header fields:**

| Field | Description |
|-------|-------------|
| `#PS4G` | File format identifier (required first line) |
| `#version` | PS4G format version (2.0 for current specification) |
| `#Command` | Full CLI command used to generate the file |
| `#TotalUniqueCounts` | Sum of all unique position counts in the file |
| `#gamete` | Reference panel gamete identifier (format: `SampleName` or `SampleName:gameteIndex`) |
| `gameteIndex` | Zero-based integer index assigned to each gamete |
| `count` | Total number of observations supporting this gamete across all positions |

### Data Section
After the header, the column names are defined, followed by data rows:

```
gameteSet	refContig	refPosBinned	count
<comma-separated gamete indices>	<contig>	<binned position>	<count>
```

**Data fields:**

| Field | Type | Description |
|-------|------|-------------|
| `gameteSet` | String | Comma-separated list of gamete indices (from header) that are supported at this position |
| `refContig` | String | Reference contig/chromosome identifier |
| `refPosBinned` | Integer | Binned genomic position (actual position divided by 256) |
| `count` | Integer | Number of reads/variants supporting this gamete set at this position |

### Example
```
#PS4G
#version=2.0
#Command: phg convert-rm2ps4g-file --read-mapping-file input.txt --hvcf-dir /path/to/hvcfs --output-dir output/
#TotalUniqueCounts: 1234
#gamete	gameteIndex	count
#LineA	0	853
#LineB	1	381
#Ref	2	100
gameteSet	refContig	refPosBinned	count
0	chr1	1000	853
0,1	chr1	2000	24
1	chr2	500	15
0,1,2	chr2	1500	5
```

In this example:

- Row 1: Gamete 0 (`LineA`) is supported at chr1 binned position 1000 (actual position ~256,000 bp) by 853 reads
- Row 2: Gametes 0 and 1 (`LineA` and `LineB`) are both supported at chr1 position 2000 (~512,000 bp) by 24 reads
- Row 3: Gamete 1 (`LineB`) is supported at chr2 position 500 (~128,000 bp) by 15 reads
- Row 4: All three gametes are supported at chr2 position 1500 (~384,000 bp) by 5 reads

## Position Binning

To reduce file size and provide efficient storage, genomic positions are binned into 256 bp windows:

**Binning process:**

The `refPosBinned` field stores the genomic position divided by 256 (integer division):

```
refPosBinned = genomicPosition / 256
```

**Converting back to approximate genomic position:**

```
approximateGenomicPosition = refPosBinned * 256
```

!!! note "Resolution and Compression"
    The 256 bp binning provides a balance between positional resolution and file compression.
    The actual genomic position is rounded down to the nearest 256 bp boundary during binning.
    This resolution is suitable for chromosome-scale imputation and pathfinding algorithms.

!!! example "Binning Example"
    - Genomic position 256,000 bp → Binned position 1,000
    - Genomic position 256,255 bp → Binned position 1,000 (same bin)
    - Genomic position 256,256 bp → Binned position 1,001 (next bin)

## Generation Methods

PS4G files can be generated from three different sources:

### 1. From Read Mapping Files (`convert-rm2ps4g-file`)

Converts PHG read mapping output (from k-mer or RopeBWT mapping) to PS4G format.

**Input:** Read mapping file with format:
```
HapIds	count
<comma-separated haplotype IDs>	<count>
```

**Process:**

1. Maps haplotype IDs to reference ranges
2. Identifies gametes at each reference range
3. Determines position from reference range start
4. Aggregates counts by gamete set and position

**Command:**
```shell
phg convert-rm2ps4g-file \
    --read-mapping-file mapping.txt \
    --hvcf-dir /path/to/hvcfs \
    --output-dir output/
```

### 2. From RopeBWT BED Files (`convert-ropebwt2ps4g-file`)

Converts RopeBWT3 maximal exact match (MEM) alignments to PS4G format using spline-based coordinate transformation.

**Input:** BED file from RopeBWT3 with MEM alignments

**Process:**

1. Loads spline knots for coordinate transformation from assembly to reference coordinates
2. Groups MEMs by read
3. Filters by minimum MEM length and maximum hits
4. Uses spline interpolation to map assembly positions to reference positions
5. Creates consensus position from multiple MEMs
6. Aggregates gamete support by position

**Command:**
```shell
phg convert-ropebwt2ps4g-file \
    --ropebwt-bed alignments.bed \
    --spline-knot-dir /path/to/splines \
    --output-dir output/ \
    --min-mem-length 148 \
    --max-num-hits 50
```

**Parameters:**

- `--min-mem-length`: Minimum MEM length to consider (default: 148)
- `--max-num-hits`: Maximum number of haplotype hits allowed (default: 50)
- `--sort-positions`: Sort output by genomic position (default: true)

### 3. From VCF Files (`convert-vcf2ps4g-file`)

Converts variant calls to PS4G format using a reference panel for gamete matching.

**Input:**

- Sample VCF file (to be imputed)
- Reference panel VCF file

**Process:**

1. Builds allele-to-gamete lookup from reference panel
2. For each variant in sample VCF:
   - Matches alleles to reference panel gametes
   - Records gamete support at that position
3. Aggregates counts by gamete set and position

**Command:**
```shell
phg convert-vcf2ps4g-file \
    --to-impute-vcf sample.vcf \
    --ref-panel-vcf reference_panel.vcf \
    --output-dir output/
```

## Data Interpretation

### Gamete Sets
The `gameteSet` field contains indices of gametes that share evidence at a position. Multiple gametes in a set indicate:

- **From read data**: Reads mapping ambiguously to multiple haplotypes
- **From VCF data**: Shared alleles across multiple reference samples

### Position Accuracy
Due to the 256 bp binning:

- Positions represent approximate genomic locations
- Multiple nearby variants/reads may contribute to the same bin
- Suitable for chromosome-scale imputation, not for fine-scale variant calling

### Count Interpretation
The `count` field represents:

- **Read mapping**: Number of reads supporting this gamete combination
- **VCF conversion**: Number of variants matching this pattern
- Higher counts indicate stronger evidence for those gametes

## Use Cases

PS4G files are designed for:

1. **Machine learning-based imputation**: Provide feature vectors for ML models to predict haplotypes
2. **Pathfinding algorithms**: Inform hidden Markov models about gamete support across the genome
3. **Quality control**: Assess read mapping quality and reference panel coverage
4. **Comparative analysis**: Compare imputation results across different methods

## File Naming Convention

PHG generates PS4G files with the naming pattern:
```
<input_basename>_<sampleGamete>_ps4g.txt
```

Examples:

- `LineA_1_readMapping_ps4g.txt` - From read mapping
- `sample_alignments_ps4g.txt` - From RopeBWT
- `input_vcf_ps4g.txt` - From VCF conversion

## Related Commands

- `rope-bwt-chr-index` - Create RopeBWT index for alignment
- `build-spline-knots` - Generate spline knots for coordinate transformation
- `map-reads` - Align reads to generate mapping files
- See [Imputation using Machine Learning](imputation_ml.md) for complete workflow

## Specification Notes

### Version History

- **v2.0** (2025-10-24): Major format update - removed position encoding, split position into separate `refContig` and `refPosBinned` columns for improved readability and flexibility
- **v1.0** (2025-02-19): Initial complete specification
- **v0.1** (2025-02-19): Draft specification

### Implementation
PS4G files are generated by PHGv2 commands in the `net.maizegenetics.phgv2.pathing.ropebwt` package:

- `ConvertRm2Ps4gFile.kt` - Read mapping conversion
- `ConvertRopebwt2Ps4gFile.kt` - RopeBWT conversion
- `ConvertVcf2Ps4gFile.kt` - VCF conversion
- `PS4GUtils.kt` - Shared utilities for file writing and data formatting
