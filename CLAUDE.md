# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

PHGv2 (Practical Haplotype Graph v2) is a Kotlin/JVM bioinformatics tool for representing and querying pangenomes, optimized for plant breeding. It stores haplotypes across samples in TileDB-VCF databases and supports imputation via k-mer mapping + HMM pathfinding.

Current version: `2.5.2` (see `version.properties`).

## Build & Test Commands

Requires JDK 21 (jvmToolchain is pinned to 21; CI uses Temurin 21).

Some `brapi` tests read `src/main/resources/application.conf` for a `TILEDB_URI`. If those tests fail locally, run `src/scripts/update_applications_conf.sh` first (it prepends `TILEDB_URI=$HOME/temp/phgv2Tests/tempDir/testTileDBURI/` to the conf file) — this is also a step in CI.

```bash
# Build
./gradlew build

# Build without tests
./gradlew build -x test

# Run all tests
./gradlew test

# Run a single test class
./gradlew test --tests "net.maizegenetics.phgv2.cli.CreateRangesTest"

# Run a specific test method
./gradlew test --tests "net.maizegenetics.phgv2.cli.CreateRangesTest.testCreateRanges"

# Clean build
./gradlew clean build

# Code coverage report (Kover, minimum 70% line coverage enforced)
./gradlew koverXMLReport

# Create distribution tarball
./gradlew distTar
```

Tests use JUnit 5 (`useJUnitPlatform()`). Test heap is set to 50G in `build.gradle.kts` — reduce if needed locally.

## Architecture

### Entry Point

`src/main/kotlin/net/maizegenetics/phgv2/cli/Phg.kt` — the main CLI class using [Clikt](https://github.com/ajalt/clikt). All subcommands are registered here via `.subcommands(...)`.

### Package Structure

| Package | Purpose |
|---|---|
| `cli/` | All CLI commands (one file per command). Each extends `CliktCommand`. |
| `api/` | Core data structures: `HaplotypeGraph`, `ReferenceRange`, `SampleGamete`, `TileDBVcfReader` |
| `pathing/` | K-mer indexing, read mapping, HMM-based path finding (imputation), QC |
| `pathing/ropebwt/` | RopeBWT3-based alignment workflow for imputation (alternative to k-mer) |
| `utils/` | Shared utilities: file I/O, VCF/BED/GFF parsing, TileDB queries, logging |
| `brapi/` | BrAPI server (Ktor-based): `Application.kt`, `api/`, `model/`, `service/`, `utilities/` |
| `rphg/` | Utilities supporting the R package rPHG2 |
| `scaffolding/` | Scaffolding-related tools |
| `kmers/` | K-mer related utilities |

### Core Data Flow

1. **Build**: Assemblies aligned to reference → MAF files → HVCF files → loaded into TileDB
2. **Impute**: Export HVCFs → build RopeBWT/k-mer index → map reads → `FindPaths` (Viterbi HMM) → imputed HVCFs
3. **Resequencing**: Imputed path → align reads to path genome → GATK haplotype calling

### Key Data Structures

- **`HaplotypeGraph`** (`api/HaplotypeGraph.kt`): In-memory graph built from `.h.vcf` files. Core lookup structure mapping `ReferenceRange × SampleGamete → haplotype checksum`.
- **`ReferenceRange`** (`api/ReferenceRange.kt`): Genomic interval (contig + start + end).
- **HVCF**: VCF files where ALT fields encode haplotype checksums per the VCF 4.2 ALT haplotype spec. Extensions: `.h.vcf` or `.h.vcf.gz`.
- **TileDB**: Two databases — one for haplotypes (HVCF), one for variants (GVCF).

### TileDB Dependency

The project uses a local `.jar` for TileDB-VCF (`repo/tiledb-vcf-java-0.37.0-1-g03553439.jar`). The `build.gradle.kts` has commented alternatives for Mac Intel/ARM. **The active jar is the Linux build** — swap comments in `build.gradle.kts` for other platforms.

### Conda Environments

Two conda environments are required for full functionality:
- `phgv2-conda` (`src/main/resources/phg_environment.yml`): AnchorWave, bcftools, samtools, AGC, ropebwt3, minimap2
- TileDB conda env (`src/main/resources/phg_tiledb_environment.yml`): TileDB-VCF CLI tools

Setup: `./phg setup-environment` (calls conda internally).

### CI/CD

CI runs on PRs to `main` targeting `src/**` changes (`.github/workflows/phgv2_ci.yml`). It runs `./gradlew clean build` which includes all tests. Code coverage is reported to Codecov via Kover (minimum 70% line coverage).

## Adding a New CLI Command

1. Create a new file in `src/main/kotlin/net/maizegenetics/phgv2/cli/` extending `CliktCommand`.
2. Register it in the `main()` function in `Phg.kt`.
3. Add a corresponding test in `src/test/kotlin/net/maizegenetics/phgv2/cli/`.
4. Test data goes in `data/test/`.
