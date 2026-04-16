!!! warning "Legacy Documentation - PHG Version 1"

    This section contains documentation for **PHG version 1**, which is
    no longer actively developed. It is preserved here for archival and
    historical reference only. If you are looking to use the Practical
    Haplotype Graph, please refer to the [PHG v2 documentation](../../index.md),
    which reflects the current version of the software.

# Quick Start

Aligning assemblies to the reference genome and adding the resulting haplotypes to the PHG database 
is can be done by the AssemblyHaplotypesMultiThreadPlugin. It can be run separately using the following command
with the maximum memory size (-Xmx option), the path to the config file, and the log file set correctly:

```
docker run --rm  maizegenetics/phg /tassel-5-standalone/run_pipeline.pl \
-Xmx50G -debug -configParameters /phg/config.txt \
-AssemblyHaplotypesMultiThreadPlugin -endPlugin > /path/to/logfile.txt
```

Note that the log file path is outside the docker. If the log file is not supplied, then the log information
will be written to the terminal (sysout).

The config file must contain the following entries to run this step. The first 4 entries pertain to database access, the remaining entries are for processing your assembly data. Replace the entries with your own parameter values.  The entries listed here assume you have setup your directory structure using the MakeDefaultDirectoryPlugin.


```
#!java
host=localHost
user=sqlite
password=sqlite
DB=/phg/phgAssemblyPath.db
DBtype=sqlite

#liquibase output directory/general output dir
outputDir=/phg/outputDir
liquibaseOutdir=/phg/outputDir/

asmMethodName=mummer4
asmKeyFile=/phg/load_asm_genome_key_file.txt

AssemblyHaplotypesMultiThreadPlugin.outputDir=/phg/outputDir/align/
AssemblyHaplotypesMultiThreadPlugin.keyFile=/phg/load_asm_genome_key_file.txt

```


# Details

This class runs mummer4 alignment of chrom-chrom fastas in a parallel fashion.
It expects a directory for both reference and assembly fastas, split by
chromosome, as well as a keyFile.

The keyfile provides information on the location and name of the full assembly genome file that was used for alignment.  This data lives in the "AssemblyServerDir" and "AssemblyGenomeFasta" columns.  The PHG mummer4 code aligns the ref and assembly on a per-chromosome basis.  The location of the directory and file names for the ref and assembly per-chromosome fastas  are denoted in the RefDir/RefFasta columns (for the reference) and AssemblyDir/AssemblyFasta columns (for the assembly).  The chromosome name and assembly name for loading
to the database are also specified in the keyFile.  The file contains 8 columns and should be tab-delimited format.  The header columns should be as below:

`AssemblyServerDir RefDir RefFasta AssemblyDir AssemblyGenomeFasta AssemblyFasta   AssemblyDBName  Chromosome`

An example key file that has information for a single assembly genome is shown below.  Note all fields are tab-delimited and because this is run in a docker, the reference and assembly directories are docker paths.  The AssemblyServerDir should be an accessible location outside the docker.

```
AssemblyServerDir       RefDir  RefFasta        AssemblyDir     AssemblyGenomeFasta     AssemblyFasta   AssemblyDBName  Chromosome
https://download.maizegdb.org/Zm-CML103-REFERENCE-NAM-1.0/  /phg/inputDir/reference/        B73chr1.fa  /phg/inputDir/assemblies/       Zm-CML103-REFERENCE-NAM-1.0.fa.gz        CML103chr1.fa        CML103_Assembly  1
https://download.maizegdb.org/Zm-CML103-REFERENCE-NAM-1.0/  /phg/inputDir/reference/        B73chr2.fa  /phg/inputDir/assemblies/       Zm-CML103-REFERENCE-NAM-1.0.fa.gz        CML103chr2.fa        CML103_Assembly  2
https://download.maizegdb.org/Zm-CML103-REFERENCE-NAM-1.0/  /phg/inputDir/reference/        B73chr3.fa  /phg/inputDir/assemblies/       Zm-CML103-REFERENCE-NAM-1.0.fa.gz        CML103chr3.fa        CML103_Assembly  3

```

An example key file that has per-chromosome fastas for multiple species is shown in the next example.  Note that the assembly chromosome fasta files must exist in the same directory.

```
AssemblyServerDir       RefDir  RefFasta        AssemblyDir     AssemblyGenomeFasta     AssemblyFasta   AssemblyDBName  Chromosome
https://download.maizegdb.org/Zm-CML103-REFERENCE-NAM-1.0/  /phg/inputDir/reference/        B73chr1.fa  /phg/inputDir/assemblies/       Zm-CML103-REFERENCE-NAM-1.0.fa.gz        CML103chr1.fa        CML103_Assembly  1
https://download.maizegdb.org/Zm-CML103-REFERENCE-NAM-1.0/  /phg/inputDir/reference/        B73chr2.fa  /phg/inputDir/assemblies/       Zm-CML103-REFERENCE-NAM-1.0.fa.gz        CML103chr2.fa        CML103_Assembly  2
https://download.maizegdb.org/Zm-CML103-REFERENCE-NAM-1.0/  /phg/inputDir/reference/        B73chr3.fa  /phg/inputDir/assemblies/       Zm-CML103-REFERENCE-NAM-1.0.fa.gz        CML103chr3.fa        CML103_Assembly  3
https://download.maizegdb.org/Zm-Ki3-REFERENCE-NAM-1.0/  /phg/inputDir/reference/        B73chr1.fa  /phg/inputDir/assemblies/       Zm-Ki3-REFERENCE-NAM-1.0.fa.gz        Ki3chr1.fa        Ki3_Assembly  1
https://download.maizegdb.org/Zm-Ki3-REFERENCE-NAM-1.0/  /phg/inputDir/reference/        B73chr2.fa  /phg/inputDir/assemblies/       Zm-Ki3-REFERENCE-NAM-1.0.fa.gz        Ki3chr2.fa        Ki3_Assembly  2
https://download.maizegdb.org/Zm-Ki3-REFERENCE-NAM-1.0/  /phg/inputDir/reference/        B73chr3.fa  /phg/inputDir/assemblies/       Zm-Ki3-REFERENCE-NAM-1.0.fa.gz        Ki3chr3.fa        Ki3_Assembly  3
```


The mummer4 alignment scripts will be kicked off in parallel based on
the number of threads the user specifies.  Default is 3, which is essentially
single threaded and 2 of the 3 are reserved for I/O.

Note that running the alignments in parallel will take not just threads, but
memory.  This should be considered when the user decides on the number of threads.  It has also been found that running with too many threads may cause contention for DB writing access.  If problems arise when using a large number of threads, consider reducing the number of threads.

## Parameters

Run the following command to print a current list of parameters for this plugin:

```
docker run --rm  maizegenetics/phg /tassel-5-standalone/run_pipeline.pl -AssemblyHaplotypesMultiThreadPlugin
```

Running the command will print the following usage statement:

```
Usage:
AssemblyHaplotypesMultiThreadPlugin <options>
-outputDir <Output Directory> : Output directory including trailing / for writing files (required)
-keyFile <keyFile> : Name of the Keyfile to process.  Must have columns RefDir, RefFasta, AssemblyDir, AssemblyFasta, AssemblyDBName, and Chromosome.   (required)
-mummer4Path <Mummer4 Executables Path> : Path where mummer4 executables live: nucmer, delta-filter, show-snps, show-coords   (Default: /mummer/bin/)
-clusterSize <Mummer4 Nucmer Cluster Size > : Cluster size to use with mummer4 nucmer script.  (Default: 250)
-gvcfOutputDir <GVCF Output Dir> : Directory for gvcf files to be output for later use
-entryPoint <Assembly Entry Point> : Where to begin processing. All runs everything.  Refilter means run from the re-filtering of the coords files. hapltypes runs just the haplotypes processing. (Default: all)
-minInversionLen <Minimum Inversion Length> : Minimum length of inversion for it to be kept as part of the alignment. Default is 7500 (Default: 7500)
-loadDB <true | false> : Boolean: true means load haplotypes to db, false means do not populate the database.  Defaults to true (Default: true)
-assemblyMethod <Assembly Method Name> : Name to be stored to methods table for assembly method, default is mummer4   (Default: mummer4)
-numThreads <Num Threads> : Number of threads used to process assembly chromosomes.  The code will subtract 2 from this number to have the number of worker threads.  It leaves 1 thread for IO to the DB and 1 thread for the Operating System. (Default: 3)
```

Description of required parameters:

* keyfile: The keyfile must be a tab-delimited file with 1 row for each per-chromosome assembly fasta you wish to align and load to the database.  The required columns for the keyfile are AssemblyServerDir,RefDir, RefFasta,  AssemblyDir, AssemblyGenomeFasta, AssemblyFasta, AssemblyDBName and Chromosome.  The oolumn descriptions are below.
    * AssemblyServerDir: machine and directory where the full genome fasta for the assembly may be found.  This should be data on a publicly available server, not a docker specific path.
    * RefDir: full path to the directory containing the reference fasta file.
    * RefFasta: name of the reference fasta file to be aligned to the assembly fasta file in this row.  This should be a fasta file at the chromosome level.  It is the reference file used for the ref-chromosome to assembly-chromosome alignment
    * AssemblyDir: full path to the directory containing the assembly fasta files.
    * AssemblyGenomeFasta:  name of the full genome fasta for this assembly that can be found at the location specified in the AssemblyServerDir column.  This is the file stored to the genome_file_data table.
    * AssemblyFasta:  name of the assembly fasta file to be aligned to the reference fasta file in this row.  This should be a fasta file at the chromosome level.  It is the file that will be used for the ref-chromosome to assembly-chromosome alignment.
    * AssemblyDBName: the name for this assembly to add as the "line_name" to the database's genotype table.
    * Chromosome: the chromosome name that is being processed.  This name must match a reference chromosome that was loaded with the LoadGenomeIntervals.sh script.
* configFile:  DB Config File containing properties host,user,password,DB and DBtype where DBtype is either "sqlite" or "postgres". This file is expected to live in the folder mounted to Docker tempFileDir/data/. This file may also contain values for optional parameters.

The optional parameters that may be defined in the configuration file are  mummer4Path, clusterSize, entryPoint, minInversionLen, loadDB, assemblyMethod and numThreads.
If you wish to set a value other than the default for these parameters, they
should be included in the configuration file in the format shown below (with your own values replacing the default values shown here):

* AssemblyHaplotypesMultiThreadPlugin.mummer4Path=/mummer/bin
* AssemblyHaplotypesMultiThreadPlugin.clusterSize=250
* AssemblyHaplotypesMultiThreadPlugin.minInversionLen=7500
* AssemblyHaplotypesMultiThreadPlugin.loadDB=true
* AssemblyHaplotypesMultiThreadPlugin.assemblyMethod=mummer4
* AssemblyHaplotypesMultiThreadPlugin.numThreads=3
* AssemblyHaplotypesMultiThreadPlugin.entryPoint=all

The "entryPoint" parameter is provided for users who wish to speed up processing by running the mummer4 alignments on multiple machines. You may do this, then gather all the mummer files to a single machine.  If you do this, set the loadDB flag to FALSE, run all your alignments, then gather them to a single machine for loading into the database.

Description of optional parameters:

* mummer4Path: If the mummer4 executables exist in a path other than /mummer/bin, then provide that path via this parameter.  If you are running in the PHG docker, the path will be /mummer/bin and this parameter is not necessary.
* clusterSize:  This is a parameter to mummer4's nucmer alignment program.  It sets the minimum length for a cluster of matches.  We have found the value of 250 provides good coverage with acceptable speed.  If you wish a different value, set this parameter.
* minInversionLen: Minimum length of inversion for it to be kept as part of the alignment. Default is 7500
* Boolean: true means load haplotypes to db, false means do not populate the database.  Useful if the user is running mummer4 alignments on multiple machines and plans later to put all mummer4 files on 1 machine to load to a single database.  Defaults to true
* numThreads: Number of threads used to process assembly chromosomes.  The code will subtract 2 from this number to get the number of worker threads.  It leaves 1 thread for IO to the DB and 1 thread for the Operating System.
* entryPoint:  This parameter indicates at which point assembly processing should begin.  If a step other than "all" is chosen, files normally created from the previous steps must be present in a sub-directory named "align" in your output directory.  They must be named using the format shown below for the software to recognize them:
    * all: Run the entire pipeline.  The software creates all necessary files.
    * deltafilter:  Assumes the nucmer aligning step has already been performed. Processing starts with mummer4 delta-filter step includes running mummer4 show-coords on both the original alignment and the filtered alignment as well as mummer4 show-snps.  The delta file must be available in a sub-directory named "align" in the directory mounted to docker's /phg/outputDir/.  The file name must be in the format:
        * ref_<assemblyDBName from config file>_<chromosome from config file>_c<clusterSize>.delta
        * ex:  ref_W22_1_c250.delta
    * refilter: Assumes delta-filter and show-coords for both the delta and the delta-filter steps have been run. This step runs additional filtering to process overlapping alignments, remove embedded alignments and find a longest-path solution.Your output/align directory  must have the files for deltafilter plus:
        * ref_<assemblyDBName from config file>_<chromosome from config file>_c<clusterSize>.delta_filtered
        * ref_<assemblyDBName from config file>_<chromosome from config file>_c<clusterSize>.coords_orig
        * ref_<assemblyDBName from config file>_<chromosome from config file>_c<clusterSize>.coords_filtered
    haplotypes: This step takes the mummer4 output files and from them creates haplotype sequence fro the database.  If starting from this step, you must have the files from above plus:
        * ref_<assemblyDBName from config file>_<chromosome from config file>_c<clusterSize>.coords_filteredPlus_noEmbedded
        * ref_<assemblyDBName from config file>_<chromosome from config file>_c<clusterSize>.coords_final
        * ref_<assemblyDBName from config file>_<chromosome from config file>_c<clusterSize>.snps_prefiltered
        * ref_<assemblyDBName from config file>_<chromosome from config file>_c<clusterSize>.snps_final