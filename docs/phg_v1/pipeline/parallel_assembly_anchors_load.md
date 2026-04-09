!!! warning "Legacy Documentation - PHG Version 1"

    This section contains documentation for **PHG version 1**, which is
    no longer actively developed. It is preserved here for archival and
    historical reference only. If you are looking to use the Practical
    Haplotype Graph, please refer to the [PHG v2 documentation](../../index.md),
    which reflects the current version of the software.

The ParallelAssemblyAnchorsLoad.sh script takes 3 parameters:  a configuration file containing database connection information, a keyfile containing specifics on reference and assembly fastas files to be aligned and turned into haplotypes for loading into a PHG database, and an output directory.  

This script calls the AssemblyHaplotypesMultiThreadPlugin plugin, which takes reference and assembly fastas, split by chromosome, and aligns them with mummer4 scripts.  It performs the same functions as are executed
when the LoadAssemblyAnchors.sh script is called, but will run the mummer4 alignment programs in parallel.  The number of processes run concurrently is dependent on a user supplied thread count.

When determining the number of threads to run, please keep in mind the alignment of each ref-assembly
fasta pair will take RAM commensurate with the size of the genomes.  Consider the RAM limitations
of your machine as well as the number of cores available.

The required input to this script is a configuration file with database connection information,
a keyfile with fasta information and an output directory.  Both the configuration file and the keyfile must live in a local machine directory that is mounted to the docker's /tempFileDir/data directory.

When ParallelAssemblyAnchorsLoad.sh is called as part of a Docker script, the Docker script expects these mount points:

* Mount localMachine:/pathToDataDir/ to docker:/tempFileDir/data/
* Mount localMachine:/pathToOutputDir/ to docker:/tempFileDir/outputDir/

The script will create a subdirectory called "align" under the mounted output directory where mummer4 alignment files will be stored.

Description of required parameters:

* keyfile: The keyfile must be a tab-delimited file with 1 row for each per-chromosome assembly fasta you wish to align and load to the database.  The required columns for the keyfile are RefDir, RefFasta,  AssemblyDir, AssemblyFasta, AssemblyDBName and Chromosome.  The oolumn descriptions are below.
    * RefDir: full path to the directory containing the reference fasta file.
    * RefFasta: name of the reference fasta file to be aligned to the assembly fasta file in this row.  This should be a fasta file at the chromosome level.
    * AssemblyDir: full path to the directory containing the assembly fasta files.
    * AssemblyFasta:  name of the assembly fasta file to be aligned to the reference fasta file in this row.  This should be a fasta file at the chromosome level.
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

The "entryPoint" parameter is provided for users who wish To speed up processing by running the mummer4 alignments on multiple machines. You may do this, then gather all the mummer files to a single machine.  If you do this, set the loadDB flag to FALSE, run all your alignments, then gather them to a single machine for loading into the database.

Description of optional parameters:

* mummer4Path: If the mummer4 executables exist in a path other than /mummer/bin, then provide that path via this parameter.  If you are running in the PHG docker, the path will be /mummer/bin and this parameter is not necessary.
* clusterSize:  This is a parameter to mummer4's nucmer alignment program.  It sets the minimum length for a cluster of matches.  We have found the value of 250 provides good coverage with acceptable speed.If you wish a different value, set this parameter.
* minInversionLen: Minimum length of inversion for it to be kept as part of the alignment. Default is 7500
* Boolean: true means load haplotypes to db, false means do not populate the database.  Useful if the user is running mummer4 alignments on multiple machines and plans later to put all mummer4 files on 1 machine to load to a single database.  Defaults to true
* numThreads: Number of threads used to process assembly chromosomes.  The code will subtract 2 from this number to get the number of worker threads.  It leaves 1 thread for IO to the DB and 1 thread for the Operating System.
* entryPoint:  This parameter indicates at which point assembly processing should begin.  If a step other than "all" is chosen, files normally created from the previous steps must be present in a sub-directory named "align" in your output directory.  They must be named using the format shown below for the software to recognize them:
    * all: Run the entire pipeline.  The software creates all necessary files.
    * deltafilter:  Assumes the nucmer aligning step has already been performed. Processing starts with mummer4 delta-filter step includes running mummer4 show-coords on both the original alignment and the filtered alignment as well as mummer4 show-snps.  The delta file must be available in a sub-directory named "align" in the directory mounted to docker's /tempFileDir/outputDir/.  The file name must be in the format:
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

An example Docker script that would call this shell script is below.  Note: you need to change the mounted user directories to match your own directory configuration.  The docker directories (right side of the count pair beginning with "tempFileDir") need to remain as shown here.  The file names, e.g. "config.txt" should match the names of your files.  And the "config.txt" and "parallelKeyFile_B204.txt" must live in the directory mounted to /tempFileDir/data/.

```
echo "Starting taxa B104"

docker1 run --name phg_assembly_container_B104 --rm \
        -v /workdir/${USER}/phg_assemblyParallel_test/DockerOutput/:/tempFileDir/outputDir/ \
        -v /workdir/${USER}/phg_assemblyParallel_test/DataFolders/:/tempFileDir/data/ \
        -v /workdir/${USER}/phg_nam_assemblies/B73/:/tempFileDir/referenceFastas/ \
        -v /workdir/${USER}/phg_nam_assemblies/B104/:/tempFileDir/assemblyFastas/ \
        -t phgrepository_test:latest \
        /ParallelAssemblyAnchorsLoad.sh config.txt \
                parallelKeyFile_B104.txt


echo "Finished taxa 104 "


```

An example key file the has information for a single assembly genome is shown below.  Note all fields are tab-delimited and because this is run in a docker, the reference and assembly directories are docker paths.

```
RefDir  RefFasta        AssemblyDir     AssemblyFasta   AssemblyDBName  Chromosome
/tempFileDir/referenceFastas/   B73chr1.fa      /tempFileDir/assemblyFastas/    B104chr1.fa     B104_Assembly   1
/tempFileDir/referenceFastas/   B73chr2.fa      /tempFileDir/assemblyFastas/    B104chr2.fa     B104_Assembly   2
/tempFileDir/referenceFastas/   B73chr3.fa      /tempFileDir/assemblyFastas/    B104chr3.fa     B104_Assembly   3

```

An example key file that has per-chromosome pastas for multiple species is shown in the next example.  Note that the assembly chromosome fasta files must exist in the same directory.

```
RefDir  RefFasta        AssemblyDir     AssemblyFasta   AssemblyDBName  Chromosome
/tempFileDir/referenceFastas/   B73chr1.fa      /tempFileDir/assemblyFastas/    B104chr1.fa     B104_Assembly   1
/tempFileDir/referenceFastas/   B73chr2.fa      /tempFileDir/assemblyFastas/    B104chr2.fa     B104_Assembly   2
/tempFileDir/referenceFastas/   B73chr3.fa      /tempFileDir/assemblyFastas/    B104chr3.fa     B104_Assembly   3
/tempFileDir/referenceFastas/   B73chr1.fa      /tempFileDir/assemblyFastas/    B97chr1.fa      B97_Assembly    1
/tempFileDir/referenceFastas/   B73chr2.fa      /tempFileDir/assemblyFastas/    B97chr2.fa      B97_Assembly    2
/tempFileDir/referenceFastas/   B73chr3.fa      /tempFileDir/assemblyFastas/    B97chr3.fa      B97_Assembly    3
```


The ParalleleAssemblyAnchorsLoad.sh script uses the parameter cache option to load the configuration
file into memory.  This is the -configParameters <configFile> directive given to the tassel run_pipeline.pl script used in the Docker container. if you run this plugin manually at the command line, you will need to include this yourself with a command similar to the one below:

```
/tassel-5-standalone/run_pipeline.pl -Xmx100G -debug -configParameters myConfigFile.txt -AssemblyHaplotypesMultiThreadPlugin  -keyFile myKeyFile.txt  -outputDir /workdir/user/phg_assemblies/DockerOutput/ -endPlugin
```