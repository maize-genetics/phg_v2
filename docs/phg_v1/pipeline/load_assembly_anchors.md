!!! warning "Legacy Documentation - PHG Version 1"

    This section contains documentation for **PHG version 1**, which is
    no longer actively developed. It is preserved here for archival and
    historical reference only. If you are looking to use the Practical
    Haplotype Graph, please refer to the [PHG v2 documentation](../../index.md),
    which reflects the current version of the software.

The LoadAssemblyAnchors.sh script uses mummer4 scripts to align a single chromosome reference fasta file to a single matching chromosome assembly fasta file.  For example:  run maize b73 chromosome 1 against maize ph207 chromosome 1, or b73 chromosome 2 against ph207 chromosome 2. The Mummer4 nucmer, delta-filter, show-coords and show-snps scripts are run to create the necessary data.  This data is then processed and put in a format necessary for loading to the PHG data base.  View details of this process in the description of the [AssemblyHaplotypesPlugin](assembly_haplotypes_plugin.md), which performs this functionality.

When the files are split by chromosome, the fasta id line should contain the chromosome name.  If additional data is desired on this line, there must be a space after the chromosome name and before the additional data.  The chromosome names should either be just a number, or an alpha numeric string.  Note that any leading "chr" or "chromosome" will be stripped off.  This could cause inconsistencies in names between the name stored in the database and the name in the fasta files.  It is recommended you do NOT precede your chromosome name with "chr" or "chromosome"


The reference and assembly chromosome names must be consistent in the respective fasta files.

When LoadAssemblyAnchors is called as part of a Docker script, the Docker script expects these mount points:

* Mount localMachine:/pathToAssemblyFastas/ to docker:/tempFileDir/data/assemblyFasta
* Mount localMachine:/pathToDataDir/ to docker:/tempFileDir/data/
* Mount localMachine:/pathToOutputDir/ to docker:/tempFileDir/outputDir/
* Mount localMachine:/pathToReferenceFastaDir/ to docker:/tempFileDir/data/reference

Note the script will create a subdirectory named "align" under the /tempFileDir/outputDir/ directory. This is where the mummer4 programs output will be written.

The parameters expected by this shell script are:

* DB config File: DB Config File containing properties host,user,password,DB and DBtype where DBtype is either "sqlite" or "postgres". This file is expected to live in the folder mounted to Docker tempFileDir/data/. A sample config file can be found here:[Config File Wiki](https://bitbucket.org/bucklerlab/practicalhaplotypegraph/wiki/DockerPipeline/ConfigFile)
* Reference Fasta File:  The reference fasta file with only 1 chromosome of data.  This file is used when running Mummer4 alignment scripts.  This file is expected to live in the folder mounted to Docker /tempFileDir/data/reference.
* Assembly Fasta File:  The assembly fasta file with only 1 chromosome of data.  This file is used to align against the respective reference chromosome when the mummer4 scripts are run.  It is expected this file is in the folder mounted to Docker path /tempFileDir/data/assemblyFasta.
* Assembly name:  name to append to the output files that identifies the assembly processed. This is also the name that will be stored in the database genotypes table.  To differentiate the assembly from other instances of taxa with this name, you may want to name it <taxa>_Assembly.
* chrom:  the chromosome that is processed. This must match the chromosome name as stored on the ID line in the fasta file, or can be just the chromosome number if the chromosome is listed as "chr1" or "chromosome1".  Note if leading 0's are added to the chromosome name, they must be consistently added in both assembly fasta, reference fasta, and this shell parameter.
* clusterSize:  The minimum length for a cluster of matches used when running mummer4 nucmer program.  Nucmer default is 65 but we have found when running maize assemblies 250 is faster and provides good results.  That is the default we show below.

All output from running this script will be written to files in the directory mounted to Docker /tempFileDir/outputDir/align.

The assembly entries in the genotypes table will be defaulted to the following values:

* line_name:  the assembly name
* line_data: the mummer script parameters used with this version of PHG
* Ploidy:  1
* is_reference:  false
* isPhasedAcrossGenes: true
* isPhasesAcrossChromosomes: true

The assembly entries in the gametes table will be defaulted to the following values:

* hapNumber: 0
* Confidence:  1

The assembly entries in the methods table will be defaulted to the following values:

* method_type: DBLoadingUtils.MethodTYpe.ASSEMBLY_HAPLOTYPES
* name:  mummer4
* description:  the mummer script parameters used with this version of PHG ( same as genotypes:line_data

An example Docker script that would call this shell script is below.  Note: you need to change the mounted user directories to match your own directory configuration.  The docker directories (right side of the count pair beginning with "tempFileDir" need to remain as shown here.  The file names, e.g. "config.txt" should match the names of your files.

This example runs through all 10 chromosomes for a particular assembly.  In this example, the assembly files are broken into fasta files name ph207chr1.fasta, ph297chr2.fasta, etc.  The reference fastas are broken by chromosome and are named b73chr1.fasta, b73chr2.fasta, etc.  You need to change those names to match your fasta file names.

```
#!python

#!/bin/bash
chromList=(1 2 3 4 5 6 7 8 9 10)

for chrom in "${chromList[@]}"
do

echo "Starting chrom ${chrom} "

docker run --name phg_assembly_container_chr${chrom} --rm \
        -v /workdir/lcj34/mummer4_testing/output:/tempFileDir/outputDir/ \
        -v /workdir/lcj34/mummer4_testing/data:/tempFileDir/data/ \
        -v /workdir/lcj34/mummer4_testing/refFastas/:/tempFileDir/data/reference/ \
        -v /workdir/lcj34/mummer4_testing/assemblyFastas/:/tempFileDir/data/assemblyFasta/ \
        -v /workdir/lcj34/mummer4_testing/mummer_from_shellDocker/:/tempFileDir/outputDir/align/ \
        -t phgrepository_test:latest \
        /LoadAssemblyAnchors.sh configSQLite_docker.txt \
                b73chr${chrom}.fasta \
                ph207chr${chrom}.fasta \
                PH207_Assembly \
                ${chrom} \
                250


echo "Finished chrom  ${chrom} "
done

```

Note: use docker1 if running on cbsu.

The --name parameter provides a name for the container.  This is optional.

The --rm parameter indicates the container should be deleted when the program finishes executing.  This is optional.

The -v directives are used to mount data from the user machine into the Docker.  The path preceding the ":" is the path on the user machine.  The directory path following the ":" are the paths inside the Docker where the user home directories will be mounted.

The "-t" directive indicates the Docker image of which this container will be an instance.  The last line tells the Docker container to run the LoadAssemblyAnchors.sh script which is found in the root directory.  The items following are the parameters to the LoadAssemblyAnchors.sh script.