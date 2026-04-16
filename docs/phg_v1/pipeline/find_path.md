!!! warning "Legacy Documentation - PHG Version 1"

    This section contains documentation for **PHG version 1**, which is
    no longer actively developed. It is preserved here for archival and
    historical reference only. If you are looking to use the Practical
    Haplotype Graph, please refer to the [PHG v2 documentation](../../index.md),
    which reflects the current version of the software.

The FindPath.sh script is used to determine the best path through the haplotype graph.  The script uses fastq files to generate haplotype counts using a Haplotype Graph created from the database. The counts are then used with a Hidden Markov Model (HMM) algorithm to find the best path through the graph for each taxon. The list haplotype node ids comprising the path is written to a text file.  The haplotype counts are written to the database table haplotype_counts, the paths data is written to the database table paths.

This shell script takes 8 parameters in this order. These parameters are passed to the called plugins.

* name of the PHG database from which sequences will be pulled
* the name of a database configuration file, for connecting to the database in later steps. A sample config file can be found here:[Config File Wiki](https://bitbucket.org/bucklerlab/practicalhaplotypegraph/wiki/DockerPipeline/ConfigFile).
* The name of the method (consensus type or other) used to create the haplotype sequences to be used when processing the skim sequences.
* The name of a reference fasta file used when processing the db haplotype sequences. (i.e. Zea_mays.AGPv4.dna.toplevel.fa)
* A method name to store for this haplotype counting instance.
* A method name to store for the paths table for this execution instance.
* The name of a file containing the reference ranges to keep (this is optional)

When invoked, the script first looks to see if a haplotype fasta file of the given name (first parameter) exists.  If not, one is created from the graph and indexed (bwa indexed).  The [FastqToHapCountPlugin](fastq_to_hap_count_plugin.md) is then run for every file in the fastq directory. This plugin determines the number of reads from the fastq file that match the sequence of the haplotype nodes. It scores which haplotypes are identical, excluded, or unresolved relative to a perfect hit GenotypeMap.  For every taxon represented by a fastq file, the haplotype counts are stored to the database and written to files in the fastq_hap_count directory.

Once the haplotype counts have been calculated, [HapCountBestPathToTextPlugin](hap_count_best_path_to_text_plugin.md) is run to choose the best path through the graph for each taxon.  The best paths are written to files in the hap_count_best_path directory.

When FindPath.sh is run as part of a Docker container script, the Docker script expects the following directory mount points:

* Mount localMachine:/pathToInputs/FastQFiles/ to docker:/tempFileDir/data/fastq
* Mount localMachine:/pathToInputs/refFolder/ to docker:/tempFileDir/data/reference/
* Mount localMachine:/pathToOutputs/ to docker:/tempFileDir/outputDir/
* Mount localMachine:/pathToInputs/config.txt to docker:/tempFileDir/data/config.txt

It is expected the database is stored in the User specified outputDir that is mounted below and the config.txt specifies the database name and login parameters.

An example Docker script to run the FindPath.sh shell script is:

```
#!python

docker1 run --name cbsu_phg_container_findPath --rm \
        -v /workdir/user/DockerTuningTests/InputFiles/Reference/:/tempFileDir/data/reference/ \
        -v /workdir/user/DockerTuningTests/InputFiles/GBSFastq/:/tempFileDir/data/fastq/ \
        -v /workdir/user/DockerTuningTests/DockerOutput/:/tempFileDir/outputDir/ \
        -v /workdir/user/DockerTuningTests/InputFiles/config.txt:/tempFileDir/data/config.txt \
        -t phgrepository_test:latest \
        /FindPath.sh panGenome12Taxa config.txt CONSENSUS Zea_mays.AGPv4.dna.toplevelMtPtv3.fa HAP_COUNT_METHOD  PATH_METHOD
```
The --name parameter provides a name for the container.  This is optional.

The --rm parameter indicates the container should be deleted when the program finishes executing.  This is optional.

The -v directives are used to mount data from the user machine into the Docker.  The path preceding the ":" is the path on the user machine.  The directory path following the ":" are the paths inside the Docker where the user home directories will be mounted.

The "-t" directive indicates the Docker image of which this container will be an instance.  The last line tells the Docker container to run the FindPath.sh script which is found in the root directory.  The items following are the parameters to the FindPath.sh script.