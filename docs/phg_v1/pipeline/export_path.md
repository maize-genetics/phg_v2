!!! warning "Legacy Documentation - PHG Version 1"

    This section contains documentation for **PHG version 1**, which is
    no longer actively developed. It is preserved here for archival and
    historical reference only. If you are looking to use the Practical
    Haplotype Graph, please refer to the [PHG v2 documentation](../../index.md),
    which reflects the current version of the software.

This script allows paths that have been exported to text files from the [HapCountBestPathToTextPlugin](hap_count_best_path_to_text_plugin.md) to be imported via [ImportHaplotypePathFilePlugin](import_haplotype_path_file_plugin.md) and associated with nodes in the haplotype graph, then exported as a VCF file via the [PathsToVCFPlugin](paths_to_vcf_plugin.md).  This method needs to be updated to take paths either from files or from the PHG DB paths table.  The end result of this script is a VCF file.

The parameters to this shell script are:

* Data base config file:  Configuration file containing properties host, user, password, DB and DBtype where DBtype is either sqlite or postgres. A sample config file can be found here:[Config File Wiki](https://bitbucket.org/bucklerlab/practicalhaplotypegraph/wiki/DockerPipeline/ConfigFile)
* Consensus method: The Consensus Method used to create haplotypes in graph.
* Output file: Output VCF file
* path method: Used to grab the correct paths from the DB.  This needs to match what was used in FindPathMinimap2.sh
 
When ExportPath.sh runs in a Docker container, the following mounted directories are expected

* Mount localMachine:/pathToOutputs/FindPathDir/ to docker:/tempFileDir/outputDir/. To make this work correctly, the DB must also be here.  If the DB is defined in a different place within the config file, you will need to make a new mount point accordingly.
* Mount localMachine:/pathToInputs/config.txt to docker:/tempFileDir/data/config.txt

In addition, it is expected the database is stored in the User FindPathDir that is mounted below

An example Docker script to run the ExportPath.sh shell script is:

```
#!python

docker run --name cbsu_phg_container_exportPath --rm \
        -v /workdir/user/DockerTuningTests/DockerOutput/FindPathDir/:/tempFileDir/outputDir/ \
        -v /workdir/user/DockerTuningTests/InputFiles/config.txt:/tempFileDir/data/config.txt \
        -t maizegenetics/phg:latest \
        /ExportPath.sh config.txt CONSENSUS testOutput1.vcf PATH_METHOD
```
The --name parameter provides a name for the container.  This is optional.

The --rm parameter indicates the container should be deleted when the program finishes executing.  This is optional.

The -v directives are used to mount data from the user machine into the Docker.  The path preceding the ":" is the path on the user machine.  The directory path following the ":" are the paths inside the Docker where the user home directories will be mounted.

The "-t" directive indicates the Docker image of which this container will be an instance.  The last line tells the Docker container to run the ExportPath.sh script which is found in the root directory.  The items following are the parameters to the ExportPath.sh script.