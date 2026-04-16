!!! warning "Legacy Documentation - PHG Version 1"

    This section contains documentation for **PHG version 1**, which is
    no longer actively developed. It is preserved here for archival and
    historical reference only. If you are looking to use the Practical
    Haplotype Graph, please refer to the [PHG v2 documentation](../../index.md),
    which reflects the current version of the software.

The CreateConsensi.sh script is used to create Consensus Haplotypes for each genome interval using a specific set of raw haplotypes.  Basically this script will create a haplotype graph and attempt to create consensus haplotypes for each anchor.  Once a set of consensus haplotypes are created, the script will upload the all consensus haplotypes to the Database.  Once this process is done, a graph can be made using the consensus haplotypes.

This processing is done via the [RunHapCollapsePlugin](run_hap_collapse_plugin.md).

# Mount Points For use with the PHG Docker #
* Mount localMachine:/pathToInputs/config.txt:/tempFileDir/data/config.txt - (required)
* Mount localMachine:/pathToInputs/refFolder/ to docker:/tempFileDir/data/reference/ - (required)
* Mount localMachine:/pathsToDB/dbFolder/ to docker /tempFileDir/outputDir/ - (required)
### There is one optional mount point for getting out the VCF files created from combining the GVCFs ###
* This is used in conjunction with the exportMergedVCF config file parameter.  If this parameter is set in the config file, you can mount a local directory to that docker directory to get out the intermediate VCFs.  If this is not set, the code should run faster as there is less IO. 

# CreateConsensi.sh Parameters #
* configFile: Path to config file containing DB parameters host, user, password, DB, type. Used for making the database connection. Type must be wither "sqlite" or "postgres" to identify db type for connection. This file will also contain parameters related to creating the consensus sequences and is described in Relevant Config File Parameters.  In the below example this is /tempFileDir/data/config.txt.  A sample configuration file can be found here: [Config File Wiki](https://bitbucket.org/bucklerlab/practicalhaplotypegraph/wiki/DockerPipeline/ConfigFile)
* Reference file: Reference fasta file.  In the below example this is reference.fa and is contained in the mount point /tempFileDir/data/reference/
* Haplotype Method Name: Name of the raw haplotypes which Consensus should be made.
* Consensus Method Name: Name to be stored in the database denoting this run of CreateConsensi.

# Pipeline Steps #
1. Create a Haplotype Graph
2. For Each Anchor do the following:
    1. If the graph does not have the variants, extract the variants from the DB for the raw Haplotypes.
    2. Take all the GVCF variants and combine them into a single TASSEL GenotypeTable object.
    3. Remove any Positions covering an indel.
    4. Create Consensus clusters and create GenotypeTables for each cluster.
    5. Reintroduce indels which were filtered out if they agree.  If they do not agree, it will follow the indelMergeRule defined in the Configuration File.
3. Collect up to 1000 Consensus haplotypes and upload them to the DB.

# Relevant Config File Parameters.  Must be in the form param=value #
Sample file is here: [Config File Wiki](https://bitbucket.org/bucklerlab/practicalhaplotypegraph/wiki/DockerPipeline/ConfigFile)
### DB Parameters ###
* host - host name of the db
* user - db user name
* password - password for db
* DB - name of the db
* DBtype - sqlite or postgres

### Java Arguments ###
* Xmx - max heap space for the pipeline.  This is similar to Java's -Xmx argument.

### Graph Building Parameters ###
* includeVariants - true. This needs to be true to work correctly.  Given the new Variant Storage, the graphs variants can fit into memory much better.


### Consensus Parameters ###
* mxDiv - maximum amount of divergence allowed when clustering.  The default is 0.01. 
* seqErr - Error rate for calling a het vs a homozygote in the consensus haplotype.  The default is 0.01.
* minSites - minimum number of non-missing sites to be clustered.  Any haplotypes with fewer number of non-missing calls will be ignored.  The default is 20.
* minTaxa -  minimum number of taxa for a cluster to be created. The default is 2.
* rankingFile - path to a ranking file.  If your haplotypes were created from Assemblies, you need to produce a method ranking the taxon. This file will take the form: taxon\tscore where higher scores mean we trust that taxon more highly.  Do no include a header line.  When clustering assemblies, when we have a cluster of similar haplotypes, we choose whichever taxon in that group which has the higher ranking score.  To break ties, be sure to give each taxon a different score.  One simple way to score things is to count the number of haplotypes covered by each taxon in the DB and use that count as the score.  Any other arbitrary ranking can be used.
* clusteringMode - either upgma(default) or upgma_assembly.  If running WGS haplotypes, either use the default or specify upgma.  If the haplotypes are assemblies, use upgma_assembly.  The differences between the two are that upgma just builds a tree based on a pairwise calculated distance matrix and then will try to merge haplotypes together.  upgma_assembly will also do a pairwise distance matrix and then will select the better haplotype in the group based on the ranking file specified by the "rankingFile" parameter.

# Example Run Scripts #
An example Docker script to run the CreateConsensi.sh shell script is:

```
#!bash

docker run --name phg_container_consensus --rm \
        -v /workdir/user/DockerTuningTests/Reference/:/tempFileDir/data/reference/ \
        -v /workdir/user/DockerTuningTests/DockerOutput/:/tempFileDir/outputDir/ \
        -v /workdir/user/DockerTuningTests/DataFolders/LoadRefDataDocker/config.txt:/tempFileDir/data/config.txt \
        -t maizegenetics/phg:latest \
        /CreateConsensi.sh /tempFileDir/data/config.txt reference.fa GATK_PIPELINE CONSENSUS
```

The --name parameter provides a name for the container.  This is optional.

The --rm parameter indicates the container should be deleted when the program finishes executing.  This is optional.

The -v directives are used to mount data from the user machine into the Docker.  The path preceding the ":" is the path on the user machine.  The directory path following the ":" are the paths inside the Docker where the user home directories will be mounted.

The "-t" directive indicates the Docker image of which this container will be an instance.  The last line tells the Docker container to run the CreateConsensi.sh script which is found in the root directory.  The items following are the parameters to the CreateConsensi.sh script.