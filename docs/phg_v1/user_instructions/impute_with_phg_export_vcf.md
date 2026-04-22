!!! warning "Legacy Documentation - PHG Version 1"

    This section contains documentation for **PHG version 1**, which is
    no longer actively developed. It is preserved here for archival and
    historical reference only. If you are looking to use the Practical
    Haplotype Graph, please refer to the [PHG v2 documentation](../../index.md),
    which reflects the current version of the software.

# Export variants as VCF

## Quick Start

1. Change `MethodName`, `Database`,  and `output.vcf` in your config file to match file paths on your computer.
2. Run `phg pathsToVCF [config.txt]`

## Details

This Process exports a VCF corresponding to haplotype calls from the mapping and path plguins. It looks for the indicated path methodName and export variants for all taxa form this method. The name of the output VCF must be indicated in the config file.  

## Kitchen Sink

The PathsToVCFPlugin is called after paths have been found and stored in the PHG database ([see imputation pipeline](impute_with_phg_main.md)) . It is used to create VCF files from haplotype paths. The user may specify certain reference ranges to be included in the exported file. Paths will have been added to the database and/or exported to text files with the previous step, which ran the FastqToMappingPlugin and the BestHaplotypePathPlugin. 

If haplotypes were not stored in the database and were only output to files, this step also requires that the ImportHaplotypePathFilePlugin is called before the PathsToVCFPlugin. 

The output from this step is a VCF file.

There are 4 parameters used in this step:

* CONFIG_FILE: Configuration file containing properties host, user, password, DB and DBtype where DBtype is either sqlite or postgres. A sample config file can be found here:Master Config File
* CONSENSUS_METHOD: The Consensus Method used to create haplotypes in graph.
* OUTPUT_FILE: Output VCF file
* PATH_METHOD_NAME: Used to grab the correct paths from the DB.  This needs to match what was used in FindPathMinimap2.sh

### *Details on running this step with wrapper scripts*

When running this step on the command line, all file paths and parameters are set in the config file. The only call that needs to be run in the terminal is `phg pathsToVCF /path/to/config.txt`. If you would like to overwrite the parameters set in the config file, you can do that by setting the parameters on the command line directly.

For example, to ignore the config file PATH_METHOD_NAME level and set one directly, you could run:
```
phg findPaths -configFile /path/to/config.txt -PATH_METHOD_NAME MyNewPath1
```

### *Details on running this step through docker*

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

The -t directive indicates the Docker image of which this container will be an instance.  The last line tells the Docker container to run the ExportPath.sh script which is found in the root directory.  The items following are the parameters to the ExportPath.sh script.

### *Files*

**Config file**

An example can be found here: Master config file

**Haplotype Path files** (optional)

By default, haplotype paths are added to the database. However, in the previous step you may have chosen to output haplotype paths as files instead of adding them to the database. If you output haplotype files in the previous step, you will have a set of files with one haplotype ID per line. Each file corresponds to the predicted path for one taxon.

### *Plugins*

#### **PathsToVCFPlugin**

The PathsToVCFPlugin is used to create VCF files from haplotype paths.  The user may specific specific reference ranges to be included in the exported file.  This plugin can be chained with the HaplotypeGraphBuilderPlugin and the ImportHaplotypePathFilePlugin. It takes input from both and uses this to create the requested VCF files.

The parameters to this plugin are:

* -outputFile <Output file name> Output VCF File Name. (Default=null)(OPTIONAL)
* -refRangeFileVCF <Reference Range File> Reference Range file used to further subset the paths for only specified regions of the genome. (OPTIONAL)
* -positions <Position List> Genotype file (i.e. VCF, Hapmap, etc.), bed file, or json file containing the requested positions. (OPTIONAL)
* -ref <Reference Genome> Reference Genome. (OPTIONAL)
* -outputAllSNPs <true | false> Whether to output all SNPs known by haplotype graph.

Here is an example of how to chain the HaplotypeGraphBuilderPlugin to the PathsToVCFPlugin.
```
perl /tassel-5-standalone/run_pipeline.pl $fullXmx -debug \
-configParameters ${DATABASE_CONFIG} \
-HaplotypeGraphBuilderPlugin \
	-configFile $DATABASE_CONFIG \
	-methods $CONSENSUS_METHOD \
	-includeVariantContexts true \
	-endPlugin \
-PathsToVCFPlugin \
	-outputFile $OUTPUT_FILE \
	-endPlugin
```

#### **ImportHaplotypePathFilePlugin**

Path information is stored to the PHG database by default, so you are unlikely to need to run this. However, if you chose in the previous step to output a file with haplotype paths and not store paths to the database, you will need to run the ImportHaplotypePathFilePlugin before the PathsToVCFPlugin. 

This plugin expects a PHG Haplotype Graph to be sent as an incoming parameter, which can be aquired by chaining the HaplotypeGraphBuilderPlugin to the ImportHaplotypePathFilePluign. The output from the ImportHaplotypePathFilePlugin is the input to the PathsToVCFPlugin.

The ImportHaplotypePathFilePlugin takes the following parameters;

* -inputFileDirectory <Input File Directory> A directory containing the paths text files to be imported and used to create the VCF files. (Default=null) (REQUIRED)


Here is an example of how to chain all three of these plugins to read in haplotype path files and export a VCF.
```
perl /tassel-5-standalone/run_pipeline.pl $fullXmx -debug \
-configParameters ${DATABASE_CONFIG} \
-HaplotypeGraphBuilderPlugin \
	-configFile $DATABASE_CONFIG \
	-methods $CONSENSUS_METHOD \
	-includeVariantContexts true \
	-endPlugin \
-ImportHaplotypePathFilePlugin \
	-pathMethodName ${PATH_METHOD_NAME} \
	-endPlugin \
-PathsToVCFPlugin \
	-outputFile $OUTPUT_FILE \
	-endPlugin
```

## Troubleshooting



[Return to Step 3 pipeline](impute_with_phg_main.md)

[Return to Wiki Home](../home.md)
