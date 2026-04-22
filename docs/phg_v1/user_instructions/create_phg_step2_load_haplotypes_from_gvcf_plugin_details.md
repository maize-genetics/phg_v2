!!! warning "Legacy Documentation - PHG Version 1"

    This section contains documentation for **PHG version 1**, which is
    no longer actively developed. It is preserved here for archival and
    historical reference only. If you are looking to use the Practical
    Haplotype Graph, please refer to the [PHG v2 documentation](../../index.md),
    which reflects the current version of the software.

# LoadHaplotypesFromGVCFPlugin Details

When you have a GVCF file (created from the MAFToGVCFPlugin or other means), the next step is to load
the database with the information from this GVCF.  Please note you need to index the GVCF files.  The MAFToGVCFPlugin code will
do this for you.  If you have created the GVCF files from other means, you can run htsjdk's tabix on the files to create the indexes.

The LoadHapltoypesFromGVCFPlugin takes a keyfile containing
a list of GVCF files for uploading.  It parses and formats the data, then loads this data as haplotypes to a PHG database defined by the user configuration file.

This step can be run using a shell script similar to the one detailed below:


```
WORKING_DIR=/workdir/lcj34/anchorwaveTesting/
DOCKER_CONFIG_FILE=/phg/config.txt

# The LoadHaplotypesFromGVCFPlugin parameters will be obtained via the config file
docker1 run --name load_haplotypes_container --rm \
    -v ${WORKING_DIR}/:/phg/ \
    -t maizegenetics/phg:latest \
    /tassel-5-standalone/run_pipeline.pl -Xmx100G -debug -configParameters ${DOCKER_CONFIG_FILE} \
    -LoadHaplotypesFromGVCFPlugin -endPlugin
```


It may also be run directly from tassel-5-standalone as per this command.  Parameters are explicitly defined in this example, but they may alternately be pulled from a config file as with the example above.

```
/workdir/lcj34/tassel-5-standalone/run_pipeline.pl -Xmx100G -debug -LoadHaplotypesFromGVCFPlugin -wgsKeyFile /pathToDir/keyFile.txt -gvcfDir /workdir/lcj34/GVCF_files/ -referenceFasta ref.fasta -bedFile /pathToDir/mybedfile.bed -hapltoypteMethodName gvcfHaplotypes > output_loadHaplotypeFromGVCF.txt
```


If running with a config file, the parameters listed below should be defined(replace values with appropriate values for your own run).  Optional parameters are not shown in the example below, but may be included as desired.

```
LoadHaplotypesFromGVCFPlugin.referenceFasta=/phg/inputDir/reference/Zm-B73-REFERENCE-NAM-5.0.fa
LoadHaplotypesFromGVCFPlugin.wgsKeyFile=/phg/gvcfKeyFile.txt
LoadHaplotypesFromGVCFPlugin.gvcfDir=/phg/inputDir/loadDB/gvcfs/
LoadHaplotypesFromGVCFPlugin.bedFile=myBedFile.bed
LoadHaplotypesFromGVCFPlugin.haplotypeMethodName=assembly_by_anchorwave
LoadHaplotypesFromGVCFPlugin.haplotypeMethodDescription="files aligned with anchorwave, then turned to gvcf with plugin"

```


# Details

The GVCF files that will be loaded are defined in a keyfile that provides the location of the gvcfs, the sample names to use with each file, and information relating to gene phasing for the sample.  All data will be processed into haplotypes to be stored in the PHG database haplotypes table.  These haplotypes will be associated with the method defined by the plugin parameters.

## Parameters

Run the following command to print a current list of parameters for this plugin:

```
docker run --rm  maizegenetics/phg /tassel-5-standalone/run_pipeline.pl -LoadHaplotypesFromGVCFPlugin
```

## Parameter Descriptions

* wgsKeyFile: (required) Name of the keyfile to process.  Details on the keyfile are in a separate section below.
* gvcfDir: (required) Full path to the directory holding all the gvcf files.  If running with PHG Version 0.0.40 or lower, these are files with *.gvcf extension.  If running with PHG version 1.0 or higher, these are bgzipped and tabix'd versions of the gvcf files.  For each genome, there should be a *.gz and corresponding *.gz.tbi file.
* referenceFasta: (required) Full path to the reference fasta file
* sampleName: (required ) Sample name to write to the GVCF file.  This is also the name that will later be written to the database.
* bedFile: (required) Path to a bed formatted file containing the database reference range information.
* haplotypeMethodName: (required) The method name to use for the haplotypes being uploaded. This is a user defined method name that should be descriptive of the method used to create the haplotypes.  Examples would be "anchorwave_NAM_assemblies" or "anchorwave_released_assemblies" or "WGS_haplotypes".  The name may be as long as you wish but should be meaningful to you.
* haplotypeMethodDescription: (optional: no default) Description of the method used for the haplotypes being uploaded. This may describe methods of alignment or any other information specific to the haplotype gvcf creation.
* numThreads: (optional: default=3) The number of CPU threads to use when processing the data.  The GVCF upload will subtrace 2 from this number, elaving 1 thread for IO to the DB and 1 thread for the operating system.
* maxNumHapsStaged: (optional: default=1000) Number of haplotype instances to stage per chromosome before a DB write is triggered.  Lower this number if you are running into RAM issues.  It will take longer to process, but should help balance the load.
* mergeRefBlocks: (optional: default=false)  Merge consecutive GVCF ReferenceBlocks together.  If there is at least 1 bp between two gvcf refBlock records, the records will not be merged
* queueSize: (optional: default=30) Size of Queue used to pass information to the DB writing thread.  Increase this number to have better thread utilization at the expense of RAM. If you are running into Java heap Space/RAM issues and cannot use a bigger machine, decrease this parameter.
* isTestMethod: (optional: default=false) Indication if the data is to be loaded against a test method. Data loaded with test methods are not cached with the PHG ktor server

## Keyfile
The LoadHaplotypesFromGVCFPlugin keyfile is used to specify the inputs to the database.  It must have columns sample_name, files, type, chrPhased, genePhased, and phasingConf. If running with PHG version 1.0 or greater, it must also have a column named gvcfServerPath.  Optionally you can have sample_description, assemblyServerPath, assemblyLocalPath. Ordering of the columns does not matter.

If present, the "sample_description" column  becomes the description stored in the genotypes table "description" field associated with the sample taxon.  The "assemblyServerPath" and "assemblyLocalPath" columns are used when the gvcf files represent data from assemblies. The information in these 2 columns is used to populate the genome_file_data table.  The "assemblyServerPath" value is stored in this PHG table to preserve information on where this assembly may be found in the public domain.  The "assemblyLocalPath" value is used to access the genome fasta locally for creating a hash value to store with this file data.

The keyfile should be in tab-delimited format.  The headers should look as below: (this example includes the optional columns, as well as the gvcfServerPath required for PHG version 1.0 or greater)

```
type    sample_name     sample_description      mafFile files   chrPhased       genePhased      phasingConf     assemblyServerPath      assemblyLocalPath  gvcfServerPath
```

An example keyfile for this plugin is below. Note all fields are tab-delimited and because this is run in a docker, the assemblyLocalPath value is a docker path.  The assemblyServerPath should be an accessible location outside the docker, as should be the gvcfServerPath.  The gvcfServerPath values must be of the format <server>;<path on server>, with a semi-colon separating the server and path values

```
type    sample_name     sample_description      mafFile files   chrPhased       genePhased      phasingConf     assemblyServerPath      assemblyLocalPath       gvcfServerPath
GVCF    LineA_Assembly  smallSeq description for assembly LineA /phg/outputDir/align//LineA.maf LineA.gvcf.gz   true    true    0.9     irods:AssemblyFilePath/LineA.fa /phg/inputDir/assemblies//LineA.fa      198.42.78.6;/mypath/remoteGvcfs/
GVCF    LineB_Assembly  smallSeq description for assembly LineB /phg/outputDir/align//LineB.maf LineB.gvcf.gz   true    true    0.9     irods:AssemblyFilePath/LineB.fa /phg/inputDir/assemblies//LineB.fa      198.42.78.6;/mypath/remoteGvcfs/
```

Column details:

* sample_name: name of the taxon you are uploading
* files: a comma-separated list of file names(without path) providing the names of the GVCF files to be uploaded.  If there are multiple files, gamete_ids will be assigned in order.
* type: the type of files.   This plugin will only process keyfile lines with type = "GVCF".
* chrPhased: must be either 'true' or 'false' and represent whether or not the chromosome is phased.
* genePhased: must be either 'true' or 'false' and represent whether or not the genes are phased.
* phasingConf: must be a number between 0.0 and 1.0 representing the confidence in the phasing.
* sample_description:  a short description of the sample name to be uploaded.  If this column is not specified an empty description will be used.
* assemblyServerPath (optional column): an external path where the genome fasta file for this assembly may be found. It should be directory and file.  This column is optional and generally included when the GVCF files are the result of assembly aligning with the anchorwave program.
* assemblyLocalPath (optional column): the path to the assembly genome fasta file that can be accessed from this program.  This file is used to create am MD5 hash for the genome_file_data table.  This column is optional and generally included when the GVCF files are the result of assembly aligning with the anchorwave program.
* gvcfServerPath (required when running PHG versions 1.0 or greater) an external server and path to which this gvcf file will be uploaded.

[Return to Step 2 pipeline, version 0.0.40 or earlier](create_phg_step1_2_main.md)

[Return to Step 2 pipeline, version 1.0 or later](create_phg_step2_assembly_and_wgs_haplotypes.md)
