!!! warning "Legacy Documentation - PHG Version 1"

    This section contains documentation for **PHG version 1**, which is
    no longer actively developed. It is preserved here for archival and
    historical reference only. If you are looking to use the Practical
    Haplotype Graph, please refer to the [PHG v2 documentation](../../index.md),
    which reflects the current version of the software.

# MergeFastqPlugin Detailed Documentation

This plugin will take a key file and a directory full of single end GBS-like fastq files and will concatenated the reads together into batched fastq files. It will also output a Grouping File which holds records of which reads are in which batched fastq file and a series of alignment script templates which can be used to run minimap2.

## Example Command

```bash
time ./tassel-5-standalone/run_pipeline.pl -Xmx100G -debug -MergeFastqPlugin \
						-fastqDir fastqDir/ \
						-outputDir outputDir/ \
						-outputBAMDir outputBAMDir/ \
						-makeBAMDir true \
						-outputGroupingFile outputGroupingFile.txt \
						-numToMerge 50 \
						-scriptTemplate minimap2Script \
						-numberOfScriptOutputs 2 \
						-numThreads 20 \
						-minimapLocation minimap2 \
						-minimap2IndexFile phgIndexFile.mmi \
						-keyFile keyFile.txt  -endPlugin > mergeFastq.log

```
In the above command, MergeFastqPlugin will make batches of 50 fastq files per batch and will be saved in -outputDir.  It will also make 2(based on -numberOfScriptOutputs) scripts labeled minimap2Script_1.sh and minimap2Script_2.sh which can be run to execute the Minimap2 commands.  In those scripts, it will use the -outputBAMDir to setup the correct BAM directory in the script, the -numThreads to set the -t parameter in minimap2, the -minimap2Location to specify where the minimap2 executable actually is and the -minimap2IndexFile to specify where the PHG index file is stored.  If all these parameters are correctly set, the scripts should be able to be run without any modification. 

## Parameter documentation

This plugin has the following Parameters available:

* -fastqDir(required): Name of the Fastq Directory to Process.  Must be an existing directory on the machine and must be filled with Single End GBS-like fastq files.
* -outputDir(required): Directory to write out the Merge Fastq files.
* -outputBAMDir(default: <bamFolder>):  This is the expected BAM directory written in the alignment scripts.  If left as the default <bamFolder>, an easy find and replace can update it to the real BAM directory.
* -makeBAMDir(default: true):  Option to add in a mkdir command to the alignment scripts.  If the directory already exists, this can be set to false.
* -outputGroupingFile(required): Output file to keep track of how the Fastq files were merged together.
* -numToMerge(default: 50): The number of fastq files to merge per batch.  Larger batches will require more processing time per batch, but overall less time as the minimap2 index only needs to be loaded once.
* -scriptTemplate(default: runMinimapTemp.sh):  This is the first portion of the output script name.  When this is run, a _1.sh, _2.sh ... will be appended to the end of the provided name.
* -numberOfScriptOutputs(default: 1): This sets the number of output scripts to write.  The Plugin will append _1.sh,_2.s ... to the end of the -scriptTemplate parameter.
* -numThreads(default: 20): This sets the number of threads parameter for minimap2 to use in the output scripts. 
* -minimap2Location(default: minimap2): This sets the minimap2 executable location in the output script.  Consider changing this if Minimap2 is not stored on the PATH.
* -minimap2IndexFile(default: <refIndex>): This sets the Index filename written to the output scripts.  By default it will write out <refIndex> for an easy find and replace.
* -keyFile(required) : File name for the genotyping keyfile.  This tab delimited keyfile must require the following columns: cultivar, flowcell_lane, and filename.  It must also have the headers as well. 
* -minimapN(default: 60): Integer which sets the minimap2 -N parameter in the output alignment scripts.  Please refer back to the minimap2 documentation to change this appropriately.
* -minimapf(default: "5000,6000"): String which sets the minimap2 -f parameter in the output alignment scripts.  Please refer back to the minimap2 documentation to change this appropriately. 
* -outputSams(default: false): If set to true, the output alignment scripts will output SAM files instead of BAM files.