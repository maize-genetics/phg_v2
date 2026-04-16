!!! warning "Legacy Documentation - PHG Version 1"

    This section contains documentation for **PHG version 1**, which is
    no longer actively developed. It is preserved here for archival and
    historical reference only. If you are looking to use the Practical
    Haplotype Graph, please refer to the [PHG v2 documentation](../../index.md),
    which reflects the current version of the software.

# Step 3a: Create ReadMappings using the MultisampleBAM Pipeline

This step can replace the read mapping step of the impute pipeline for Single End GBS Fastq files.  The idea is that this pipeline is faster as it involves batching up the Fastq files and running minimap2 on each batch individually.  This allows minimap2 to only have to load in the PHG index file once per batch(instead of once per sample).  Then after alignment is run, we run MultisampleBAMToMappingPlugin to split out the SAM/BAM file into the original taxa and upload them to the DB for use in path finding.

# Quick Start

1. Write a config file and fill in the unassigned values.
2. Create a keyfile with information about reads to be imputed.
3. Run MergeFastqPlugin
4. Run the alignment scripts
5. Run MultisampleBAMToMappingPlugin

# Details

## Writing a config file
Before running the pipeline, you need to create a config file and put it in "baseDir" or some other directory accessible by the Docker. The config file name does not have to be config.txt. It also can be in any subdirectory of baseDir as long as the same file name location is used for -configParameters in the full pipeline command. 
Two sample config files are provided, one for imputing haplotypes from sequence in [fastq or SAM files (_link here_)](../files/sample_fastq_imputation_config.txt) and another for imputing haplotypes from variants in a [VCF file (_link here_)](../files/sample_vcf_imputation_config.txt). 
Values of **UNASSIGNED** in a config file must be filled in before using it. See the section [Setting parameters in the config file](#markdown-header-setting-parameters-in-the-config-file) for more information about what those values should be.

## Creating a Key File

## MergeFastqPlugin

This plugin will take a key file and a directory full of single end GBS-like fastq files and will concatenated the reads together into batched fastq files.  It will also output a Grouping File which holds records of which reads are in which batched fastq file and a series of alignment script templates which can be used to run minimap2.

To run please use the following command:


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

For more detailed documentation please look [here](merge_fastq_plugin_detailed_docs.md):

## Run Alignment Scripts

The next step is to run the alignment scripts created by the MergeFastqPlugin.  Each script will run a number of alignment commands to align each batched fastq file to the PHG minimap2 index.  It is recommended that you run these in parallel as it will greatly speed up processing time.

To Execute the script run the following command

```bash

time ./minimap2Script_1.sh
```

Once all the alignments have been processed, we can move onto the MultisampleBAMToMappingPlugin.

## Run MultisampleBAMToMappingPlugin

This plugin will take the BAM files along with the grouping file made in MergeFastqPlugin and will re-separate out each alignment into the original source fastq.  It will then be converted into a ReadMapping file and then uploaded to the DB for use in Path Finding. 

To run this step run the following command:


```bash

time tassel-5-standalone/run_pipeline.pl -debug -Xmx100G -configParameters config.txt \
			-HaplotypeGraphBuilderPlugin \
					-configFile config.txt \
					-methods HAPLOTYPE_METHOD \
					-includeVariantContexts false -endPlugin \
			-MultisampleBAMToMappingPlugin \
					-keyFile keyFile.txt \
					-fastqGroupingFile /outputGroupingFile.txt \
					-samDir outputBAMDir/ \
					-methodName GBS_PATH_METHOD \
					-methodDescription GBS_PATH_METHOD \
					-debugDir readMappingDebugDir/ \
					-outputSecondaryStats false \
					-isTestMethod true -endPlugin > multisampleBAM.log
```

To have this work correctly, you will need to set HAPLOTYPE_METHOD to match what you used to create the PHG Fasta file.  The -keyFile, -fastqGroupingFile and the -samDir need to match the parameters in MergeFastqPlugin. -methodName is the name of this PathMethod to be stored and the -methodDescription is the description to be stored in the DB.  Once this is done running you are able to run PathFinding to impute the haplotypes.