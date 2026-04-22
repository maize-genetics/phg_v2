!!! warning "Legacy Documentation - PHG Version 1"

    This section contains documentation for **PHG version 1**, which is
    no longer actively developed. It is preserved here for archival and
    historical reference only. If you are looking to use the Practical
    Haplotype Graph, please refer to the [PHG v2 documentation](../../index.md),
    which reflects the current version of the software.

# SAMToMappingPlugin

## Overview

This plugin allows you to create ReadMappings from an existing SAM/BAM file and upload them to the DB.  This uses a KeyFile similarly to FastqToMappingPlugin, but will skip the Minimap2 Alignment step.

### Steps:
1. Loop through each SAM file in the key file
2. Apply a set of filters to the mappings
3. Push the readMapping to the DB
4. Generate a Path KeyFile for use with Path finding.

## Example Run Command

```
#!bash

time /tassel-5-standalone/run_pipeline.pl -debug -Xmx100G -configParameters [CONFIG_FILE] -HaplotypeGraphBuilderPlugin -configFile [CONFIG_FILE] -methods [HAPLOTYPE_METHOD] -includeVariantContexts false -includeSequences false  -endPlugin -SAMToMappingPlugin -endPlugin

```




## Config Parameters

```
#!bash

SAMToMappingPlugin.keyFile = Name of the Keyfile to process.  Must have columns cultivar, flowcell_lane, filename, and PlateID.

SAMToMappingPlugin.samDir = Name of the SAM/BAM dir to process.

SAMToMappingPlugin.maxRefRangeErr = Maximum allowed error when choosing best reference range to count.  Error is computed 1 - (mostHitRefCount/totalHits)

SAMToMappingPlugin.lowMemMode = true(default) or false.  If true, will run in a memory efficient manner with minimal performance hit.  If the reads are not ordered in the SAM file, you must use false and the whole file must fit in RAM

SAMToMappingPlugin.methodName = Name of the ReadMappings to be stored in the DB.

SAMToMappingPlugin.methodDescription = Description of the method

SAMToMappingPlugin.debugDir = (Optional)Directory on the file system that will write the ReadMapping files to disk.

SAMToMappingPlugin.outputSecondaryStats = True or false(default).  If true this will out put some additional mapping statistics which can be used to debug.  The files will be written to the current working directory.

```
All of these parameters can also be set on the command line as well.  A Path keyfile will be exported which can then be used to find paths.