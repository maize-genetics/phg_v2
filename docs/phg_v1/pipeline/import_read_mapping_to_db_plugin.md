!!! warning "Legacy Documentation - PHG Version 1"

    This section contains documentation for **PHG version 1**, which is
    no longer actively developed. It is preserved here for archival and
    historical reference only. If you are looking to use the Practical
    Haplotype Graph, please refer to the [PHG v2 documentation](../../index.md),
    which reflects the current version of the software.

# ImportReadMappingToDBPlugin

## Overview

This plugin allows the user to *merge* the ReadMappings stored in one DB into another.  One DB will be used as the final DB and all the other DBs needing to be merged will have the ReadMappings extracted and then written to the final DB.

### Steps:
1. Loop through each DB in -inputMappingDir(or readMapping file if -loadFromDB = false)
2. Write the ReadMapping to final DB.
3. (Optional) Generate a Path KeyFile for use with Path finding.

## Example Run Command

```
#!bash

tassel-5-standalone/run_pipeline.pl -debug -Xmx100G -ImportReadMappingToDBPlugin -configFileForFinalDB final_DB_config.txt -loadFromDB true -inputMappingDir /path/to/config/files/ -readMappingMethod READ_MAPPING_METHOD -outputKeyFile /path/to/output/pathKeyfile.txt -endPlugin
```

If you ran FastqToMappingPlugin with the -debug flag set, you can also upload those ReadMapping text files.


```
#!bash
tassel-5-standalone/run_pipeline.pl -debug -Xmx100G -ImportReadMappingToDBPlugin -configFileForFinalDB final_DB_config.txt -loadFromDB false -inputMappingDir /path/to/readMapping/files/ -readMappingMethod READ_MAPPING_METHOD -outputKeyFile /path/to/output/pathKeyfile.txt -endPlugin

```


## Plugin Parameters

```
#!bash

-configFileForFinalDB = the config file that holds the DB connection information for the 'final' DB
-loadFromDB = (true or false) if true, the plugin expects config files in -inputMappingDir, otherwise, it expects Read Mapping files.
-inputMappingDir = Directory holding all the config or ReadMapping files to be uploaded to the DB
-readMappingMethod = This must match the Read Mapping Method used in the 'final' DB.  Otherwise a new method will be written
-outputKeyFile = an automatically generated path key file representing the new ReadMapping IDs written to the DB.  This allows you to immediately run BestHaplotypePathPlugin after this is done.
```

Once this is done, you can run BestHaplotypePathPlugin to get Paths.