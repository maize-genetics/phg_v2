!!! warning "Legacy Documentation - PHG Version 1"

    This section contains documentation for **PHG version 1**, which is
    no longer actively developed. It is preserved here for archival and
    historical reference only. If you are looking to use the Practical
    Haplotype Graph, please refer to the [PHG v2 documentation](../../index.md),
    which reflects the current version of the software.

# Separate genome intervals into different groups


## Details

Once you have created the initial set of genome intervals, you can group these intervals in various ways. This step lets you add multiple reference range groups after the original reference range loadings. It does not add different reference ranges; only combines the existing reference ranges in different groups. To add different reference range coordinates you need to create an entirely new database.

## Kitchen Sink

The AddRefRangeGroupPlugin method creates groupings for the reference ranges defined in your database. You can create as many sets of groupings in the database as you would like. 

There are 4 parameters used by this step:

* configFile: Path to file containing database access information. The config file contains separate lines for host=X, user=X, password=X, DB=X, and DBtype=X where X is defined by the user, and DBtype is either sqlite or postgres. (required)
* methodName: Method name for this reference range group.  (required)
* methodDetails: Description for this group of reference ranges. (required)
* ranges: This should be a tab-delimited bed formatted file with a subset of reference ranges currently existing in the DB. The file should contain chromosome, reference range start position, and reference range end position for the intervals you want to include in this group. There should not be a header line. (required)


### *Details on running this step through docker*

When AddRefRangeGroupPlugin is run as part of a Docker container script, the Docker script expects the following directory mount points.  It assumes the directory structure used matches that created from the MakeDefaultDirectoryPlugin, and that both the config.txt and ranges.txt files live at the top level in the folder mounted to /phg:

* mount point 1: folder that contains the configFile and ranges file.


An example Docker script to run the AddRefRangeGroupPlugin plugin is:
```
docker run --name group_intervals --rm \
	-v /localDir:/phg
	-t maizegenetics/phg:latest \
	/tassel-5-standalone/run_pipeline.pl -debug -AddRefRangeGroupPlugin -configFile /phg/config.txt -methodName ${methodName} -methodDetails ${methodDetails} -ranges ${rangesFile}
```
The --name parameter provides a name for the container.  This is optional.

The --rm parameter indicates the container should be deleted when the program finishes executing.  This is optional.

The -v directives are used to mount data from the user machine into the Docker.  The path preceding the ":" is the path on the user machine.  The directory path following the ":" are the paths inside the Docker where the user home directories will be mounted.

The -t directive indicates the Docker image of which this container will be an instance.  

The last line tells the Docker container to run the GroupGenomeIntervals.sh script which is found in the root directory.  The items following are the parameters to the GroupGenomeIntervals.sh script.


### *Files*

**Config file**

An example can be found here: Master config file

**Ranges file**

This is a tab-delimited bed file with a subset of the reference ranges used to create the database. The file should contain chromosome, reference range start position, and reference range end position for the intervals you want to include in this group. There should not be a header line. There is no limit to the number of ranges that can be included in this file, but reference ranges in this file must match ranges already in the database (added in the LoadGenomeIntervals step). A reference range can be added to multiple groups by running this step multiple times with different method names.


#### **AddRefRangeGroupPlugin**

The plugin takes a bed file with a list of reference ranges in the database, a method name for that group of ranges, a description of that group of ranges, and the configuration file indicating database name and type. The plugin reads the bed file into an object and all reference ranges from the database. It verifies that all entries from the bed file exist in the reference_ranges table (it will throw an exception if they are not all present). Next, the plugin verifies that the user-defined method name does not already exist in the methods table. If it does not, the plugin will add that method name to the database and add the specified reference ranges to the ref_range_ref_range_methods table with a new method ID.

Database tables updated via this method are:

* methods
* ref_range_ref_range_method

The parameters to this plugin are:

* -methodName <Method Name> Method name for this reference range group. (REQUIRED)
* -methodDetails <Method Details> Description for this group of reference ranges. (REQUIRED)
* -ranges <Ranges File> Tab-delimited, BED Formatted file containing chrom, ref range start position, ref range end position. No header line. (REQUIRED)
* -configFile <DataBase Configuration File> Path to file containing database access information, separate lines for host=X, user=X, password=X, DB=X, DBtype=X where X is user defined, and DBtype is either sqlite or postgres. (REQUIRED)

This plugin expects a database connection as input. This can be obtained by chaining the GetDBConnectionPlugin with the AddRefRangeGroupPlugin when called from the command line with the TASSEL run_pipeline.pl script.

An example of chaining these plugins is below:
```
#!python

/tassel-5-standalone/run_pipeline.pl -Xmx10G -debug \
-GetDBConnectionPlugin \
	-config ${dbConfigFile} \
	-create true  \
	-endPlugin \
-AddRefRangeGroupPlugin \
	-config ${dbConfigFile} \
	-ranges ${rangesFile} \
	-methodName ${method_name}  \
	-methodDetails ${method_details} \
	-endPlugin
```

## Troubleshooting

1. Check file paths and names for mistakes.
2. Double check that all reference ranges in the ranges file are present in the database (an easy way to check this is by verifying these ranges are in your original reference range intervals file). 
3. Make sure the methodName you are trying to use is unique and not already in the database methods table.


[Return to Step 1 pipeline version 0.0.40 or earlier](create_phg_step1_2_main.md)

[Return to Step 1 pipeline version 1.0 or later](create_phg_step1_create_db_load_ref.md)

[Return to Wiki Home](../home.md)
