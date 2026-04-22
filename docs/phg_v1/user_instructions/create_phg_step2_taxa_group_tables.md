!!! warning "Legacy Documentation - PHG Version 1"

    This section contains documentation for **PHG version 1**, which is
    no longer actively developed. It is preserved here for archival and
    historical reference only. If you are looking to use the Practical
    Haplotype Graph, please refer to the [PHG v2 documentation](../../index.md),
    which reflects the current version of the software.

# Separate taxa  into different groups


## Details

Once you have taxa loaded to the database, you can group these taxa in various ways. This step lets you create multiple taxa groups after the original taxa has been loaded.  For example, for a maize database you may want to group taxa into NAM or or 282.  Individual taxons may belong to multiple groups.

## Kitchen Sink

The AddTaxaToTaxaGroupPlugin method creates groupings for the taxa defined in the genotypes table of your database. You can create as many sets of groupings in the database as you would like. 

There are 4 parameters used by this step:

* configFile: This file comes from the cached configuration parameters.  It is a path to file containing database access information. The config file contains separate lines for host=X, user=X, password=X, DB=X, and DBtype=X where X is defined by the user, and DBtype is either sqlite or postgres. (required) 
* groupName: Group name for this group of taxa.  (required)
* taxa: Comma separated list of taxa.  Taxa names must match the names existing in the PHG database genotypes table. (required)


The last 2 parameters may be defined on the command line or in the configuration file.

### *Details on running this step through docker*

When AddTaxaToTaxaGroupPlugin is run as part of a Docker container script, the Docker script expects the following directory mount points.  It assumes the directory structure used matches that created from the MakeDefaultDirectoryPlugin, and that both the config.txt and ranges.txt files live at the top level in the folder mounted to /phg:

* mount point 1: folder that contains the configFile and 


An example Docker script to run the AddTaxaToTaxaGroupPlugin plugin is:
```
docker run --name group_intervals --rm \
	-v /localDir:/phg
	-t maizegenetics/phg:latest \
	/tassel-5-standalone/run_pipeline.pl -Xmx10G -debug \
          -GetDBConnectionPlugin \
            -configParameters ${dbConfigFile} \
            -create false \
            -endPlugin \
          -AddTaxaToTaxaGroupPlugin -groupName ${methodName} -taxa ${taxa}
```
The --name parameter provides a name for the container.  This is optional.

The --rm parameter indicates the container should be deleted when the program finishes executing.  This is optional.

The -v directives are used to mount data from the user machine into the Docker.  The path preceding the ":" is the path on the user machine.  The directory path following the ":" are the paths inside the Docker where the user home directories will be mounted.

The -t directive indicates the Docker image of which this container will be an instance.  

The last line tells the Docker container to run the GroupGenomeIntervals.sh script which is found in the root directory.  The items following are the parameters to the GroupGenomeIntervals.sh script.


#### **AddTaxaToTaxaGroupPlugin**

The plugin takes a taxa exsiting in the database, a group name for that group of taxa and the configuration file indicating database name and type. The plugin reads the list of taxa, verifying that all exist in the db.  If any do not, an exception is thrown and the missing taxa are printed.  Next, the plugin creates the group if it doesn't already exist, and adds the taxa.  If the group exists and any of the taxa are already in the group, they are ignored.

Database tables updated via this method are:

* taxa_groups
* taxa_groups_genoid

The parameters to this plugin are:

* -groupName <Group  Name> name for this taxa group. (REQUIRED)
* -taxa <Taxa> Comma separated list of taxa.  Taxa names must match the names existing in the PHG database genotypes table. (REQUIRED)
* -configFile <DataBase Configuration File> Path to file containing database access information, separate lines for host=X, user=X, password=X, DB=X, DBtype=X where X is user defined, and DBtype is either sqlite or postgres. NOTE: - this file should come from the command line parameter cache when obtaining a database connection.

This plugin expects a database connection as input. This can be obtained by chaining the GetDBConnectionPlugin with the AddTaxaToTaxaGroupPlugin when called from the command line with the TASSEL run_pipeline.pl script.

An example of chaining these plugins is below:
```
#!python

/tassel-5-standalone/run_pipeline.pl -Xmx10G -debug \
-GetDBConnectionPlugin \
	-configParameters ${dbConfigFile} \
	-create true  \
	-endPlugin \
-AddTaxaToTaxaGroupPlugin \
	-ranges ${rangesFile} \
	-methodName ${method_name}  \
	-methodDetails ${method_details} \
	-endPlugin
```

## Troubleshooting

1. Check file paths and names for mistakes.
2. Double check that all taxa in the taxa list  are present in the database (check the value of field "line_name" in your database's genotypes table). 


[Return to Step 2 pipeline](create_phg_step1_2_main.md)

[Return to Wiki Home](../home.md)
