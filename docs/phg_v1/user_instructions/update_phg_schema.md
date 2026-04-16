!!! warning "Legacy Documentation - PHG Version 1"

    This section contains documentation for **PHG version 1**, which is
    no longer actively developed. It is preserved here for archival and
    historical reference only. If you are looking to use the Practical
    Haplotype Graph, please refer to the [PHG v2 documentation](../../index.md),
    which reflects the current version of the software.

# Update an existing PHG database to latest schema with liquibase

## Quick Start

Running one of the pipeline plugins automatically uses liquibase to make sure the software and database versions are compatible. Running individual plugins, however, does not result in the running of liquibase. For these cases the liquibase check must be run independently as follows:

Download the latest version of the PHG Docker. If an earlier version of the PHG Docker is used, the liquibase check will make sure that the PHG database is compatible with that version. However, it is not possible to roll back to an earlier version of the database.  



WHen pulling the phg docker, it is best to indicate the specific version to pull.  While "latest" is an available tag, it may be confusing if later you need to know which specifig PHG version was run.  In the example below, the tag is "1.2".  Your tag may be different.

```
docker pull maizegenetics/phg:1.2

```

Run the liquibase plugin. You may need to change the path to the config file or the output directory. At a minimum the config file needs to contain the database access parameters. Also if you have an outputDir= parameter setting in the config file, it does not need to be included on the command line.  Note the tag version in the command below matches the version from the "docker pull" command above.

```
WORKING_DIR=local/directory/where/files/will/go/
docker run --name create_directory --rm \
	-v ${WORKING_DIR}/:/phg/ \
	-t maizegenetics/phg:1.2 \
	/tassel-5-standalone/run_pipeline.pl -debug -configParameters /phg/config.txt \
        -CheckDBVersionPlugin -outputDir /phg/outputDir -endPlugin \
        -LiquibaseUpdatePlugin -outputDir /phg/outputDir -endPlugin
```

After running the plugin check for the `run_yes.txt` file that will be output. This means your database is recent enough to allow an update to be applied.  If instead you see a `run_no.txt` file, your database cannot be updated to the current PHG db schema.  In this case, you can try using an older version of the PHG, one that matches the version used when the db was last created/updated. If it was created with a variants table, the latest version that can be run against it is the PHG 0.0.40 version.


If the update was attempted, check the liquibase_output.log and liquibase_error.log files to verify the update has completed successfully.

## Details
Because the PHG is under active development, the database schema may change periodically. In most cases when this happens it will be possible to update the existing database with the [liquibase migration tool](https://www.liquibase.org), and the database will not need to be created from scratch.  Note that an update is not possible between 0.0.40 and earlier versions to PHG versions 1.0 and greater.

## Kitchen Sink

The PHG uses the liquibase data migration tool to update user databases. This tool works on both sqlite and postgresql databases.

When [creating a new database](create_phg_step1_2_main.md), the db is current and there is no need to run the liquibase Docker. However, the liquibase check should be run whenever a user is pulling a new PHG image to run with an existing database. Users may pull the most recent version or any preceding version, though a more recent database cannot be rolled back to an earlier version. 

If a "run_no.txt" file is created, it means your database is too old to be updated. In this scenario, you must start the database fresh. The run_no.txt file will contain a message indicating what schema is missing from the database, which should give an indication of the database versioning.

If the run_yes.txt file is present, then liquibase attempted to perform the database migrations.  Check the liquibase_error.log and liquibase_output.log files in your mounted output directory to verify the migrations ran without error.

The *Quick Start* command calls the LiquibaseUpdatePlugin. This plugin runs the liquibase "update" command. The software identifies database changes that have already been run against the database and skips them. 

## Troubleshooting

1. If you try to update your PHG database and get a `run_no.txt` file output, your database is likely too old to update. In this case you will need to create a new database because your current database is not guaranteed to work in subsequent imputation steps.

[Return to PHG version 1.0 or later home](../home_variants_in_gvcf_files.md)

[Return to PHG version 0.0.40 or earlier home](../home_variant_tables.md)

[Return to Wiki Home](../home.md)