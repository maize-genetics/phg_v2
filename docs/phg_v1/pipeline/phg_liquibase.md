!!! warning "Legacy Documentation - PHG Version 1"

    This section contains documentation for **PHG version 1**, which is
    no longer actively developed. It is preserved here for archival and
    historical reference only. If you are looking to use the Practical
    Haplotype Graph, please refer to the [PHG v2 documentation](../../index.md),
    which reflects the current version of the software.

PHG uses the liquibase data migration tool to update user databases.  This tool works on both sqlite and postgresql databases.

Liquibase is run via the maizegenetics/phg_liquibase Docker.  When creating a new database from LoadGenomeIntervals.sh, the db is current and there is no need to run the liquibase Docker. The data migration docker should be run whenever a user is pulling a new PHG image to run with an existing database.  In this case, both the maizegenetics/phg and the maizegenetics/phg_liquibase docker images should be pulled to the user system.  Users may pull the "latest" version or any preceding version, as long as both the phg and phg_liquibase dockers are the same version.

A sample shell script for running phg_liquibase docker on your db is:
```
#!python

 docker run --name liquibase_update --rm \
  -v /workdir/lcj34/liquibase/data/:/tempFileDir/data/ \
  -v /workdir/lcj34/liquibase/output/:/tempFileDir/outputDir \
  -t maizegenetics/phg_liquibase \
  /liquibase/RunLiquibaseUpdates.sh configFile.txt
```

In the above, mount the directory containing your configuration file to :/tempFileDir/data/
In the above, mount the directory where you would like output written to :/tempFileDir/outputDir
In the above, "configFile.txt" should be replaced with the name of your configuration file that lives in the directory you mounted to /tempFileDir/data.

If you are running with an SQLite database, the database file must live in the output directory you have mounted to /tempFileDir/outputDir.

If you are running with a postgreSQL database, the liquibase docker will access the database using information from the configuration file. Please see the "Rules of Thumb for Connecting to a PHG PostgreSQL docker image" section of the "PHG Postgres Docker" section for details on specifying the postgreSQL host and port parameters.

The liquibase migration tool will write 3 files to the specified output directory:
  liquibase_error.log
  liquibase_output.log
  run_yes.txt or run_no.txt

If a "run_no.txt" file is created, it means your database is too old to be updated.  In this scenario, you must start the database fresh.  The run_no.txt file will contain a message indicating what schema is missing from the database, which should give an indication of the database versioning.

If the run_yes.txt file is present, then liquibase attempted to perform the database migrations.  Check the liquibase_error.log and liquibase_output.log files in your mounted output directory to verify the migrations ran without error.

The LiquibaseUpdatePlugin called from the liquibase docker will run the liquibase "update" command.  The software identifies database changes that have already been run against the database and skips them. 

While there are not always database changes each time there is a PHG version change, a new liquibase docker image will be created each time a PHG image is created.  This is to remove confusion as to which liquibase docker matches a particular PHG version.  Always run the same version maizegenetics/phg and maizegenetics/phg_liquibase dockers together.