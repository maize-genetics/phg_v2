!!! warning "Legacy Documentation - PHG Version 1"

    This section contains documentation for **PHG version 1**, which is
    no longer actively developed. It is preserved here for archival and
    historical reference only. If you are looking to use the Practical
    Haplotype Graph, please refer to the [PHG v2 documentation](../../index.md),
    which reflects the current version of the software.

# Creating a Simple Example Database

Creating a simple database with simulated data provides a convenient way to test whether you are able to run the 
PHG software successfully. A simple reference genome along with fastq files of a few taxa can be created using the 
CreateSmallGenomesPlugin. The database can be created either using the PHG Docker or using the standalone version of TASSEL.
The following TASSEL command can be used to create the data after replacing "/path/to/basedir/" with a directory name. 

```
./run_pipeline.pl -debug -Xmx1G -CreateSmallGenomesPlugin \
   -baseDir /path/to/basedir -geneLength 3000 -interGeneLength 3500 \
   -numberOfGenes 10 -endPlugin
```

The command can be run without the -baseDir parameter, in which case the baseDir default will be *userDir*/temp/phgSmallSeq. 
If you use docker to run the command (below), userDir = (root). In that case using the default value, the data will be 
written to root not to /phg/ and the data will disappear when the container is removed. In short, when running 
CreateSmallGenomesPlugin using the docker **always** use the -baseDir parameter and start it with /phg/.

When using a docker to run the command, a directory called "basedir" will be created in the directory specified by 
"/Path/to" if it does not already exist. To run the command using the docker

```
docker run --name smallGenomes_container --rm -v /Path/to/:/phg/ -t maizegenetics/phg \
   /tassel-5-standalone/run_pipeline.pl -debug -Xmx1G \
   -CreateSmallGenomesPlugin -baseDir /phg/basedir \
   -geneLength 3000 -interGeneLength 3500 \
   -numberOfGenes 10 -endPlugin > /path/to/smallGenomes_log.txt
```

The plugin runs [MakeDefaultDirectoryPlugin](make_default_directory.md) and populates some of the resulting directories 
with simulated data representing a single chromosome with 10 (-numberOfGenes) genes. Specifically, it creates two sub-directories, "answer" and "dockerBaseDir"
in the folder specified by the -baseDir parameter. The "answer" directory holds data that can be used to validate imputation
results from running the PHG. The "dockerBaseDir" directory holds data that can be used to populate the PHG, map reads, 
and impute data. It will also hold files produced by running PHG. For the docker commands that follow, "dockerBaseDir" is the 
directory that should be mounted as "/phg/".

# Using the Example Database


The four pipeline plugins used for the PHG are

* [MakeDefaultDirectoryPlugin](make_default_directory.md)
* [MakeInitialPHGDBPipelinePlugin](make_initial_phgdb_pipeline.md)
* [PopulatePHGDBPipelinePlugin](populate_phgdb_pipeline.md)
* [ImputePipelinePlugin](impute_with_phg_main.md)

If you ran CreateSmallGenomesPlugin one of the things it did was to run MakeDefaultDirectoryPlugin so that is already 
done. Following are the docker commands and links to config files needed to run the other plugins. The first step of 
each pipeline is to make a call to liquibase to  make sure that the version of the PHG software and the PHG database are
 the same. As a result you need to have pulled the same versions of both maizegenetics/phg and maizegenetics/phg_liquibase. 
 
## MakeInitialPHGDBPipelinePlugin
 
The config file needed here was written to dockerBaseDir when CreateSmallGenomesPlugin was run. The docker command to run this pipeline is

```
docker run --name makeDb_container --rm -v /path1/to/dockerBaseDir:/phg/ -t maizegenetics/phg \
   /tassel-5-standalone/run_pipeline.pl -debug -Xmx1G -configParameters /phg/configSQLiteDocker.txt \
   -MakeInitialPHGDBPipelinePlugin -endPlugin > /path2/to/makedb_log.txt
```

Replace "/path/to" with the path that is correct for your computer. The path to the log file is on your computer and 
outside of docker. "maizegenetics/phg" is the name of the docker image. If the plugin ran correctly, "phgSmallSeq.db" will 
have been created in dockerBaseDir. Also, you should look at the log file to make sure there are no error messages. 

## PopulatePHGDBPipelinePlugin
 
The docker command to run this pipeline, after replacing both instances of "/path/to", is

```
docker run --name populateDb_container --rm -v /path/to/dockerBaseDir:/phg/ \
   -t maizegenetics/phg /tassel-5-standalone/run_pipeline.pl -debug -Xmx1G \
   -configParameters /phg/configSQLiteDocker.txt -PopulatePHGDBPipelinePlugin -endPlugin \
   > /path/to/populatedb_log.txt

```


This step may take several minutes. It is aligning fastq files to reference and extracting WGS based haplotypes using the 
script CreateHaplotypesFromFastq.groovy. It will populate the database with haplotypes and write some files to folders in 
dockerBaseDir in the process. As usual, check the log file to make sure there were no errors and that the plugin ran as expected.
 
## ImputePipelinePlugin
 
 The config file needed for this step is [imputeConfig.txt](../files/imputeConfig.txt). You will need to put a copy of that in dockerBaseDir. The docker command 
 to run this pipeline, after replacing both instances of "/path/to" and copying the config file to baseDockerDir, is

```
docker run --name impute_container --rm -v /path/to/dockerBaseDir:/phg/ -t maizegenetics/phg \
   /tassel-5-standalone/run_pipeline.pl -debug -Xmx1G -configParameters /phg/imputeConfig.txt \
   -ImputePipelinePlugin -imputeTarget pathToVCF -endPlugin > /path/to/impute_log.txt
```

The final result is a vcf file, dockerBaseDir/outputDir/output.vcf.

## Examine the data

An excellent way to look at the database contents is to use [rPHG](https://bitbucket.org/bucklerlab/rphg/wiki/Home).