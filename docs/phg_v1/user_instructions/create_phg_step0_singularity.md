!!! warning "Legacy Documentation - PHG Version 1"

    This section contains documentation for **PHG version 1**, which is
    no longer actively developed. It is preserved here for archival and
    historical reference only. If you are looking to use the Practical
    Haplotype Graph, please refer to the [PHG v2 documentation](../../index.md),
    which reflects the current version of the software.

# Use Singularity to run PHG

## **IMPORTANT**

Because Singularity loads the user $HOME folder, it also picks up any environment variables the user has set.   When running any scripts that load liquibase, this can be a problem if the JAVA_HOME variable is set and it is different from what is in the container.  To avoid errors it is best to unset this environment variable.  To do this on a linux machine, please type the following before running your Singularity containers:


**unset JAVA_HOME**

Also, a note on **Singularity binding**.  Singularity allows the user to bind using 2 different syntaxes.  The first specifies just a source destination with the assumption there exists a folder with the same name in the singularity container.  The second option is to bind using the source:destination syntax.  This syntax explicitly specifies the container-specific destination folder where the source folder will be mounted.  When no colon is present Singularity assumes the source and destination are the same.  

While the no-colon method works for plugins, it doesn't play well with Liquibase, which expects folders to be relative to /Libuibase.  Because of this, we recommend users always bind with the source:destination method as seen in the example below:

** singularity exec -B ${WORKING_DIR}:/phg **

## Quick Start

1. Follow the [Quick Start](https://sylabs.io/guides/3.5/user-guide/quick_start.html) installation steps to install Singularity on your machine.
2. Run the following commands to create a singularity image based on the PHG docker images. The command below downloads the version 0.0.22 phg image from docker hub and wraps it in a singularity image named "phg_22.simg".  You should replace the phg version tag with the tag you wish to download.
```
singularity build phg_22.simg docker://maizegenetics/phg:0.0.22

```

## Details
You can use [Singularity](https://sylabs.io/docs/) containers to run docker containers without installing Docker. The benefit of using Singularity is that it does not require root access on your machine.

Singularity and Docker work well together, and you can run the PHG Docker container through singularity without needing Docker installed. 

Singularity communicates with [Docker Hub](https://hub.docker.com) to pull a Docker image into the Singularity Docker Registry.

As with Docker, Singularity containers give users the ability to work in a uniform and consistent computing environment which facilitates reproducible science.

Here are example Singularity PHG commands that work on a TACC machine.  Note the "export PATH=${PATH}:/sbin" command may not be needed.  At the time of this writing, there was a singularity issue on TACC machines that required this step.


Also note: "/scratch1/07005/lynnjo"  should be replaced by your own scratch directory path, which can be found by doing a "cds" followed by "pwd".

For long singularity runs (e.g. loading haplotypes to the db, or imputing values) a slurm job should be executed.  See TACC instructions for further details.

```
#!python

> cds
> mkdir -p phg_scratch/phg_run1
> export PATH=${PATH}:/sbin
> module load tacc-singularity
> cd phg_scratch
> singularity build phg_22.simg docker://maizegenetics/phg:0.0.22
> 
> singularity exec -B /scratch1/07005/lynnjo/phg_scratch/phg_run1/:/phg/ /scratch1/07005/lynnjo/phg_scratch/phg_22.simg /tassel-5-standalone/run_pipeline.pl -debug -Xmx1G -MakeDefaultDirectoryPlugin -workingDir /phg/ -endPlugin




```




## Troubleshooting

1. If you are having trouble downloading the Docker image with Singularity, check the Singularity [troubleshooting guide](https://sylabs.io/guides/2.5/user-guide/troubleshooting.html).


[Return to Step 0 pipeline](create_phg_step0_main.md)

[Return to PHG version 1 instructions](../home_variants_in_gvcf_files.md)

[Return to Wiki Home](../home.md)