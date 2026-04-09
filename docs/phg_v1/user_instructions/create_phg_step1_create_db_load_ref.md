!!! warning "Legacy Documentation - PHG Version 1"

    This section contains documentation for **PHG version 1**, which is
    no longer actively developed. It is preserved here for archival and
    historical reference only. If you are looking to use the Practical
    Haplotype Graph, please refer to the [PHG v2 documentation](../../index.md),
    which reflects the current version of the software.

# Create and Populate a PHG database with haplotypes

Before creating the db, you need to setup the PHG directory structure.  We recommend using the default directory structure, which can be created as described in step 1A below. 
When the default directory is created, two of the files created are config.txt and load_genome_data.txt. Those files need to have values added or modified by the user before other steps can be run.

A sample config file can be viewed [here](https://bitbucket.org/bucklerlab/practicalhaplotypegraph/wiki/UserInstructions/SampleHaplotypeConfigFile)

There are two steps to building and populating a PHG database. The first step creates the PHG database and 
populates it with the reference genome. The second step populates the database with additional haplotypes that have been aligned to the
reference genome reference intervals.

# Step 1: Create database and load the reference genome and genome intervals

![CreateAndPopulatePHGDatabase_step1.png](../images/CreatePHG_step1.png)

This part of the pipeline populates the database with the reference genome, broken into "reference ranges" based on user-defined genome intervals. These same reference ranges are used to create haplotype blocks in later steps. The user can define different groups of reference ranges for different analyses.

# Run Step 1
# Use the links below to work through each step in the Step 1 flow chart

A. [Create the default directory structure](make_default_directory.md)

B. [Create a PHG database](make_initial_phgdb_pipeline.md)

C. [Create bed file to define genome intervals](create_phg_step1_bedfile.md)

D. [Optionally, set up additional groups of PHG intervals](create_phg_step1_group_genome_intervals.md)

Step 1D can be done at any time after the initial genome intervals have been loaded.


[Return to PHG version 1 Home ](../home_variants_in_gvcf_files.md)

[Return to Wiki Home](../home.md)
