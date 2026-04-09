!!! warning "Legacy Documentation - PHG Version 1"

    This section contains documentation for **PHG version 1**, which is
    no longer actively developed. It is preserved here for archival and
    historical reference only. If you are looking to use the Practical
    Haplotype Graph, please refer to the [PHG v2 documentation](../../index.md),
    which reflects the current version of the software.

# Populate PHG DB Pipeline Plugin

PopulatePHGDBPipelinePlugin is the second step to create a PHG DB.

If you provide assembly haplotype files(and a keyfile), the assemblies will be aligned to the reference genome.  There is a choice of using either mummer4 as an aligner (the original code) or using the newer anchorage aligner.  The anchorwave alignment has been shown to be more sensitive aligner than other aligners.  But it is a more time intensive step. Using mummer4 will be quicker but not as sensitive as anchorwave.  The pipeline uses mummer4 as the default.  

You can change the default assembly aligner by setting the "asmAligner" parameter to "anchor4" when invoking the PopulatePHGDBPipelinePlugin plugin, however it is not recommended to used anchorwave with this particular plugin.  The initial alignment step is slow and takes alot of memory.  We suggest, when using anchorwave, that the [AssemblyMAFFromAnchorWavePlugin](create_phg_step2_assembly_anchor_wave_plugin_details.md) be run on your assemblies.  The resulting MAF files should then be processed into GVCF files using the [MAFToGVCFPlugin](create_phg_step2_maf_to_gvcf_plugin_details.md), and then loaded to the database using the [LoadHaplotypesFromGVCFPlugin](create_phg_step2_load_haplotypes_from_gvcf_plugin_details.md).

> [!WARNING]
> The MAFToGVCFPlugin currently available with the PHG only supports MAF files created from anchorwave version 1.2.2
> or earlier.  Later versions of the anchorwave aligner include changes to reverse strand processing
> that have not yet been incorporated into the PHG software.
> Please bear this in mind when using the MAFToGVCFPlugin from PHG or you may have erroneous results.

If WGS fastq and a wgs keyfile are provided, the sequence will be aligned to reference with bwa mem and the resulting haplotypes stored in the database. Alternatively, aligned sequence can be loaded from BAM or GVCF files (and a keyfile). 

Both haplotype sequences and variants called from reference are stored in the DB.

Once some haplotype information has been loaded into the DB, similar haplotypes can be clustered to create a set of consensus haplotypes. This is computationally and biologically useful and is highly recommended.

## Quick Start Using the PHG Docker
The value of WORKING_DIR should be the same as that used to run MakeDefaultDirectory.
```
WORKING_DIR=/local/directory/where/MakeDefaultDirectory/was/run/
DOCKER_CONFIG_FILE=/phg/config.txt

docker run --name load_haplotypes --rm \
    -v ${WORKING_DIR}:/phg/ \
    -t maizegenetics/phg:latest \
    /tassel-5-standalone/run_pipeline.pl -Xmx100G -debug -configParameters ${DOCKER_CONFIG_FILE} \
    -PopulatePHGDBPipelinePlugin -endPlugin
```

## Details

This pipeline focuses mostly on creating haplotypes and loading them into the DB as well as merging similar haplotypes during the Create Consensus step.  Depending on the input data, different plugins will be run.
It makes use of the following plugins/Docker Scripts:

```
AssemblyHaplotypesMultiThreadPlugin
AssemblyMAFFromAnchorWavePlugin
CreateHaplotypesFromFastq.groovy
CreateHaplotypesFromBAM.groovy
CreateHaplotypesFromGVCF.groovy
HaplotypeGraphBuilderPlugin
LoadHaplotypesFromGVCFPlugin
MAFToGVCFPlugin
RunHapConsensusPipelinePlugin
```
Further information about the plugins and scripts along with detailed parameter descriptions can be found here:

[Align assemblies using anchorwave aligner](create_phg_step2_assembly_via_anchorwave.md)

[Align assemblies using mummer4 aligner](create_phg_step2_assembly_haplotype_parameters.md)

[WGS Haplotype Creation](create_phg_step2_add_haps_from_fastq.md)

[Making Consensus Haplotypes](create_phg_step2_consensus.md)

This pipeline will make use of parameters stored in a config file.  There are a lot of optional parameters for these plugins.  Please refer back to the links above for more information.

A sample config file with a basic set of pipeline parameters follows. The first 5 values are for database access.  They should be changed to match values specific to your database.  The database type must be either "sqlite" or "postgres".  Not all optional parameters are included:

```
host=localHost
user=sqlite
password=sqlite
DB=/phg/phg_db_name.db
DBtype=sqlite

numThreads=5

#Liquibase output directory
liquibaseOutdir=/phg/outputDir

# Haplotype creation parameters
referenceFasta=/phg/inputDir/reference/Ref.fa

# Assembly Loading Parameters
# Set the desired assembly aligner when running via PopulatePHGDBPopelinePlugin
PopulatePHGDBPipelinePlugin.asmAligner=anchorwave

#Assembly Loading Parameters when running the mummer4 aligner
asmMethodName=anchorwave
asmKeyFile=/phg/asm_keyFile.txt
outputDir=/phg/outputDir/align/
gvcfOutputDir=/phg/outputDir/align/gvcfs/

#Assembly Loading Parameters when running the anchorwave aligner
# Note running the anchorwave aligner from PopulatePHGDBPipelinePlugin involves
# the execution of 3 plugins.  PopulatePHGDBPipelinPlugin needs the AssemblyMAFFromAnchorWavePlugin 
# parameters to be defined in the config file.  It will infer the parameters needed for the other 3
# plugins from the AssemblyMAFFromAnchorWavePlugin parameters. It will also create the keyfile needed
# for both the MAFToGVCFPlugin and LoadHaplotypesFromGVCFPlugin classes.
AssemblyMAFFromAnchorWavePlugin.outputDir=/phg/outputDir
AssemblyMAFFromAnchorWavePlugin.keyFile=/phg/Ia453_K0326Y_keyfile.txt
AssemblyMAFFromAnchorWavePlugin.gffFile=/phg/inputDir/reference/Zm-B73-REFERENCE-NAM-5.0_Zm00001e.1.gff3
AssemblyMAFFromAnchorWavePlugin.refFasta=/phg/inputDir/reference/Zm-B73-REFERENCE-NAM-5.0.fa
AssemblyMAFFromAnchorWavePlugin.threadsPerRun=4
AssemblyMAFFromAnchorWavePlugin.numRuns=2


# WGS Loading Parameters
wgsMethodName=GATK_PIPELINE
haplotypeMethodName=GATK_PIPELINE
refRangeMethods=FocusRegion,FocusComplement
wgsKeyFile=/phg/keyFile.txt
gvcfDir=/phg/inputDir/loadDB/gvcf/
extendedWindowSize=0

# WGS Filtering Parameters
GQ_min=50
DP_poisson_min=.01
DP_poisson_max=.99
filterHets=true

# CONSENSUS Parameters
inputConsensusMethods=GATK_PIPELINE
consensusMethodName=CONSENSUS
minFreq=0.05
maxClusters=30
minSites=30
minCoverage=0.1
minTaxa=1
mxDiv=0.001
```

This example config file follows the directory structure defined by [MakeDefaultDirectoryPlugin](make_default_directory.md) but is the Docker specific version of the paths.  Docker paths all start with "/phg/". If running outside of the docker, please make the necessary adjustments.

### Controlling Which Steps Are Run

Running the assembly loading step requires an assembly key file. To skip loading haplotypes from assemblies, comment(#) 
or delete the line in the config file that starts with "asmKeyFile". 
Running the WGS loading step requires a WGS key file. To skip loading haplotypes from WGS, comment (#) 
or delete the line in the config file that starts with "wgsKeyFile".
Running the step that creates consensus haplotypes requires a consensus method name. To avoid running consensus, 
comment(#) or delete the line starting with "consensusMethodName".
The parameter "inputConsensusMethods" determines which haplotypes are used as the basis of the consensus haplotypes. 
In the sample config file, that is set to GATK_PIPELINE, which creates the consensus from the WGS haplotypes. 
Setting inputConsensusMethods=mummer4 would create consensus haplotypes from the assembly haplotypes.

If you have already created/loaded haplotypes to the database and wish to only run consensus, then comment or delete both the "asmKeyFile" and "wgsKeyFile" lines.
