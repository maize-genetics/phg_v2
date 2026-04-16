!!! warning "Legacy Documentation - PHG Version 1"

    This section contains documentation for **PHG version 1**, which is
    no longer actively developed. It is preserved here for archival and
    historical reference only. If you are looking to use the Practical
    Haplotype Graph, please refer to the [PHG v2 documentation](../../index.md),
    which reflects the current version of the software.

# MAFToGVCFPlugin Details

After aligning your assemblies and obtaining a MAF file from anchorwave, the next step is to convert
the MAF file to a GVCF file.  It is this latter format which will be used for loading
the alignment data as haplotypes to the PHG database.

> [!WARNING]
> The MAFToGVCFPlugin currently available with the PHG only supports MAF files created from anchorwave version 1.2.2
> or earlier.  Later versions of the anchorwave aligner include changes to reverse strand processing
> that have not yet been incorporated into the PHG software.
> Please bear this in mind when using the MAFToGVCFPlugin from PHG or you may have erroneous results.

This step  can be run using a shell script similar to the one detailed below:


```
WORKING_DIR=/workdir/lcj34/anchorwaveTesting/
DOCKER_CONFIG_FILE=/phg/config.txt

# The the config file must contain the necessary MAFToGVCFPLugin parameters
docker run --name anchorwave_assemblies_maf --rm \
    -v ${WORKING_DIR}/:/phg/ \
    -t maizegenetics/phg:latest \
    /tassel-5-standalone/run_pipeline.pl -Xmx100G -debug -configParameters ${DOCKER_CONFIG_FILE} \
    -MAFToGVCFPlugin -endPlugin
```


The plugin may also be run directly from tassel-5-standalone as per the command below.  In this example, parameters are explicitly defined, but they may alternately be pulled from a config file as with the example above.

```
/workdir/lcj34/tassel-5-standalone/run_pipeline.pl -Xmx100G -debug -MAFToGVCFPlugin -referenceFasta /workdir/lcj34/phg_nam_assemblies/B73/Zm-B73-REFERENCE-NAM-5.0.fa -mafFile /workdir/lcj34/MAF_files/-CML108_B73.maf -sampleName CML108_anchorwave -gvcfOutput /workdir/lcj34/GVCF_output/CML108ToB73.gvcf -fillGaps false > CML108_outputMafToGVCF.txt
```


If running with a config file, the parameters listed below should be defined (replace values with appropriate values for your own run).  Note the MAFToGVCFPlugin does not require database access.  This plugin neither reads not writes to a database.

```
MAFToGVCFPlugin.referenceFasta=/phg/inputDir/reference/Zm-B73-REFERENCE-NAM-5.0.fa
MAFToGVCFPlugin.mafFile=/phg/outputDir/align/CML108_B73.maf
MAFToGVCFPlugin.gvcfOutput=/phg/inputDir/loadDB/gvcfs/CML108ToB73.gvcf
MAFToGVCFPlugin.sampleName=CML108_anchorwave
MAFToGVCFPlugin.fillGaps=false
MAFToGVCFPlugin.twoGvcfs=false

```


# Details

This plugin takes as input a UCSC Multiple Alignment Format (MAF) file and creates a GVCF file from the MAF data.  It keeps track of SNP, indel and reference positions.  When available, it also keeps track of the assembly chromosome, the assembly start and end positions, and the assembly strand from the alignment.  For MAF files created from diploid alignments, it will create 2 gvcfs.  This is indicated when the user sets the "twoGvcfs" flag to "true".  When creating 2 gvcf files, _1 and _2 will be appended to the filename.

While this plugin should work to create a GVCF file from any correctly formatted MAF file, its intended use is for creating GVCF files from anchorwave created MAF files.  These GVCF files would then be loaded into a PHG database.  Because of this, the code requires that only a reference and a single sample are included in the MAF file.  It will ignore any MAF block lines except the "a" denoting the start of a MAF block, and the first 2 "s" lines, denoting the reference and single sample sequence lines.  Any MAF block "i", "e" and "q" lines are ignored.

## Parameters

Run the following command to print a current list of parameters for this plugin:

```
docker run --rm  maizegenetics/phg /tassel-5-standalone/run_pipeline.pl -MAFToGVCFPlugin
```

## Parameter Descriptions

* gvcfOutput: (required) Output GVCF file. If generating two output files, _1 and _2 will be appended to the filename.
* mafFile: (required) Input MAF file.  Please note that this needs to be a MAF file with 2 samples. The first will be assumed to be the Reference and the second will be the assembly.
* referenceFasta: (required) Full path to the reference fasta file
* sampleName: (required ) Sample name to write to the GVCF file.  This is also the name that will later be written to the database.
* fillGaps: (optional: default=false) When true, if the maf file does not fully cover the reference genome any gaps in coverage will be filled in with reference blocks. This is necessary if the resulting GVCFs are to be combined.
* twoGvcfs: (optional: default=false) The input maf was created from a diploid alignment and should be used to create two separate gVCF files.
* outputJustGT (optional: default=false) Output just the GT flag.  If set to false(default) will output DP, AD and PL fields.
* outputType (optional: default=gvcf) Output GVCF typed files. If set to gvcf(default) it will send all REFBlocks, SNPs and Indels.  If set to vcf, it will only output SNPs, Indels and missing values(if fillGaps = true.
