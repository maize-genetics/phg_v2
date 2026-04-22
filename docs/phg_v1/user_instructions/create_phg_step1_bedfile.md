!!! warning "Legacy Documentation - PHG Version 1"

    This section contains documentation for **PHG version 1**, which is
    no longer actively developed. It is preserved here for archival and
    historical reference only. If you are looking to use the Practical
    Haplotype Graph, please refer to the [PHG v2 documentation](../../index.md),
    which reflects the current version of the software.

# Create bedfile of intervals for PHG reference ranges

## Intervals file documentation

The input to the loadGenomes pipeline is a tab-delimited file in bedfile format (positions are 0-based, start position is inclusive, end position is exclusive). The file must contain at least four columns. The first column should indicate chromosome, the second column indicates start position,  the third column indicates end position and the fourth column indicates the group name to which the interval belongs. The group name column should be a single word text.  The database will create reference range groups containing the intervals associated with each group from this fourth column.

Other columns may be present but will be ignored. Header lines must begin with # and will be skipped. Documentation for the bed file format can be found [here](https://m.ensembl.org/info/website/upload/bed.html).

If column headers are present, the first four must be named "chrom","chromStart", "chromEnd", and "name" as per bed file format.

You can break the reference genome into any set of intervals you want. In later steps these intervals can also be organized into different groups. For example, you can use a bed file created based on gff gene regions and organize these regions into genic and intergenic intervals in your database. If you have an existing bed file that includes only the regions on which you want to focus, but would like the non-represented regions of the genome to be included as well, you may run this bed file through the CreateValidIntervalsFilePlugin TASSEL plugin either within or outside of that docker.

In later steps, the PHG loadGenomeIntervals step will load reference genome and assembly sequence for both the intervals defined in the provided intervals file. WGS data is only added for the intervals defined in the bed file. If you would like WGS data added for the entire genome, you must create a bed file with no gaps between the intervals.

## Things to consider when creating a PHG bed file

- How large is the genome you're working with?
- What genome intervals are you most interested in studying?
- How diverse is the species you're working with? Do you expect large IBD regions?

## CreateValidIntervalsFilePlugin

The CreateValidIntervalsFilePlugin allows the user a convenient way to both check their intervals file for overlaps, and to create intervals for regions of the reference genome not covered by the initial intervals file.  The output from this script is a properly formatted intervals file that can be used as input when loading reference ranges to the PHG database. 


The CreateValidIntervalsFilePlugin has 6 parameters: 
 
* referenceFasta: path to reference fasta file (required)
* IntervalsFile: original bed file with intervals, id official bedfile format (required)
* generatedFile: name with path for intervals file created by the plugin (required)
* mergeOverlaps: true or false indicating whether to merge overlapping intervals (optional - defaults to false)
* userRegionsGroupName: name to give group of reference intervals provided in user file (optional - defaults to FocusRegion)
* otherRegionsGroupName: name to give for any intervals that are created representing regions not in the original intervals file (optional - defaults to FocusComplement)

An example of a script to run this plugin in a docker is below.  Note the paths assume the user has run the MakeDefaultDirectoryPlugin to create the directory structure:


```
#!java
WORKING_DIR=/workdir/user/phg_new
DOCKER_CONFIG_FILE=/phg/config.txt

docker1 run --name test_assemblies --rm  \
    -v ${WORKING_DIR}/:/phg/ \
    -t maizegenetics/phg:0.0.24 \
    /tassel-5-standalone/run_pipeline.pl -Xmx100G -debug -configParameters ${DOCKER_CONFIG_FILE} \
    -CreateValidIntervalsFilePlugin -intervalsFile /phg/inputDir/reference/original_intervals.bed \
    -referenceFasta /phg/inputDir/reference/ref.fa \
    -mergeOverlaps true \
    -generatedFile /phg/validBedFile.bed -endPlugin
```

If not running in docker, the command line plugin CreateValidIntervalsFilePlugin may be invoked from TASSEL in a similar manner.
If the user chooses not to have overlapping intervals merged and overlapping intervals are found, the program will end with an error message.  A list of the overlapping intervals will be printed.

Here is example command to call the CreateValidIntervalsFilePlugin from the command line outside a docker. You may, in addition, provide the optional parameters "-userRegionsGroupName" and "-otherRegionsGroupName".  The default values for these groups are "FocusRegion" and "FocusComplement".  By adding these parameters you may change the group names.

```
#!python

tassel-5-standalone/run_pipeline.pl -Xmx50G -debug -CreateValidIntervalsFilePlugin -intervalsFile myInitialBedFile.bed -referenceFasta ref.fa -mergeOverlaps true -generatedFile validBedFile.bed -mergeOverlaps true -endPlugin > plugin_output.txt
```

On a successful run, the output file will contain columns for chrom, chromStart, chromEnd and name.  Any overlapping intervals will have been merged, and intervals will have been created for areas of the reference genome not covered by the original intervals file.  All created intervals by this script will have the group name as provided in the group name parameter.

The user may then edit the provided file, removing any intervals they do not wish to be included.  It is strongly recommended that the whole genome be represented by these intervals.  When regions of the genome are missing, you risk losing assembly sequence that may be conserved, but did not align to ranges specified in your intervals file.

## Troubleshooting

1. Bed file formatting convention dictates that the bed file is 0-based, with an inclusive start position and exclusive end position. This differs from VCF file format, which is 1-based. Make sure to keep these coordinates in mind when creating your bed file.  Also, all header lines in the bedfile must begin with #.
2. The PHG cannot deal with overlapping reference range intervals. Gene models from gff files often overlap; we recommend merging all overlapping intervals into a single interval when creating the PHG reference range interval bed file. It is recommended the CreateValidIntervalsFile.sh script or CreateValidIntervalsFilePlugin (TASSEL command line) be used to create a valid file from your original intervals file.  Other option: [Bedtools](https://bedtools.readthedocs.io/en/latest/) offers one set of useful command line tools to create and manipulate bed files.


[Return to Step 1 pipeline version 0.0.40 or earlier](create_phg_step1_2_main.md)

[Return to Step 1 pipeline version 1.0 or later](create_phg_step1_create_db_load_ref.md)

[Return to Wiki Home](../home.md)