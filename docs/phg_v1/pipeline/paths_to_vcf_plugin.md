!!! warning "Legacy Documentation - PHG Version 1"

    This section contains documentation for **PHG version 1**, which is
    no longer actively developed. It is preserved here for archival and
    historical reference only. If you are looking to use the Practical
    Haplotype Graph, please refer to the [PHG v2 documentation](../../index.md),
    which reflects the current version of the software.

The PathsToVCFPlugin is used to create VCF files from haplotype paths.  The user may specific specific reference ranges to be included in the exported file.  This plugin would likely be chained with the HaplotypeGraphBuilderPlugin and the [ImportHaplotypePathFilePlugin](import_haplotype_path_file_plugin.md).  It takes input from both and uses this to create the requested VCF files.

The parameters to this plugin are:

* -outputFile <Output File> The file name for storing the VCF data. (Default=null) (REQUIRED)
* -RefRangeFileVCF <Ref Range File VCF> Reference Range file used to further subset the paths for only specified regions of the genome.  It not provided, all references ranges are used.  (Default=null) (OPTIONAL)
