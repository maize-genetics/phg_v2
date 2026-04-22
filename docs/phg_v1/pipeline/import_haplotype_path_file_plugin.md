!!! warning "Legacy Documentation - PHG Version 1"

    This section contains documentation for **PHG version 1**, which is
    no longer actively developed. It is preserved here for archival and
    historical reference only. If you are looking to use the Practical
    Haplotype Graph, please refer to the [PHG v2 documentation](../../index.md),
    which reflects the current version of the software.

The ImportHaplotypePathFilePlugin takes the haplotypePath files produced by [HapCountBestPathToTextPlugin](hap_count_best_path_to_text_plugin.md) and imports them to a data structure to be used by the [PathsToVCFPlugin](paths_to_vcf_plugin.md).  The latter produces VCF files from the paths.

This plugin expects a PHG Haplotype Graph to be sent as an incoming parameter.  The plugin requires a PHG Haplotype Graph, which can be adquired by chaining the running of HaplotypeGraphBuilderPlugin and this plugin.  The output from the HaplotypeGraphBuilderPlugin will be input to the ImportHaplotypePathFilePlugin.  In addition, the ImportHaplotypePathFilePlugin provides the input to the PathsToVCFPlugin.  All 3 plugins are likely chained when called via the command line or a shell script.  See [ExportPath.sh](export_path.md),

The ImporthaplotypePathFilePlugin  takes the following parameters:

* -inputFileDirectory <Input File Directory> A directory containing the paths text files to imported and used to create the vcf files. (Default=null) (REQUIRED)


