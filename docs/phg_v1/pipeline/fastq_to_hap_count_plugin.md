!!! warning "Legacy Documentation - PHG Version 1"

    This section contains documentation for **PHG version 1**, which is
    no longer actively developed. It is preserved here for archival and
    historical reference only. If you are looking to use the Practical
    Haplotype Graph, please refer to the [PHG v2 documentation](../../index.md),
    which reflects the current version of the software.

The FastqToHapCountPlugin executes code to count the number of reads which align to a haplotype node. It scores which haplotypes are identical, excluded, or unresolved relative to a perfect hit GenotypeMap. The plugin requires a PHG Haplotype Graph, which can be acquired by chaining the running of HaplotypeGraphBuilderPlugin and this plugin.  The output from the HaplotypeGraphBuilderPlugin will be input to the FastqToHapCountPlugin.  For example:

```
#!java
     /tassel-5-standalone/run_pipeline_jbwa.pl -Xmx5g -debug -HaplotypeGraphBuilderPlugin 
          -configFile $CONFIGFILE -onlyAnchors true -consensusMethod $CONSENSUS_METHOD 
          -includeVariantContexts true -endPlugin 
          -FastqToHapCountPlugin -taxon $file -configFile $CONFIGFILE -haplotypeFile $HAPLOTYPE_FASTA 
          -method $HAP_COUNT_METHOD -version $VERSION -refFile $REFERENCE_FILE 
          -rawReadFile ${FASTQ_DIR}/${file}.fastq -exportHaploFile ${FASTQ_HAP_COUNT_DIR}/${file}.txt
```

The parameters for this plugin are:

* -configFile <Config File> DB Config File containing properties host,user,password,DB and DBtype where DBtype is either sqlite or postgres (Default=null) (REQUIRED) A sample config file can be found here:[Config File Wiki](https://bitbucket.org/bucklerlab/practicalhaplotypegraph/wiki/DockerPipeline/ConfigFile).
* -haplotypeFile <Haplotype File> Fasta file and associated BWA indices for haplotypes stored in the PHG db. (Default=null) (REQUIRED)
* -refFile <Ref File> Reference genome file - temporary need until we can back convert coordinates. (Default=null) (REQUIRED)
* -rawReadFile <Raw Read File> Raw Read file aligned to the reference. (Default=null) (REQUIRED)
* -exportHaploFile <Export Haplo File> Text file to store haplotype scoring. (Default=null) (OPTIONAL)
* -taxon <Taxon> Taxon name to load into genotypes table from linked plugin. (Default=null) (REQUIRED)
* -debugTaxon <Debug Taxon> Used for debugging. (Defaulg=null) (OPTIONAL)
* -method <Method> Name of method used to create hap counts, for the haplotype_counts table. (Default=null) (REQUIRED)
* -version <Version> Genome intervals version name as stored in the PHG db genome_intervals_version table. (Default=null) (REQUIRED)
* -loadDB <Load DB> Whether to populate the haplotype_counts table.  Should be true unless user is testing. (Default=true) (OPTIONAL)

Once the haplotype counts have been determined, they are loaded to the database if the loadDB flag is "true"..