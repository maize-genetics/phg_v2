!!! warning "Legacy Documentation - PHG Version 1"

    This section contains documentation for **PHG version 1**, which is
    no longer actively developed. It is preserved here for archival and
    historical reference only. If you are looking to use the Practical
    Haplotype Graph, please refer to the [PHG v2 documentation](../../index.md),
    which reflects the current version of the software.

# CreateHaplotype Key File

## Specification:

The Keyfile is a tab-separated text file which is used to set up the alignment, HaplotypeCaller and GVCF Upload steps for the PHG CreateHaplotypes Scripts.

The PHG will process the following columns:

| HeaderName  | Description  | Required |
|---|---|---|
| sample_name | Name of the taxon to be processed. | Yes |
| sample_description | Short Description of the sample_name.  | No, if not specified, an empty description will be used |
| files | Comma-separated list of file names to be processed. **This file list should not have the full paths pre appended.  But rather it needs just the file names.**  | Yes | 
| type | Type of the files to be processed.  PHG Currently Supports FASTQ, BAM or GVCF.  | Yes |
| chrPhased | Are the Chromosomes Phased?  This  needs to be 'true' or 'false' | Yes for GVCF type | 
| genePhased | Are the Genes Phased? This  needs to be 'true' or 'false' | Yes for GVCF type |
| phasingConf | What is the confidence of the phasing?  This needs to be between 0.0 and 1.0. If working with inbreds, this can be set close to 1.0. | Yes for GVCF type | 
| libraryID | What is the library ID of the fastq files.  This is only used if running BWA during CreateHaplotypesFromFastq.groovy | Yes only for FASTQ type | 
|gvcfServerPath | remote server name or address, followed by semicolon, followed by path on server where gvcfs files will be stored. | Yes if running with PHG version 1.0 or higher |


Because the entries in the 'files' column are comma separated, the PHG can do the following depending on the type:

* FASTQ : pairwise or single ended alignment using bwa mem.
* BAM : Run GATK/Sentieon HaplotypeCaller on all the BAM files specified in the list to create a single GVCF.
* GVCF : upload haplotypes for taxon with ploidy > 1.  Each file in the list will create a new haplotype.  If using Heterozygous material, we expect you to phase the GVCF file prior to running CreateHaplotypesFromGVCF.groovy.



## Sample File:
Note this file has the "gvcfServerPath" column required for PHG version 1.0 and above, but not required for  prior versions.

```
#!txt

sample_name     sample_description      files   type    chrPhased       genePhased      phasingConf     libraryID       gvcfServerPath
Ref     Ref line aligned        Ref_R1.fastq    FASTQ   true    true    .99     dummyLib1       localhost;/Users/lcj34/temp/gvcfRemote/
LineA   LineA line aligned      LineA_R1.fastq  FASTQ   true    true    .99     dummyLib1       localhost;/Users/lcj34/temp/gvcfRemote/
LineB   LineB line aligned      LineB_R1.fastq  FASTQ   true    true    .99     dummyLib1       localhost;/Users/lcj34/temp/gvcfRemote/
RefA1   RefA1 line aligned      RefA1_R1.fastq  FASTQ   true    true    .99     dummyLib1       localhost;/Users/lcj34/temp/gvcfRemote/
LineA1  LineA1 line aligned     LineA1_R1.fastq FASTQ   true    true    .99     dummyLib1       localhost;/Users/lcj34/temp/gvcfRemote/
LineB1  LineB1 line aligned     LineB1_R1.fastq FASTQ   true    true    .99     dummyLib1       localhost;/Users/lcj34/temp/gvcfRemote/

RefA1   RefA1 Aligned using BWA RefA1_dummyLib1_srt_dedup.bam   BAM     true    true    .99     null    localhost;/Users/lcj34/temp/gvcfRemote/
Ref     Ref Aligned using BWA   Ref_dummyLib1_srt_dedup.bam     BAM     true    true    .99     null    localhost;/Users/lcj34/temp/gvcfRemote/
LineB1  LineB1 Aligned using BWA        LineB1_dummyLib1_srt_dedup.bam  BAM     true    true    .99     null    localhost;/Users/lcj34/temp/gvcfRemote/
LineA   LineA Aligned using BWA LineA_dummyLib1_srt_dedup.bam   BAM     true    true    .99     null    localhost;/Users/lcj34/temp/gvcfRemote/
LineB   LineB Aligned using BWA LineB_dummyLib1_srt_dedup.bam   BAM     true    true    .99     null    localhost;/Users/lcj34/temp/gvcfRemote/
LineA1  LineA1 Aligned using BWA        LineA1_dummyLib1_srt_dedup.bam  BAM     true    true    .99     null    localhost;/Users/lcj34/temp/gvcfRemote/
RefA1   RefA1 Aligned using BWA RefA1_haplotype_caller_output_filtered.g.vcf.gz GVCF    true    true    .99     null    localhost;/Users/lcj34/temp/gvcfRemote/
Ref     Ref Aligned using BWA   Ref_haplotype_caller_output_filtered.g.vcf.gz   GVCF    true    true    .99     null    localhost;/Users/lcj34/temp/gvcfRemote/
LineB1  LineB1 Aligned using BWA        LineB1_haplotype_caller_output_filtered.g.vcf.gz        GVCF    true    true    .99     null    localhost;/Users/lcj34/temp/gvcfRemote/
LineA   LineA Aligned using BWA LineA_haplotype_caller_output_filtered.g.vcf.gz GVCF    true    true    .99     null    localhost;/Users/lcj34/temp/gvcfRemote/
LineB   LineB Aligned using BWA LineB_haplotype_caller_output_filtered.g.vcf.gz GVCF    true    true    .99     null    localhost;/Users/lcj34/temp/gvcfRemote/
LineA1  LineA1 Aligned using BWA        LineA1_haplotype_caller_output_filtered.g.vcf.gz        GVCF    true    true    .99     null    localhost;/Users/lcj34/temp/gvcfRemote/

```


[Return to Step 2 pipeline, version 1.0+](create_phg_step2_assembly_and_wgs_haplotypes.md)

[Return to Step 2 pipeline, version 0.0.40-](create_phg_step1_2_main.md)

[Return to Wiki Home](../home.md)