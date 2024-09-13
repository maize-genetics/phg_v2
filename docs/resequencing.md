# Rare Alleles Pipeline

In this document, we will discuss a proposed pipeline to identify rare alleles in a population. 
This pipeline is invoked once a user has successfully run the PHG imputation pipeline and has
secured a haplotype VCF (hVCF) file. 

From the imputation h.vcf file, a composite fasta can be created.  The rare alleles pipeline begins
by re-aligning the WGS reads used in the imputation pipeline against the composite fasta derived above.
The BAM file created by this second alignment is then run through a variant caller.  We suggest using
DeepVariant (https://github.com/google/deepvariant) or Octopus (https://github.com/luntergroup/octopus)
for this purpose.  The output of the variant caller is a VCF file that contains the variants called
against the composite fasta.

Filtering is then performed to select variants based on VCF quality scores, depth or other measures. This 
final VCF provides a list of SNPs that may be translated to haplotype coordinates.  TO perfom this
transformation, the user can use the craete-haplotype-vcf command, detailed below.

** NEED examples here of running DeepVariant and Octopus from Docker ***

## Quick start

* Create a vcf in haplotype coordinates
  ```shell
  phg create-haplotype-vcf \
      --path-hvcf /path/to/imputation.hvcf \
      --variant-vcf /path/to/vcf/from/DeepVariant \
      --sample-name nameForNewSample \
      --output-file /path/to/output.vcf
  ```