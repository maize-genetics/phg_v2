# Rare Alleles Pipeline

In this document, we will discuss a proposed pipeline to identify rare alleles in a population. 
This pipeline is invoked once a user has successfully run the PHG imputation pipeline and has
secured a haplotype VCF (hVCF) file. 

From the imputation h.vcf file, a composite fasta can be created.  The rare alleles pipeline begins
by re-aligning the WGS reads used in the imputation pipeline against the composite fasta derived above.
The BAM file created by this second alignment is filtered, then run through a variant caller.  We suggest using
DeepVariant (https://github.com/google/deepvariant) or Octopus (https://github.com/luntergroup/octopus)
for this purpose.  (BRANDON - SHOULD WE LINK THE PAGE TO THE VARIANT CALLER COMPARISONS HERE ??) The output of the variant caller is a VCF file that contains the variants called
against the composite fasta.

Additional filtering is then performed to select variants based on VCF quality scores, depth or other measures. This 
final VCF provides a list of SNPs that may be translated to haplotype coordinates.  To perfom this
transformation, the user can use the craete-haplotype-vcf command.


## Quick start

* From the imputation h.vcf file, create a composite fasta file
```shell
phg create-fasta-from-hvcf \
  --hvcf-file my/imputation/h.vcf \ 
  --fasta-type composite \ 
  -o /path/to/output_folder
```
* Align the reads you used in the imputation pipeline against the composite fasta using minimap2 or other aligner
```shell
# Run minimap2 to align the reads with the composite genome
FASTA_DIR=/workdir/lcj34/phg_v2/testing/fastas
FASTQ_DIR=/workdir/lcj34/phg_v2/testing/fastqFiles
WORK_DIR=/workdir/lcj34/phg_v2/testing

time minimap2 -a -x sr -t 20 --secondary=no   --eqx -R '@RG\tID:READ_GRP1\tSM:P39' \
    ${FASTA_DIR}/P39wgs_composite.fa \
    ${FASTQ_DIR}/P39/P39_HBEN2ADXX_GTGAAA_R1.fastq.gz \
    ${FASTQ_DIR}/P39/P39_HBEN2ADXX_GTGAAA_R2.fastq.gz | \
    samtools view -b > ${WORK_DIR}/p39top39composite_bams/P39_HBEN2ADXX_GTGAAA.ori.bam
time minimap2 -a -x sr -t 20 --secondary=no   --eqx -R '@RG\tID:READ_GRP2\tSM:P39' \
    ${FASTA_DIR}/P39wgs_composite.fa \
    ${FASTQ_DIR}/P39/P39__HFFGKADXX_GTGAAA_R1.fastq.gz \
    ${FASTQ_DIR}/P39/P39__HFFGKADXX_GTGAAA_R2.fastq.gz | \
    samtools view -b > ${WORK_DIR}/p39top39composite_bams/P39_HFFGKADXX_GTGAAA.ori.bam
time minimap2 -a -x sr -t 20 --secondary=no --eqx -R '@RG\tID:READ_GRP3\tSM:P39'
    ${FASTA_DIR}/P39wgs_composite.fa \
    ${FASTQ_DIR}/P39/P39_HL5WNCCXX_GTGAAA_R1.fastq.gz \
    ${FASTQ_DIR}/P39/P39_HL5WNCCXX_GTGAAA_R2.fastq.gz | \
    samtools view -b > ${WORK_DIR}/p39top39composite_bams/P39_HL5WNCCXX_GTGAAA.ori.bam



```

* Filter the bam files

Bam files can be filtered using samtools or other tools.  We suggest running "samtools fixmate", to fill in mate coordinates,
"samtools markdup -r" to remove duplicates and "samtools view -F 4 -f 2" to exclude unmapped reads or those that are not
properly paired.  The following is an example of how to filter the bam files created above.  These commands should be
run on each bam file created by the alignment step above.

```shell
samtools fixmate -m <sample>.ori.bam <sample>.fm.bam
samtools sort --threads 8 -T ./tmp -o <sample>.sortedFM.bam <sample>.fm.bam
samtools markdup -r <sample>.sortedFM.bam <sample>.dedup.bam
samtools view -F 4 -f 2 -b <sample>.dedup.bam > <sample>.filtered.bam
samtools index <sample>.filtered.bam
```


* Merge the bam files into 1

Before running the variant caller, merge the bam files into a single bam file, then index using the commands below:
```shell
samtools merge -b /path/to/bam.list /path/to/merged.bam
samtools sort --threads 8 -T ./tmp -o /path/to/merged.sorted.bam /path/to/merged.bam
samtools index /path/to/merged.sorted.bam
```
* Run DeepVariant, Octopus or some other caller on the aligned reads

Here we show an example of running DeepVariant on the merged bam file created above.  The output of DeepVariant is a VCF file
that contains the variants called against the composite fasta.  The following command will run DeepVariant 
from a Docker container.
  ```shell
  INPUT_DIR=/path/to/input/files
  OUTPUT_DIR=/path/to/output/files
  
  docker run \
        -v ${INPUT_DIR}:/input \
        -v ${OUTPUT_DIR}:/output \
        google/deepvariant:1.6.1 \
        /opt/deepvariant/bin/run_deepvariant \
        --model_type=WGS \
        --ref=/input/Composite.fa \ // from the create-fasta-from-hvcf command above
        --reads=/input/merged.sorted.bam \
        --output_vcf=/output/merged.sorted.vcf
        --output_gvcf=/output/merged.sorted.g.vcf
        --intermediate_results_dir=/output/intermediate_results
        --num_shards=4
  ```
* Filter the VCF file created by the variant caller.

To identify rare alleles, the VCF file output by the variant caller should be filtered to select variants based on 
quality scores, depth or other measures.  The following is an example of how to filter the VCF file created above.
It filters for biallelic SNPs with a depth greater than 4 and a quality score greater than 20.

```shell
bcftools view -m2 -M2 -v snps -i 'QUAL>20 && GT="1/1" && DP>4' /path/to/merged.sorted.vcf -o /path/to/deepVariantFiltered.vcf
```

* Create a vcf in haplotype coordinates

Convert the VCF file created above to haplotype coordinates, using the create-haplotype-vcf command as shown below:
  ```shell
  phg create-haplotype-vcf \
      --path-hvcf /path/to/imputation.hvcf \
      --variant-vcf /path/to/vcf/from/DeepVariant \ // vcf created in filtered step above
      --sample-name nameForNewSample \ // samplename to use in the haplotype coordinates vcf
      --output-file /path/to/output.vcf
  ```

*** WHAT TO DO WITH THIS HAPLOTYPE VCF FILE? ***