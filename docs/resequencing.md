# Resequencing / Rare Alleles Pipeline

In this document, we will discuss the proposed pipeline to identify 
rare alleles in a PHG population. This pipeline is invoked once a user 
has successfully run the [PHG imputation pipeline](imputation.md) and 
has secured a haplotype VCF ([hVCF](hvcf_specifications.md)) file. 

From the imputation hVCF file, a composite FASTA can be created. The 
pipeline begins by re-aligning the WGS reads used in the imputation 
pipeline against the composite FASTA derived above. The BAM file 
created by this second alignment is filtered, then run through a 
variant caller. We suggest using
[DeepVariant](https://github.com/google/deepvariant) or 
[Octopus](https://github.com/luntergroup/octopus) for this purpose. 
The output of the variant caller is a VCF file that contains the 
variants called against the composite fasta.

Additional filtering is then performed to select variants based on 
VCF quality scores, depth or other measures. This final VCF provides 
a list of SNPs that may be translated to haplotype coordinates. To 
perform this transformation, the user can use the 
`composite-to-haplotype-coords` command.

## Quick start

## Detailed walkthrough

### Create a composite FASTA file

Using the [prior example outputs from the imputation section](imputation.md#find-paths),
we can first convert the hVCF output from the `find-paths` command
to a FASTA format via the `create-fasta-from-hvcf` command:

```shell
phg create-fasta-from-hvcf \
    --hvcf-file output/vcf_files_imputed/LineA_B.h.vcf \ 
    --fasta-type composite \ 
    -o output/composite_assemblies/
```

This command takes 3 necessary parameters:

* `--hvcf-file` - an hVCF file from the imputation pipeline
* `--fasta-type` - how should the output for the FASTA file look like?
    + `composite` - concatenates haplotype sequences together at the
      chromosome level (**needed for this pipeline**)
    + `haplotype` - each FASTA entry is a haplotype sequence
* `-o` or `--output-dir` - an _existing_ output directory for the 
  FASTA file

After running this command, our example directory (based on the
prior pipelines) will now have a FASTA file (`LineA_B_composite.h.vcf`) 
in the `composite_assemblies/` directory:

```
phg_v2_example/
├── data
│   └── short_reads
│   │   ├── LineA_LineB_1.fq
│   │   └── LineA_LineB_2.fq
│   ├── anchors.gff
│   ├── Ref-v5.fa
│   ├── LineA-final-01.fa
│   └── LineB-final-04.fa
├── output
│   ├── alignment_files/
│   ├── composite_assemblies *
│   │   └── LineA_B_composite.fa *
│   ├── ref_ranges.bed
│   ├── updated_assemblies/
│   ├── vcf_files/
│   │   vcf_files_imputed
│   │   └── LineA_B.h.vcf 
│   └── read_mappings/
└── vcf_dbs
    ├── assemblies.agc
    ├── gvcf_dataset/
    ├── hvcf_dataset/
    ├── hvcf_files/
    ├── reference/
    └── temp/
```


### Align short reads to composite FASTA

After creating the composite FASTA, we can re-align the WGS reads used
for the imputation steps to this new composite assembly. This can
be accomplished using a short-read aligner such as 
[minimap2](https://github.com/lh3/minimap2) (_recommended_) or other 
comparable software. Since this step is dependent on the choices made 
by the end-user or group, we will provide an example of how to run 
this via minimap2 in the following code block:

```shell
# Run minimap2 to align the reads with the composite genome
minimap2 -a -x sr -t 20 --secondary=no --eqx \
    -R '@RG\tID:READ_GRP1\tSM:LineA_B' \
    output/composite_assemblies/LineA_B_composite.fa \
    data/short_reads/LineA_LineB_1.fq \
    data/short_reads/LineA_LineB_2.fq | \
    samtools view -b > output/minimap2/LineA_B_composite_align.bam
```

In summary, this command is using the following `minimap2` parameters:

* `-a` - output **a**lignments in [SAM](https://samtools.github.io/hts-specs/SAMv1.pdf)
  format
* `-x sr` - alignment mode is **s**hort-**r**ead alignment
* `-t 20` - run this alignment using **20 threads** on our example
  machine (_this can be modified to leverage the full potential of
  your machine_)
* `--secondary=no` - suppress secondary alignments and report only
  the primary alignment
* `--eqx` - output the [CIGAR](https://timd.one/blog/genomics/cigar.php)
  string in "extended" format (`=` for matches and `X` for mismatches)
* `-R` - add read group information to the SAM file (useful if you
  have reads coming from multiple libraries)

This minimap2 procedure will create a SAM file that is piped (`|`)
into the [SAMtools](https://github.com/samtools/samtools) program which is converted (`view -b`) and
saved (`> output/minimap2/LineA_B_composite_align.bam`) as a BAM
format for more efficient storage and downstream procedures.


### Filter BAM files

BAM files can be filtered using SAMtools (_recommended_) or other 
SAM/BAM processing tools. We suggest running the following SAMtools
commands: 

* Fix potential issues in paired-end reads of a BAM file by adjusting
  and correcting mate information by ensuring each pair of reads has
  accurate information about its mate score (`-m`):

    ```shell
    samtools fixmate -m LineA_B_composite_align.bam LineA_B_composite_fm.bam
    ```

* Sort alignments by reference sequence coordinates which is
  necessary before deduplication and indexing:

    ```shell
    samtools sort --threads 8 -T ./tmp -o LineA_B_composite_sorted_fm.bam LineA_B_composite_fm.bam
    ```

* Identify and remove duplicate reads:

    ```shell
    samtools markdup -r LineA_B_composite_sorted_fm.bam LineA_B_composite_dedup.bam
    ```

* Exclude unmapped reads or those that are not properly paired:

    ```shell
    samtools view -F 4 -f 2 -b LineA_B_composite_dedup.bam > LineA_B_composite_filtered.bam
    ```

* Index alignments for fast random access:

    ```shell
    samtools index LineA_B_composite_filtered.bam
    ```



### Perform variant calling

Here we show an example of running DeepVariant on the merged bam file 
created above. The output of DeepVariant is a VCF file that contains 
the variants called against the composite fasta. The following 
command will run DeepVariant from a Docker container:

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


### Filter variants

To identify rare alleles, the VCF file output by the variant caller should be filtered to select variants based on 
quality scores, depth or other measures.  The following is an example of how to filter the VCF file created above.
It filters for biallelic SNPs with a depth greater than 4 and a quality score greater than 20.

```shell
bcftools view -m2 -M2 -v snps -i 'QUAL>20 && GT="1/1" && DP>4' /path/to/merged.sorted.vcf -o /path/to/deepVariantFiltered.vcf
```


### Create final VCF in haplotype coordinates

Convert the VCF file created above to haplotype coordinates, using the create-haplotype-vcf command as shown below:
  ```shell
  phg composite-to-haplotype-coords \
      --path-hvcf /path/to/imputation.hvcf \
      --variant-vcf /path/to/vcf/from/DeepVariant \ // vcf created in filtered step above
      --sample-name nameForNewSample \ // samplename to use in the haplotype coordinates vcf
      --output-file /path/to/output.vcf
  ```



## Extra materials

!!! todo
    Should we keep the following code blocks:

### Multiple libraries:

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


### Merge

* Merge the bam files into 1

Before running the variant caller, merge the bam files into a single 
BAM file, then index using the commands below:

```shell
samtools merge -b /path/to/bam.list /path/to/merged.bam
samtools sort --threads 8 -T ./tmp -o /path/to/merged.sorted.bam /path/to/merged.bam
samtools index /path/to/merged.sorted.bam
```