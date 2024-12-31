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

* Create a composite FASTA file

    ```shell
    phg create-fasta-from-hvcf \
        --hvcf-file output/vcf_files_imputed/LineA_B.h.vcf \ 
        --fasta-type composite \ 
        -o output/composite_assemblies/
    ```

* Index composite FASTA file

    ```shell
    samtools faidx output/composite_assemblies/LineA_B_composite.fa
    ```

* Align short reads to composite FASTA

    ```shell
    # Run minimap2 to align the reads with the composite genome
    minimap2 -a -x sr -t 20 --secondary=no --eqx \
        output/composite_assemblies/LineA_B_composite.fa \
        data/short_reads/LineA_LineB_1.fq \
        data/short_reads/LineA_LineB_2.fq | \
        samtools view -b > output/minimap2/LineA_B_composite_align.bam
    ```

* Filter BAM files

    ```shell
    samtools fixmate -m LineA_B_composite_align.bam LineA_B_composite_fm.bam
    samtools sort --threads 8 -T ./tmp -o LineA_B_composite_sorted_fm.bam LineA_B_composite_fm.bam
    samtools markdup -r LineA_B_composite_sorted_fm.bam LineA_B_composite_dedup.bam
    samtools view -F 4 -f 2 -b LineA_B_composite_dedup.bam > LineA_B_composite_filtered.bam
    samtools index LineA_B_composite_filtered.bam
    ```

* Perform variant calling

    ```shell
    PHG_DIR=phg_v2_example/

    docker run \
        -v ${PHG_DIR}:/workdir \
        google/deepvariant:1.6.1 \
        /opt/deepvariant/bin/run_deepvariant \
        --model_type=WGS \
        --ref=/workdir/output/composite_assemblies/LineA_B_composite.fa \
        --reads=/workdir/output/minimap2/LineA_B_composite_filtered.bam \
        --output_vcf=/workdir/output/vcf_files/LineA_B_composite_reseq.vcf \
        --output_gvcf=/workdir/output/vcf_files/LineA_B_composite_reseq.g.vcf \
        --intermediate_results_dir=/workdir/intermediate_results \
        --num_shards=4
    ```


* Filter variants

    ```shell
    bcftools view -m2 -M2 -v snps \
        -i 'QUAL>20 && GT="1/1" && DP>4' \
        output/vcf_files/LineA_B_composite_reseq.vcf \
        -o output/vcf_files/LineA_B_composite_reseq_filt.vcf
    ```

* Create final VCF in haplotype coordinates

    ```shell
    phg composite-to-haplotype-coords \
        --path-hvcf output/vcf_files_imputed/LineA_B.h.vcf \
        --variant-vcf output/vcf_files/LineA_B_composite_reseq_filt.vcf \
        --sample-name LineA_B \
        --output-file output/vcf_files/LineA_B_reseq_hap_coords.vcf
    ```


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
prior pipelines) will now have a FASTA file (`LineA_B_composite.fa`) 
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
│   ├── vcf_files_imputed
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

### Index composite FASTA file
Once we create the composite FASTA file, we will need to create an
[`.fai` index](http://www.htslib.org/doc/faidx.html) file for the 
subsequent DeepVariant procedure. This indexing step will significantly
decrease lookup times into our new composite reference FASTA
assembly. To create this file, we can use `faidx` command from
`samtools`:

```shell
samtools faidx output/composite_assemblies/LineA_B_composite.fa
```

This will create an `.fai` index file in the same directory where the
composite FASTA was created. In our case, this file will be named
`LineA_B_composite.fa.fai`.

!!! note
    For the DeepVariant step, this file will **need to be in the same
    directory as your composite FASTA file**. DeepVariant will assume
    that this file and the FASTA assembly are in the same directory.
    If not, an error will occur!


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
After the BAM files have been filtered, we can now run variant
calling. Since there are many different variant callers available,
this step can become "highly opinionated" for various groups and
users. In our current tests, we have found the results produced by
[DeepVariant](https://github.com/google/deepvariant) to be of high
quality.

!!! note
    For more information on comparisons between different variant
    callers we have tested, please refer to the notes found
    [on this page](variant_comparisons.md).

Here we show an example of running DeepVariant on the prior BAM file. 
The output of DeepVariant is a VCF file that contains 
the variants called against the composite fasta made in the
[first step](#create-a-composite-fasta-file). Assuming you have
retrieved the [Docker image](https://github.com/google/deepvariant/blob/r1.6.1/docs/deepvariant-quick-start.md#get-docker-image)
for DeepVariant, the following command is an example of the
recommended parameters to use for the program:

```shell
PHG_DIR=phg_v2_example/

docker run \
    -v ${PHG_DIR}:/workdir \
    google/deepvariant:1.6.1 \
    /opt/deepvariant/bin/run_deepvariant \
    --model_type=WGS \
    --ref=/workdir/output/composite_assemblies/LineA_B_composite.fa \
    --reads=/workdir/output/minimap2/LineA_B_composite_filtered.bam \
    --output_vcf=/workdir/output/vcf_files/LineA_B_composite_reseq.vcf \
    --output_gvcf=/workdir/output/vcf_files/LineA_B_composite_reseq.g.vcf \
    --intermediate_results_dir=/workdir/intermediate_results \
    --num_shards=4
```

If you have limited experience in what Docker and/or DeepVariant 
commands are, here we explain in further detail the rationale for
each of these parameters:

**Docker commands**:

* `docker run`: base Docker command to run a container
* `-v ${PHG_DIR}:/workdir`: This mounts the **host system's** (i.e., 
  your machine) working PHG directory (specified as the variable, `PHG_DIR`) 
  to the `/workdir` directory inside the container. This would be
  the path on your machine that contains the necessary files 
  (_discussed below_). Here, I am mounting it to the "example" PHG
  directory that was used previously.

**DeepVariant commands**:

* `google/deepvariant:1.6.1`: This specifies the version of the
  DeepVariant Docker image to use (`1.6.1` in this case since this is
  the latest version as of writing this documentation).
* `/opt/deepvariant/bin/run_deepvariant`: This specifies the path
  inside the Docker container to the DeepVariant executable script
  that orchestrates the variant calling pipeline.

**DeepVariant input parameters**:

* `--model_type=WGS`: Specifies the model type to use for variant
  calling. In this case, `WGS` stands for "**W**hole **G**enome 
  **S**equencing".
* `--ref=/workdir/output/composite_assemblies/LineA_B_composite.fa`:
  The reference genome in FASTA format that DeepVariant will use for
  alignment and variant calling. The reference file in this case
  **will be the composite FASTA file created earlier**.
* `--reads=/workdir/output/minimap2/LineA_B_composite_filtered.bam`:
  Specifies the prior **sorted and indexed BAM file** containing the
  aligned sequencing reads.
* `--output_vcf=/workdir/output/vcf_files/LineA_B_composite_reseq.vcf`:
  The output VCF file path where DeepVariant will write the called
  variants.
* `--output_gvcf=/workdir/output/vcf_files/LineA_B_composite_reseq.g.vcf`:
  The output gVCF file path.
* `--intermediate_results_dir=/workdir/intermediate_results`: The
  directory where DeepVariant will write intermediate files generated
  during the analysis. These files are useful for debugging or
  optimizing performance, but can be deleted after run if not needed.
* `--num_shards=4`: This specifies the number of CPU cores or threads
  to parallelize the job across. In this case, we are using `4` CPU
  threads. The higher number of shards, the faster the job, but will
  require more computational resources.



### Filter variants
To identify rare alleles, the VCF file output by the variant caller 
should be filtered to select variants based on quality scores, depth, 
or other measures. The following is an example of how to use 
[bcftools](https://github.com/samtools/bcftools) to filter the VCF 
file created above:

```shell
bcftools view -m2 -M2 -v snps \
    -i 'QUAL>20 && GT="1/1" && DP>4' \
    output/vcf_files/LineA_B_composite_reseq.vcf \
    -o output/vcf_files/LineA_B_composite_reseq_filt.vcf
```
In summary, this command is using the following bcftools parameters:

* `-m2`: Specifies that only **diploid variants** should be included.
  The value `2` represents the ploidy level (diploid).
* `-M2`: Ensures that variants with more than two alleles are
  excluded. Using this in conjunction with `-m2` ensures that only
  **biallelic diploid variants** are retained.
* `-v snps`: This option restricts the output to only **SNPs**.
* `-i 'QUAL>20 && GT="1/1" && DP>4'`: This specifies the filtering
  conditions in which variants have the qualities:
    + `QUAL>20`: only include variants with a **quality score greater
      than 20**. The quality score indicates the confidence of the
      variant call, and this filter ensures that only high-confidence
      variants are retained.
    + `GT="1/1"`: this ensures that only **homozygous alternate
      variants** are included (where the genotype is `1/1`, 
      indicating) both alleles are the alternate allele).
    + `DP>4`: retain only variants with a **depth of coverage
      greater than 4**. The depth of coverage (`DP`) indicates the
      number of reads covering the variant position.
* `output/vcf_files/LineA_B_composite_reseq.vcf`: the
  input VCF file that contains the variants to be filtered.
* `-o output/vcf_files/LineA_B_composite_reseq_filt.vcf`: 
  the output file where the filtered variants will be saved.



### Create final VCF in haplotype coordinates
Now that we have retained "high quality" variants, we can finally
convert the coordinates from the composite FASTA file to coordinates
found in the haplotypes within PHG population. To perform this, we
can use a PHGv2 command called `composite-to-haplotype-coords`:

!!! note
    The VCF coordinates created in this step will **not** be based on
    reference coordinates or coordinates relative to any single
    assembly. They will be relative to the haplotypes used to make
    the composite assembly!

```shell
phg composite-to-haplotype-coords \
    --path-hvcf output/vcf_files_imputed/LineA_B.h.vcf \
    --variant-vcf output/vcf_files/LineA_B_composite_reseq_filt.vcf \
    --sample-name LineA_B \
    --output-file output/vcf_files/LineA_B_reseq_hap_coords.vcf
```

In summary, this command takes 4 parameters:

* `--path-hvcf`: path to the imputed hVCF file created during the
  [imputation steps](imputation.md).
* `--variant-vcf`: path to the [filtered VCF file](#filter-variants) 
  created by the prior variant caller (e.g., DeepVariant).
* `--sample-name`: sample ID to include in the VCF file.
* `--output-file`: path and name of coordinate converted VCF file.



## Miscellaneous

In the following sections, we present possible alternate workflows
that may arise due to the attributes of your starting data:

### WGS reads coming from multiple libraries:

In certain cases, you may have reads that come from multiple flowcell
or sequencing libraries. To align short reads to the composite, you
can add read grouping information (`-R`) to minimap2. In the
following example, we have three WGS libraries for the sample `P39`:

```shell
# Run minimap2 to align the reads with the composite genome
FASTA_DIR=/workdir/lcj34/phg_v2/testing/fastas
FASTQ_DIR=/workdir/lcj34/phg_v2/testing/fastqFiles
WORK_DIR=/workdir/lcj34/phg_v2/testing

minimap2 -a -x sr -t 20 --secondary=no --eqx -R '@RG\tID:READ_GRP1\tSM:P39' \
    ${FASTA_DIR}/P39wgs_composite.fa \
    ${FASTQ_DIR}/P39/P39_HBEN2ADXX_GTGAAA_R1.fastq.gz \
    ${FASTQ_DIR}/P39/P39_HBEN2ADXX_GTGAAA_R2.fastq.gz | \
    samtools view -b > ${WORK_DIR}/p39top39composite_bams/P39_HBEN2ADXX_GTGAAA.ori.bam
    
minimap2 -a -x sr -t 20 --secondary=no --eqx -R '@RG\tID:READ_GRP2\tSM:P39' \
    ${FASTA_DIR}/P39wgs_composite.fa \
    ${FASTQ_DIR}/P39/P39__HFFGKADXX_GTGAAA_R1.fastq.gz \
    ${FASTQ_DIR}/P39/P39__HFFGKADXX_GTGAAA_R2.fastq.gz | \
    samtools view -b > ${WORK_DIR}/p39top39composite_bams/P39_HFFGKADXX_GTGAAA.ori.bam
    
minimap2 -a -x sr -t 20 --secondary=no --eqx -R '@RG\tID:READ_GRP3\tSM:P39'
    ${FASTA_DIR}/P39wgs_composite.fa \
    ${FASTQ_DIR}/P39/P39_HL5WNCCXX_GTGAAA_R1.fastq.gz \
    ${FASTQ_DIR}/P39/P39_HL5WNCCXX_GTGAAA_R2.fastq.gz | \
    samtools view -b > ${WORK_DIR}/p39top39composite_bams/P39_HL5WNCCXX_GTGAAA.ori.bam
```

In the above example, read groupings are added using the `-R` 
parameter found in minimap2. Each WGS library will get its own
read group:

* `ID:READ_GRP1` - `ID:READ_GRP3`: specifies that the read group
  identifier is `READ_GRP1` through `READ_GRP3`
* `SM:P39`: The sample name. In this case, it is `P39`

After alignment, the BAM files can be merged into one using the 
`merge` command from samtools:

```
samtools merge -b /path/to/bam.list /path/to/merged.bam
```

* `-b`: A list of BAM files to merge. This would be a plain-text
  list of file paths. For example:

    ```
    /path/to/sample1.bam
    /path/to/sample2.bam
    /path/to/sample3.bam
    ```

After merging, the BAM file can be sorted and indexed using the
prior guidelines:

```
samtools sort --threads 8 -T ./tmp -o /path/to/merged.sorted.bam /path/to/merged.bam
samtools index /path/to/merged.sorted.bam
```
