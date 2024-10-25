# Variant Caller Comparisons

Currently, we have investigated 3 programs for calling variants from 
BAM alignment files:

* [DeepVariant](https://github.com/google/deepvariant)
* [Octopus](https://github.com/luntergroup/octopus)
* [Clair3](https://github.com/HKU-BAL/Clair3)

While all 3 had their merits, we found **DeepVariant** and **Octopus**
were more suited to the short-read variant calling that is often used 
with PhgV2. We found DeepVariant easier to set up than Octopus, but 
was less flexible in terms of parameters. We had difficulty pulling 
Octopus [forest models](https://luntergroup.github.io/octopus/docs/guides/filtering/forest), 
but found when running in the [individual](https://luntergroup.github.io/octopus/docs/guides/models/individual) 
calling model, they were not necessary. Additionally, we discovered 
that leveraging Docker environments for these programs was more 
efficient to set up than trying to link with a Conda-based 
environment.

Pulling from our experience and web comparisons, here is a synopsis 
of each tool:

| Feature                        | DeepVariant                           | Octopus                            | Clair3                           |
|--------------------------------|---------------------------------------|------------------------------------|----------------------------------|
| **Accuracy**                   | High (short and long reads)           | High (especially in somatic)       | High (especially in long reads)  |
| **Use Cases**                  | Germline variant calling              | Germline & somatic variant calling | Long-read variant calling        |
| **Performance**                | Best with GPUs, can be slow otherwise | Fast even without GPU              | Efficient with long reads        |
| **Ease of Use**                | Easy but less customizable            | Flexible but complex setup         | Good for long-read, easier setup |
| **Supported Data Types**       | Short & long-read data                | Short & long-read data             | Long-read focused                |
| **Computational Requirements** | High if using GPUs                    | Moderate                           | Moderate (especially with GPUs)  |


For more in-depth notes, please take a look at the proceeding
sections:

## DeepVariant

### Strengths

* **Accuracy**: DeepVariant uses deep learning to call variants, 
  which leads to high accuracy in calling both SNPs and indels.
* **Performance across platforms**: Works well across different 
  sequencing technologies like Illumina, PacBio, and Oxford Nanopore, 
  making it versatile.
* **Ease of Use**: Pre-trained models for different platforms are 
  available, reducing the need for parameter tuning or custom 
  training.
* **Automation**: It is relatively easy to set up, and its automated 
  pipeline makes variant calling straightforward.
* **Active Development**: Continuous improvements are made by Google 
  and the open-source community.

### Weaknesses

* **Computational Resources**: DeepVariant requires significant 
  computational resources, especially GPUs, to achieve optimal 
  performance. However, Buckler Lab was able to run successfully 
  without GPUs on a high performance computing machine.
* **Limited Customizability**: Since it's a pre-trained deep learning 
  model, there’s less flexibility for users who want to tweak 
  algorithms or tailor the caller to niche applications.
* **Not tailored for somatic mutation calling**: It's more optimized 
  for germline variant calling, so may not perform as well on somatic 
  mutations compared to specialized callers.

## Octopus

### Strengths

* **Versatility**: Octopus supports both germline and somatic variant 
  calling, making it suitable for diverse research applications.
* **Haplotype-aware**: It uses haplotype assembly for variant calling, 
  which improves accuracy, particularly in complex regions of the 
  genome.
* **Speed**: Compared to DeepVariant, Octopus is typically faster and 
  more computationally efficient, especially on multi-core machines.
* **Supports multiple data types**: Works with short reads and 
  long-read data, making it adaptable for mixed or evolving 
  sequencing approaches.

### Weaknesses

* **Less mature**: While Octopus is powerful, it may not have as 
  large of a user base or community as DeepVariant or Clair3, meaning 
  less documentation or community support.
* **Complex setup**: Configuration and usage may be more complex due 
  to the high degree of flexibility, and tuning parameters for 
  specific datasets could be a challenge.
* **Limited GPU support**: It doesn’t fully leverage GPUs for 
  acceleration, so it might not be as fast as GPU-accelerated tools 
  like DeepVariant when working on large datasets.

## Clair3

### Strengths
* **Accuracy for Long Reads**: Clair3 is particularly designed to 
  excel with long-read sequencing data (e.g., PacBio HiFi, Oxford 
  Nanopore), showing superior performance in variant calling for such 
  data types.
* **Phasing and Indel Detection**: Clair3 shows good performance in 
  phasing and in detecting small indels, which are often problematic 
  for short-read focused callers.
* **Deep Learning-based**: Similar to DeepVariant, Clair3 also uses a 
  neural network model, but it has been particularly optimized for 
  the complexities of long-read data.
* **Low Computational Cost**: It’s more computationally efficient 
  compared to DeepVariant, especially for long-read data, which makes 
  it faster while retaining accuracy.

### Weaknesses

* **Limited Short-Read Support**: Clair3 is designed for long-read 
  sequencing and may not be as well-suited for short-read data, which 
  can make it less versatile than DeepVariant.
* **Community Size**: Clair3 is newer and has a smaller user base 
  than DeepVariant, so there may be less widespread support or 
  examples of usage in various research contexts.
* **Requires specialized hardware**: While not as resource-heavy as 
  DeepVariant, Clair3 still benefits from high-end hardware (e.g., 
  GPUs) to maximize performance.


## Example code usage

=== "DeepVariant"

    ```shell
    # Run DeepVariant from Docker, filter for chr9
    INPUT_DIR=/workdir/lcj34/deepVariant/p39toCompositeBamFastas
    OUTPUT_DIR=/workdir/lcj34/deepVariant/p39toComposite_dvOutput
    
    docker run \
        -v ${INPUT_DIR}/:/input/ \
        -v ${OUTPUT_DIR}/:/output/ \
        google/deepvariant:1.6.1 \
        /opt/deepvariant/bin/run_deepvariant \
        --model_type=WGS \
        --ref=/input/P39wgs_composite.fa \
        --reads=/input/P39toCompositeMergedBams.bam \
        --regions "chr9" \
        --output_vcf=/output/p39Compositechr9_output.vcf.gz \
        --output_gvcf=/output/p39Compositechr9_output.g.vcf.gz \
        --intermediate_results_dir /output/intermediate_results \
        --num_shards=1
    ```

=== "Octopus"

    ```shell
    # Run Octopus from Docker, filter for chr9
    INPUT_DIR=/workdir/lcj34/octopus/input
    OUTPUT_DIR=/workdir/lcj34/octopus/output
    
    docker run \
        -v ${INPUT_DIR}/:/input/ \
        -v ${OUTPUT_DIR}/:/output/ \
        dancooke/octopus:latest \
        octopus \
        -R /input/P39.fa \
        -I /input/P39toB73MergedBams.bam \
        -T chr9 \
        --sequence-error-model PCRF.X10 \
        --forest /opt/octopus/resources/forests/germline.v0.8.0.forest \
        -o /output/p39b73chr9_octopus.vcf.gz \
        --threads=16
    ```
=== "Clair3"

    ```shell
    # Run Octopus from Docker, filter for chr9
    INPUT_DIR=/workdir/lcj34/octopus/input
    OUTPUT_DIR=/workdir/lcj34/octopus/output
    
    docker run \
        -v ${INPUT_DIR}/:/input/ \
        -v ${OUTPUT_DIR}/:/output/ \
        dancooke/octopus:latest \
        octopus \
        -R /input/P39.fa \
        -I /input/P39toB73MergedBams.bam \
        -T chr9 \
        --sequence-error-model PCRF.X10 \
        --forest /opt/octopus/resources/forests/germline.v0.8.0.forest \
        -o /output/p39b73chr9_octopus.vcf.gz \
        --threads=16
    ```
