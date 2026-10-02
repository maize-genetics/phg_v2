Process for imputing from a VCF file
This is beta code that has not been fully tested

The imputation steps are the same as imputing from a fastq file except instead of calling map-reads you should call map-reads-vcf

* phg map-reads-vcf \
--hvcf-dir /my/hvcf/dir \
--index /my/index/dir/myindex.fmd \
--key-file /my/path/keyfile \
--output-dir /my/mapping/dir

Example commands

* ./phg map-reads-vcf \
--hvcf-dir output/vcf_files \
--index output/index_files/myindex.fmd \
--key-file data/key_files/read_mapping_data.txt \
--output-dir output/read_mappings

--hvcf-dir - the directory containing the hVCF files.

--key-file - a tab-delimited list of files
In the above example, I have made a keyfile and placed it in a subdirectory under the data folder called key_files.
My example keyfile would look like the following:
sampleName  filename  filename2
project reference/iwgsc_refseqv2.1_assembly.fa  vcf/TCAP90K_NAMparents_panel-strand.vcf

--output-dir - the directory to place the read-mapping files.

--index - the ropeBWT3 index file created by the ropebwt-index command.

* Detailed explanation

The VCF file should be aligned to the reference genome so that the ref column matches the reference assembly.
Multi-allelic SNPs are not supported at this time
A fastq file is generated for each sample in the VCF file and saved into a fastq directory under the output directory. The length of the fastq reads has been hard coded to 200 bases which works well with wheat. 
Then ropebwt3 is called to map each of the sample fastq files.

* Performance

Converting a 90K Illumina array VCF file for 150 samples to fastq files takes 1-2 minutes.
Mapping the fastq file takes about 1 minute per sample.

