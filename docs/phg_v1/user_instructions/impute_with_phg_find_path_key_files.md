!!! warning "Legacy Documentation - PHG Version 1"

    This section contains documentation for **PHG version 1**, which is
    no longer actively developed. It is preserved here for archival and
    historical reference only. If you are looking to use the Practical
    Haplotype Graph, please refer to the [PHG v2 documentation](../../index.md),
    which reflects the current version of the software.

# Find Path Key Files

These two keyfiles are tab-separated text files which setup fastq alignment to the pangenome and which readMappings need to be aggregated when finding Paths.


## FastqToMappingPlugin KeyFile Specification:

The PHG will process the following columns:

| HeaderName  | Description  | Required |
|---|---|---|
| cultivar | Name of the taxon to be processed. | Yes |
| flowcell_lane | Name of the flow cell this sample came from  | Yes |
| filename | Name of the Fastq file to be processed | Yes | 
| filename2 | Second Fastq file to be processed.  If set, minimap2 will operate in paired end mode  | No |
| PlateID | The id of the Plate. | No | 


Note that cultivar + flowcell_lane + PlateID must be a unique value across all samples in the keyfile.  If this is not, the PHG will not work correctly. After running through FastqToMapping, an additional column will be added denoting the ReadMappingID values in the DB.  This is purely for reference and is not used currently in the pipeline.

## Sample File:


```
#!txt

cultivar	flowcell_lane	filename	PlateID
RecLineB1RefA1gco4_wgs	wgsFlowcell	RecLineB1RefA1gco4_R1.fastq	wgs
RecLineB1RefA1gco4_gbs	gbsFlowcell	RecLineB1RefA1gco4_R1_gbs.fastq	gbs
RecRefA1LineBgco6_wgs	wgsFlowcell	RecRefA1LineBgco6_R1.fastq	wgs
RecRefA1LineBgco6_gbs	gbsFlowcell	RecRefA1LineBgco6_R1_gbs.fastq	gbs
LineB1_wgs	wgsFlowcell	LineB1_R1.fastq	wgs
LineB1_gbs	gbsFlowcell	LineB1_R1_gbs.fastq	gbs
LineA_wgs	wgsFlowcell	LineA_R1.fastq	wgs
LineA_gbs	gbsFlowcell	LineA_R1_gbs.fastq	gbs
LineB_wgs	wgsFlowcell	LineB_R1.fastq	wgs
LineB_gbs	gbsFlowcell	LineB_R1_gbs.fastq	gbs
LineA1_wgs	wgsFlowcell	LineA1_R1.fastq	wgs
LineA1_gbs	gbsFlowcell	LineA1_R1_gbs.fastq	gbs
RecLineALineB1gco3_wgs	wgsFlowcell	RecLineALineB1gco3_R1.fastq	wgs
RecLineALineB1gco3_gbs	gbsFlowcell	RecLineALineB1gco3_R1_gbs.fastq	gbs
RecLineB1LineBgco7_wgs	wgsFlowcell	RecLineB1LineBgco7_R1.fastq	wgs
RecLineB1LineBgco7_gbs	gbsFlowcell	RecLineB1LineBgco7_R1_gbs.fastq	gbs
RefA1_wgs	wgsFlowcell	RefA1_R1.fastq	wgs
RefA1_gbs	gbsFlowcell	RefA1_R1_gbs.fastq	gbs
Ref_wgs	wgsFlowcell	Ref_R1.fastq	wgs
Ref_gbs	gbsFlowcell	Ref_R1_gbs.fastq	gbs
RecLineBLineB1gco8_wgs	wgsFlowcell	RecLineBLineB1gco8_R1.fastq	wgs
RecLineBLineB1gco8_gbs	gbsFlowcell	RecLineBLineB1gco8_R1_gbs.fastq	gbs
RecLineA1LineA1gco2_wgs	wgsFlowcell	RecLineA1LineA1gco2_R1.fastq	wgs
RecLineA1LineA1gco2_gbs	gbsFlowcell	RecLineA1LineA1gco2_R1_gbs.fastq	gbs
RecLineBLineB1gco5_wgs	wgsFlowcell	RecLineBLineB1gco5_R1.fastq	wgs
RecLineBLineB1gco5_gbs	gbsFlowcell	RecLineBLineB1gco5_R1_gbs.fastq	gbs
RecLineA1RefA1gco1_wgs	wgsFlowcell	RecLineA1RefA1gco1_R1.fastq	wgs
RecLineA1RefA1gco1_gbs	gbsFlowcell	RecLineA1RefA1gco1_R1_gbs.fastq	gbs

```

## Find Path Keyfile Specification

Note: This file will likely not be needed to create.  FastqToMappingPlugin will export a sample keyfile given its inputs.  This will probably be good enough for the majority of use cases.

The PHG will process the following columns:

| HeaderName  | Description  | Required |
|---|---|---|
| sampleName | Name of the taxon to be processed. | Yes |
| ReadMappingIds | ReadMappingIds in the DB. All readMappings need to be comma separated. | Yes |
| LikelyParents | List of taxon which the user believes are the likely parents for this sample.  Currently this is not used. | No |


## Sample File:


```
#!txt

SampleName	ReadMappingIds	LikelyParents
RecLineB1RefA1gco4_wgs	1		
RecLineB1RefA1gco4_gbs	2		
RecRefA1LineBgco6_wgs	3		
RecRefA1LineBgco6_gbs	4		
LineB1_wgs	5		
LineB1_gbs	6		
LineA_wgs	7		
LineA_gbs	8		
LineB_wgs	9		
LineB_gbs	10		
LineA1_wgs	11		
LineA1_gbs	12		
RecLineALineB1gco3_wgs	13		
RecLineALineB1gco3_gbs	14		
RecLineB1LineBgco7_wgs	15		
RecLineB1LineBgco7_gbs	16		
RefA1_wgs	17		
RefA1_gbs	18		
Ref_wgs	19		
Ref_gbs	20		
RecLineBLineB1gco8_wgs	21		
RecLineBLineB1gco8_gbs	22		
RecLineA1LineA1gco2_wgs	23		
RecLineA1LineA1gco2_gbs	24		
RecLineBLineB1gco5_wgs	25		
RecLineBLineB1gco5_gbs	26		
RecLineA1RefA1gco1_wgs	27		
RecLineA1RefA1gco1_gbs	28		

```

[Return to Step 3 pipeline version 0.0.40 or older](impute_with_phg_main.md)

[Return to Step 3 pipeline version 0.1.0 or newer](../home_variants_in_gvcf_files.md)

[Return to Wiki Home](../home.md)