# PHG v2 Terminology

 Term                | Definition                                                                                          
---------------------|-----------------------------------------------------------------------------------------------------
 Reference genome    | genome used for initial alignment and base coordinates                                              
 Reference range     | segment of the reference genome                                                                     
 Haplotype           | sequence of part of an individual chromosome with its start and stop defined by the reference range 
 Reference Haplotype | haplotype from the reference genome                                                                 
 Alternate genome    | high quality genomes used to identify alternate haplotypes                                          
 Alternate haplotype | haplotype derived from a genome assembly                                                            
 Composite genome    | inferred genome based on its composite set of alternate and reference haplotypes                    
 Haplotype ID        | MD5 checksum for the haplotype sequence                                                             
 Sample              | genotype (haploid or diploid or higher), taxon, individual                                          
 Path                | phased set of haplotype ids through the pangenome graph                                             


 File Type | Use        
-----------|------------
 VCF       |  
 h.VCF     | a VCF file 
 g.VCF     | 
 BCF       | 
 BAM       |
 fasta     | 
 agc       |
 bed       | 
 tiledb    |
 h.tileDB  |
 g.tileDB  | 
