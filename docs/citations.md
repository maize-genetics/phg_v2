On this page, we highlight articles referencing the use of the
Practical Haplotype Graph system. This page is separatated into
three sections:

* [**Core Publications**](#core-publications) - publications from the
  [Buckler Lab](https://www.maizegenetics.net) detailing the use and
  creation of the PHG in various plant systems.
* [**Publications leveraging PHG resources**](#publications-leveraging-phg-resources) - publications that
  utilize the PHG or resources generated from the PHG in various
  plant systems.
* [**Publications discussing the PHG**](#publications-discussing-the-phg) - publications that discuss
  the PHG in various article formats (e.g., book chapters, reviews,
  etc.)


## Core Publications


!!! citation "The Practical Haplotype Graph, a platform for storing and using pangenomes for imputation"

    === ":octicons-telescope-16: Summary"

        * Details the overview, architecture, and use-cases of **version 1** of the Practical
          Haplotype Graph
        * Highlights data compression and imputation capabilities with low coverage sequence data
        * Highlights expandability with open source tools such as the [Breeding API](https://www.brapi.org)
          and [rPHG](https://rphg2.maizegenetics.net)

    === ":material-sprout: Plant Species"

        $\quad$ **Maize ([_Zea mays_](https://en.wikipedia.org/wiki/Maize))**

    === ":material-account-school: Citation"

        Bradbury, P J and Casstevens, T and Jensen, S E and Johnson, L C and Miller, Z R and Monier, B and Romay, M C and Song, B and Buckler, E S (2022). **The Practical Haplotype Graph, a platform for storing and using pangenomes for imputation.** *Bioinformatics*. DOI: [10.1093/bioinformatics/btac410](https://doi.org/10.1093/bioinformatics/btac410)

    === ":simple-latex: Citation (BibTeX)"

        ```
        @article{10.1093/bioinformatics/btac410,
            author = {Bradbury, P J and Casstevens, T and Jensen, S E and Johnson, L C and Miller, Z R and Monier, B and Romay, M C and Song, B and Buckler, E S},
            title = "{The Practical Haplotype Graph, a platform for storing and using pangenomes for imputation}",
            journal = {Bioinformatics},
            volume = {38},
            number = {15},
            pages = {3698-3702},
            year = {2022},
            month = {06},
            abstract = "{Pangenomes provide novel insights for population and quantitative genetics, genomics and breeding not available from studying a single reference genome. Instead, a species is better represented by a pangenome or collection of genomes. Unfortunately, managing and using pangenomes for genomically diverse species is computationally and practically challenging. We developed a trellis graph representation anchored to the reference genome that represents most pangenomes well and can be used to impute complete genomes from low density sequence or variant data.The Practical Haplotype Graph (PHG) is a pangenome pipeline, database (PostGRES \&amp; SQLite), data model (Java, Kotlin or R) and Breeding API (BrAPI) web service. The PHG has already been able to accurately represent diversity in four major crops including maize, one of the most genomically diverse species, with up to 1000-fold data compression. Using simulated data, we show that, at even 0.1× coverage, with appropriate reads and sequence alignment, imputation results in extremely accurate haplotype reconstruction. The PHG is a platform and environment for the understanding and application of genomic diversity.All resources listed here are freely available. The PHG Docker used to generate the simulation results is https://hub.docker.com/ as maizegenetics/phg:0.0.27. PHG source code is at https://bitbucket.org/bucklerlab/practicalhaplotypegraph/src/master/. The code used for the analysis of simulated data is at https://bitbucket.org/bucklerlab/phg-manuscript/src/master/. The PHG database of NAM parent haplotypes is in the CyVerse data store (https://de.cyverse.org/de/) and named/iplant/home/shared/panzea/panGenome/PHG\_db\_maize/phg\_v5Assemblies\_20200608.db.Supplementary data are available at Bioinformatics online.}",
            issn = {1367-4803},
            doi = {10.1093/bioinformatics/btac410},
            url = {https://doi.org/10.1093/bioinformatics/btac410},
            eprint = {https://academic.oup.com/bioinformatics/article-pdf/38/15/3698/49884088/btac410.pdf},
        }
        ```

!!! citation "Genome-wide imputation using the practical haplotype graph in the heterozygous crop cassava"

    === ":octicons-telescope-16: Summary"

        * Using the PHG improves genotype imputation accuracy in cassava, especially for rare and heterozygous alleles.
        * Using IBD (Identity-by-Descent) to sample phased haplotypes leads to more accurate imputation than computational methods.
        * High genomic prediction accuracy from low-depth sequencing, supporting more efficient plant breeding.

    === ":material-sprout: Plant Species"

        $\quad$ **Cassava ([_Manihot esculenta_](https://en.wikipedia.org/wiki/Cassava))**

    === ":material-account-school: Citation"

        Long, Evan M and Bradbury, Peter J and Romay, M Cinta and Buckler, Edward S and Robbins, Kelly R (2021). **Genome-wide imputation using the practical haplotype graph in the heterozygous crop cassava.** *G3 Genes|Genomes|Genetics*. DOI: [10.1093/g3journal/jkab383](https://doi.org/10.1093/g3journal/jkab383)

    === ":simple-latex: Citation (BibTeX)"

        ```
        @article{10.1093/g3journal/jkab383,
            author = {Long, Evan M and Bradbury, Peter J and Romay, M Cinta and Buckler, Edward S and Robbins, Kelly R},
            title = "{Genome-wide imputation using the practical haplotype graph in the heterozygous crop cassava}",
            journal = {G3 Genes|Genomes|Genetics},
            volume = {12},
            number = {1},
            pages = {jkab383},
            year = {2021},
            month = {11},
            abstract = "{Genomic applications such as genomic selection and genome-wide association have become increasingly common since the advent of genome sequencing. The cost of sequencing has decreased in the past two decades; however, genotyping costs are still prohibitive to gathering large datasets for these genomic applications, especially in nonmodel species where resources are less abundant. Genotype imputation makes it possible to infer whole-genome information from limited input data, making large sampling for genomic applications more feasible. Imputation becomes increasingly difficult in heterozygous species where haplotypes must be phased. The practical haplotype graph (PHG) is a recently developed tool that can accurately impute genotypes, using a reference panel of haplotypes. We showcase the ability of the PHG to impute genomic information in the highly heterozygous crop cassava (Manihot esculenta). Accurately phased haplotypes were sampled from runs of homozygosity across a diverse panel of individuals to populate PHG, which proved more accurate than relying on computational phasing methods. The PHG achieved high imputation accuracy, using sparse skim-sequencing input, which translated to substantial genomic prediction accuracy in cross-validation testing. The PHG showed improved imputation accuracy, compared to a standard imputation tool Beagle, especially in predicting rare alleles.}",
            issn = {2160-1836},
            doi = {10.1093/g3journal/jkab383},
            url = {https://doi.org/10.1093/g3journal/jkab383},
            eprint = {https://academic.oup.com/g3journal/article-pdf/12/1/jkab383/42340789/jkab383.pdf},
        }
        ```


!!! citation "A Maize Practical Haplotype Graph Leverages Diverse NAM Assemblies"

    === ":octicons-telescope-16: Summary"

        * The PHG was used to create a pangenome representation of maize by leveraging 27 high-quality
          genome assemblies, capturing substantial structural diversity in maize.
        * Successfully imputed genotypes with high accuracy. For related lines, the accuracy was over
          99%, while for unrelated lines, it reached up to 95% using whole-genome sequencing (WGS).
        * Highly space-efficient, storing data in a manner **30,000 times smaller** than typical genotype
          files, making it a powerful tool for handling large-scale genomic data.

    === ":material-sprout: Plant Species"

        $\quad$ **Maize ([_Zea mays_](https://en.wikipedia.org/wiki/Maize))**

    === ":material-account-school: Citation"

        Franco, Jose A. Valdes and Gage, Joseph L. and Bradbury, Peter J. and Johnson, Lynn C. and Miller, Zachary R. and Buckler, Edward S. and Romay, M. Cinta (2020). **A Maize Practical Haplotype Graph Leverages Diverse NAM Assemblies.** *bioRxiv*. DOI: [10.1101/2020.08.31.268425](https://doi.org/10.1101/2020.08.31.268425)

    === ":simple-latex: Citation (BibTeX)"

        ```
        @article {Franco2020.08.31.268425,
            author = {Franco, Jose A. Valdes and Gage, Joseph L. and Bradbury, Peter J. and Johnson, Lynn C. and Miller, Zachary R. and Buckler, Edward S. and Romay, M. Cinta},
            title = {A Maize Practical Haplotype Graph Leverages Diverse NAM Assemblies},
            elocation-id = {2020.08.31.268425},
            year = {2020},
            doi = {10.1101/2020.08.31.268425},
            publisher = {Cold Spring Harbor Laboratory},
            abstract = {As a result of millions of years of transposon activity, multiple rounds of ancient polyploidization, and large populations that preserve diversity, maize has an extremely structurally diverse genome, evidenced by high-quality genome assemblies that capture substantial levels of both tropical and temperate diversity. We generated a pangenome representation (the Practical Haplotype Graph, PHG) of these assemblies in a database, representing the pangenome haplotype diversity and providing an initial estimate of structural diversity. We leveraged the pangenome to accurately impute haplotypes and genotypes of taxa using various kinds of sequence data, ranging from WGS to extremely-low coverage GBS. We imputed the genotypes of the recombinant inbred lines of the NAM population with over 99\% mean accuracy, while unrelated germplasm attained a mean imputation accuracy of 92 or 95\% when using GBS or WGS data, respectively. Most of the imputation errors occur in haplotypes within European or tropical germplasm, which have yet to be represented in the maize PHG database. Also, the PHG stores the imputation data in a 30,000-fold more space-efficient manner than a standard genotype file, which is a key improvement when dealing with large scale data.Competing Interest StatementThe authors have declared no competing interest.},
            URL = {https://www.biorxiv.org/content/early/2020/09/28/2020.08.31.268425},
            eprint = {https://www.biorxiv.org/content/early/2020/09/28/2020.08.31.268425.full.pdf},
            journal = {bioRxiv}
        }
        ```

!!! citation "A sorghum practical haplotype graph facilitates genome‐wide imputation and cost‐effective genomic prediction"

    === ":octicons-telescope-16: Summary"

        * The PHG successfully imputed genome-wide variants from only 0.01x
          sequence coverage with minimal loss in accuracy, making it highly
          cost-effective for genomic prediction.
        * The PHG can combine data from different sequencing platforms,
          unifying genotype calls, and facilitating easier comparison across
          breeding programs.
        * The PHG enables cost-effective genomic selection, with prediction
          accuracies comparable to higher-cost methods like
          genotyping-by-sequencing (GBS) or rhAmpSeq, thus supporting faster
          and cheaper breeding programs.

    === ":material-sprout: Plant Species"

        $\quad$ **Sorghum ([_Sorghum bicolor_](https://en.wikipedia.org/wiki/Sorghum))**

    === ":material-account-school: Citation"

        Jensen, Sarah E. and Charles, Jean Rigaud and Muleta, Kebede and Bradbury, Peter J. and Casstevens, Terry and Deshpande, Santosh P. and Gore, Michael A. and Gupta, Rajeev and Ilut, Daniel C. and Johnson, Lynn and Lozano, Roberto and Miller, Zachary and Ramu, Punna and Rathore, Abhishek and Romay, M. Cinta and Upadhyaya, Hari D. and Varshney, Rajeev K. and Morris, Geoffrey P. and Pressoir, Gael and Buckler, Edward S. and Ramstein, Guillaume P. (2020). **A sorghum practical haplotype graph facilitates genome‐wide imputation and cost‐effective genomic prediction.** *The Plant Genome*. DOI: [10.1002/tpg2.20009](https://doi.org/10.1002/tpg2.20009)

    === ":simple-latex: Citation (BibTeX)"

        ```
        @article{https://doi.org/10.1002/tpg2.20009,
            author = {Jensen, Sarah E. and Charles, Jean Rigaud and Muleta, Kebede and Bradbury, Peter J. and Casstevens, Terry and Deshpande, Santosh P. and Gore, Michael A. and Gupta, Rajeev and Ilut, Daniel C. and Johnson, Lynn and Lozano, Roberto and Miller, Zachary and Ramu, Punna and Rathore, Abhishek and Romay, M. Cinta and Upadhyaya, Hari D. and Varshney, Rajeev K. and Morris, Geoffrey P. and Pressoir, Gael and Buckler, Edward S. and Ramstein, Guillaume P.},
            title = {A sorghum practical haplotype graph facilitates genome-wide imputation and cost-effective genomic prediction},
            journal = {The Plant Genome},
            volume = {13},
            number = {1},
            pages = {e20009},
            doi = {https://doi.org/10.1002/tpg2.20009},
            url = {https://acsess.onlinelibrary.wiley.com/doi/abs/10.1002/tpg2.20009},
            eprint = {https://acsess.onlinelibrary.wiley.com/doi/pdf/10.1002/tpg2.20009},
            abstract = {Abstract Successful management and utilization of increasingly large genomic datasets is essential for breeding programs to accelerate cultivar development. To help with this, we developed a Sorghum bicolor Practical Haplotype Graph (PHG) pangenome database that stores haplotypes and variant information. We developed two PHGs in sorghum that were used to identify genome-wide variants for 24 founders of the Chibas sorghum breeding program from 0.01x sequence coverage. The PHG called single nucleotide polymorphisms (SNPs) with 5.9\% error at 0.01x coverage—only 3\% higher than PHG error when calling SNPs from 8x coverage sequence. Additionally, 207 progenies from the Chibas genomic selection (GS) training population were sequenced and processed through the PHG. Missing genotypes were imputed from PHG parental haplotypes and used for genomic prediction. Mean prediction accuracies with PHG SNP calls range from .57–.73 and are similar to prediction accuracies obtained with genotyping-by-sequencing or targeted amplicon sequencing (rhAmpSeq) markers. This study demonstrates the use of a sorghum PHG to impute SNPs from low-coverage sequence data and shows that the PHG can unify genotype calls across multiple sequencing platforms. By reducing input sequence requirements, the PHG can decrease the cost of genotyping, make GS more feasible, and facilitate larger breeding populations. Our results demonstrate that the PHG is a useful research and breeding tool that maintains variant information from a diverse group of taxa, stores sequence data in a condensed but readily accessible format, unifies genotypes across genotyping platforms, and provides a cost-effective option for genomic selection.},
            year = {2020}
        }
        ```

## Publications leveraging PHG resources

=== "Cassava"

    * Long, Evan M and Bradbury, Peter J and Romay, M Cinta and Buckler, Edward S and Robbins, Kelly R (2021). **Genome-wide imputation using the practical haplotype graph in the heterozygous crop cassava.** *G3 Genes|Genomes|Genetics*. DOI: [10.1093/g3journal/jkab383](https://doi.org/10.1093/g3journal/jkab383)

=== "Maize"

    * Fernandes, Igor K. and Vieira, Caio C. and Dias, Kaio O. G. and Fernandes, Samuel B. (2024). **Using machine learning to integrate genetic and environmental data to model genotype-by-environment interactions.** *bioRxiv*. DOI: [10.1101/2024.02.08.579534](https://doi.org/10.1101/2024.02.08.579534)
    * Stitzer, Michelle C. and Khaipho-Burch, Merritt B. and Hudson, Asher I. and Song, Baoxing and Valdez-Franco, Jose Arcadio and Ramstein, Guillaume and Feschotte, Cedric and Buckler, Edward S. (2023). **Transposable element abundance subtly contributes to lower fitness in maize.** *bioRxiv*. DOI: [10.1101/2023.09.18.557618](https://doi.org/10.1101/2023.09.18.557618)
    * Lima, Dayane Cristina and Aviles, Alejandro Castro and Alpers, Ryan Timothy and McFarland, Bridget A. and Kaeppler, Shawn and Ertl, David and Romay, Maria Cinta and Gage, Joseph L. and Holland, James and Beissinger, Timothy and Bohn, Martin and Buckler, Edward and Edwards, Jode and Flint-Garcia, Sherry and Hirsch, Candice N. and Hood, Elizabeth and Hooker, David C. and Knoll, Joseph E. and Kolkman, Judith M. and Liu, Sanzhen and McKay, John and Minyo, Richard and Moreta, Danilo E. and Murray, Seth C. and Nelson, Rebecca and Schnable, James C. and Sekhon, Rajandeep S. and Singh, Maninder P. and Thomison, Peter and Thompson, Addie and Tuinstra, Mitchell and Wallace, Jason and Washburn, Jacob D. and Weldekidan, Teclemariam and Wisser, Randall J. and Xu, Wenwei and de Leon, Natalia (2023). **2018–2019 field seasons of the Maize Genomes to Fields (G2F) G x E project.** *BMC Genomic Data*. DOI: [10.1186/s12863-023-01129-2](https://doi.org/10.1186/s12863-023-01129-2)
    * Giri, Anju and Khaipho-Burch, Merritt and Buckler, Edward S. and Ramstein, Guillaume P. (2021). **Haplotype Associated RNA Expression (HARE) Improves Prediction of Complex Traits in Maize.** *PLOS Genetics*. DOI: [10.1371/journal.pgen.1009568](https://doi.org/10.1371/journal.pgen.1009568)
    * Franco, Jose A. Valdes and Gage, Joseph L. and Bradbury, Peter J. and Johnson, Lynn C. and Miller, Zachary R. and Buckler, Edward S. and Romay, M. Cinta (2020). **A Maize Practical Haplotype Graph Leverages Diverse NAM Assemblies.** *bioRxiv*. DOI: [10.1101/2020.08.31.268425](https://doi.org/10.1101/2020.08.31.268425)

=== "Sorghum"

    * Somegowda, Vinutha Kanuganhalli and Diwakar Reddy, S.E. and Gaddameedi, Anil and Kiranmayee, K.N.S. Usha and Naravula, Jalaja and Kavi Kishor, P.B. and Penna, Suprasanna (2024). **Genomics breeding approaches for developing Sorghum bicolor lines with stress resilience and other agronomic traits.** *Current Plant Biology*. DOI: [10.1016/j.cpb.2023.100314](https://doi.org/10.1016/j.cpb.2023.100314)
    * Jensen, Sarah E. and Charles, Jean Rigaud and Muleta, Kebede and Bradbury, Peter J. and Casstevens, Terry and Deshpande, Santosh P. and Gore, Michael A. and Gupta, Rajeev and Ilut, Daniel C. and Johnson, Lynn and Lozano, Roberto and Miller, Zachary and Ramu, Punna and Rathore, Abhishek and Romay, M. Cinta and Upadhyaya, Hari D. and Varshney, Rajeev K. and Morris, Geoffrey P. and Pressoir, Gael and Buckler, Edward S. and Ramstein, Guillaume P. (2020). **A sorghum practical haplotype graph facilitates genome‐wide imputation and cost‐effective genomic prediction.** *The Plant Genome*. DOI: [10.1002/tpg2.20009](https://doi.org/10.1002/tpg2.20009)

=== "Wheat"

    * Wang, Hongliang and Bernardo, Amy and St. Amand, Paul and Bai, Guihua and Bowden, Robert L. and Guttieri, Mary J. and Jordan, Katherine W. (2023). **Skim exome capture genotyping in wheat.** *The Plant Genome*. DOI: [10.1002/tpg2.20381](https://doi.org/10.1002/tpg2.20381)
    * Jordan, Katherine W and Bradbury, Peter J and Miller, Zachary R and Nyine, Moses and He, Fei and Fraser, Max and Anderson, Jim and Mason, Esten and Katz, Andrew and Pearce, Stephen and Carter, Arron H and Prather, Samuel and Pumphrey, Michael and Chen, Jianli and Cook, Jason and Liu, Shuyu and Rudd, Jackie C and Wang, Zhen and Chu, Chenggen and Ibrahim, Amir M H and Turkus, Jonathan and Olson, Eric and Nagarajan, Ragupathi and Carver, Brett and Yan, Liuling and Taagen, Ellie and Sorrells, Mark and Ward, Brian and Ren, Jie and Akhunova, Alina and Bai, Guihua and Bowden, Robert and Fiedler, Jason and Faris, Justin and Dubcovsky, Jorge and Guttieri, Mary and Brown-Guedira, Gina and Buckler, Ed and Jannink, Jean-Luc and Akhunov, Eduard D (2021). **Development of the Wheat Practical Haplotype Graph database as a resource for genotyping data storage and genotype imputation.** *G3 Genes|Genomes|Genetics*. DOI: [10.1093/g3journal/jkab390](https://doi.org/10.1093/g3journal/jkab390)


## Publications discussing the PHG

=== "General Discussion"

    * **Song, Baoxing and Buckler, Edward S. and Stitzer, Michelle C.** (2024). New whole-genome alignment tools are needed for tapping into plant diversity. *Trends in Plant Science*. DOI: [10.1016/j.tplants.2023.08.013](https://doi.org/10.1016/j.tplants.2023.08.013)
    * **Coletta, Rafael Della and Fernandes, Samuel B. and Monnahan, Patrick J. and Mikel, Mark A. and Bohn, Martin O. and Lipka, Alexander E. and Hirsch, Candice N.** (2023). Importance of genetic architecture in marker selection decisions for genomic prediction. *bioRxiv*. DOI: [10.1101/2023.02.28.530521](https://doi.org/10.1101/2023.02.28.530521)
    * **da Costa Lima Moraes, Aline and Mollinari, Marcelo and Ferreira, Rebecca Caroline Ulbricht and Aono, Alexandre and de Castro Lara, Letícia Aparecida and Pessoa-Filho, Marco and Barrios, Sanzio Carvalho Lima and Garcia, Antonio Augusto Franco and do Valle, Cacilda Borges and de Souza, Anete Pereira and Vigna, Bianca Baccili Zanotto** (2023). Advances in genomic characterization ofUrochloa humidicola: exploring polyploid inheritance and apomixis. *bioRxiv*. DOI: [10.1101/2023.08.31.555743](https://doi.org/10.1101/2023.08.31.555743)
    * **Ruperao, Pradeep and Gandham, Prasad and Odeny, Damaris A. and Mayes, Sean and Selvanayagam, Sivasubramani and Thirunavukkarasu, Nepolean and Das, Roma R. and Srikanda, Manasa and Gandhi, Harish and Habyarimana, Ephrem and Manyasa, Eric and Nebie, Baloua and Deshpande, Santosh P. and Rathore, Abhishek** (2023). Exploring the sorghum race level diversity utilizing 272 sorghum accessions genomic resources. *Frontiers in Plant Science*. DOI: [10.3389/fpls.2023.1143512](https://doi.org/10.3389/fpls.2023.1143512)
    * **Brown, Pat J.** (2023). Haplotyping interspecific hybrids by dual alignment to both parental genomes. *The Plant Genome*. DOI: [10.1002/tpg2.20324](https://doi.org/10.1002/tpg2.20324)
    * **Aylward, Anthony J. and Petrus, Semar and Mamerto, Allen and Hartwick, Nolan T. and Michael, Todd P.** (2023). PanKmer:k-mer based and reference-free pangenome analysis. *bioRxiv*. DOI: [10.1101/2023.03.31.535143](https://doi.org/10.1101/2023.03.31.535143)
    * **Tanaka, Ryokei and Wu, Di and Li, Xiaowei and Tibbs‐Cortes, Laura E. and Wood, Joshua C. and Magallanes‐Lundback, Maria and Bornowski, Nolan and Hamilton, John P. and Vaillancourt, Brieanne and Li, Xianran and Deason, Nicholas T. and Schoenbaum, Gregory R. and Buell, C. Robin and DellaPenna, Dean and Yu, Jianming and Gore, Michael A.** (2022). Leveraging prior biological knowledge improves prediction of tocochromanols in maize grain. *The Plant Genome*. DOI: [10.1002/tpg2.20276](https://doi.org/10.1002/tpg2.20276)
    * **Bayer, Philipp E. and Petereit, Jakob and Durant, Éloi and Monat, Cécile and Rouard, Mathieu and Hu, Haifei and Chapman, Brett and Li, Chengdao and Cheng, Shifeng and Batley, Jacqueline and Edwards, David** (2022). Wheat Panache - a pangenome graph database representing presence/absence variation across 16 bread wheat genomes. *NA*. DOI: [10.1101/2022.02.23.481560](https://doi.org/10.1101/2022.02.23.481560)
    * **Cerioli, Tommaso and Hernandez, Christopher O. and Angira, Brijesh and McCouch, Susan R. and Robbins, Kelly R. and Famoso, Adam N.** (2022). Development and validation of an optimized marker set for genomic selection in southern U.S. rice breeding programs. *The Plant Genome*. DOI: [10.1002/tpg2.20219](https://doi.org/10.1002/tpg2.20219)
    * **Tinker, Nicholas A. and Wight, Charlene P. and Bekele, Wubishet A. and Yan, Weikai and Jellen, Eric N. and Renhuldt, Nikos Tsardakas and Sirijovski, Nick and Lux, Thomas and Spannagl, Manuel and Mascher, Martin** (2022). Genome analysis in Avena sativa reveals hidden breeding barriers and opportunities for oat improvement. *Communications Biology*. DOI: [10.1038/s42003-022-03256-5](https://doi.org/10.1038/s42003-022-03256-5)
    * **Singh, Pummi and Huang, Shun-Yuan and Hernandez, Alvaro G. and Adhikari, Pragya and Jamann, Tiffany M. and Mideros, Santiago X.** (2022). Genomic regions associated with virulence in Setosphaeria turcica identified by linkage mapping in a biparental population. *Fungal Genetics and Biology*. DOI: [10.1016/j.fgb.2021.103655](https://doi.org/10.1016/j.fgb.2021.103655)
    * **Vaughn, Justin N. and Branham, Sandra E. and Abernathy, Brian and Hulse-Kemp, Amanda M. and Rivers, Adam R. and Levi, Amnon and Wechter, William P.** (2022). Graph-based pangenomics maximizes genotyping density and reveals structural impacts on fungal resistance in melon. *Nature Communications*. DOI: [10.1038/s41467-022-35621-7](https://doi.org/10.1038/s41467-022-35621-7)
    * **Schneider, Michael and Casale, Federico and Stich, Benjamin** (2022). Accurate recombination estimation from pooled genotyping and sequencing: a case study on barley. *BMC Genomics*. DOI: [10.1186/s12864-022-08701-7](https://doi.org/10.1186/s12864-022-08701-7)
    * **Wang, Shuo and Qian, Yong-Qing and Zhao, Ru-Peng and Chen, Ling-Ling and Song, Jia-Ming** (2022). Graph-based pan-genomes: increased opportunities in plant genomics. *Journal of Experimental Botany*. DOI: [10.1093/jxb/erac412](https://doi.org/10.1093/jxb/erac412)
    * **Boatwright, J. Lucas and Sapkota, Sirjan and Jin, Hongyu and Schnable, James C. and Brenton, Zachary and Boyles, Richard and Kresovich, Stephen** (2022). Sorghum Association Panel whole‐genome sequencing establishes cornerstone resource for dissecting genomic diversity. *The Plant Journal*. DOI: [10.1111/tpj.15853](https://doi.org/10.1111/tpj.15853)
    * **Zhao, Yikun and Tian, Hongli and Li, Chunhui and Yi, Hongmei and Zhang, Yunlong and Li, Xiaohui and Zhao, Han and Huo, Yongxue and Wang, Rui and Kang, Dingming and Lu, Yuncai and Liu, Zhihao and Liang, Ziyue and Xu, Liwen and Yang, Yang and Zhou, Ling and Wang, Tianyu and Zhao, Jiuran and Wang, Fengge** (2022). HTPdb and HTPtools: Exploiting maize haplotype-tag polymorphisms for germplasm resource analyses and genomics-informed breeding. *Plant Communications*. DOI: [10.1016/j.xplc.2022.100331](https://doi.org/10.1016/j.xplc.2022.100331)
    ***Schneider, Michael and Shrestha, Asis and Ballvora, Agim and Léon, Jens** (2022). High-throughput estimation of allele frequencies using combined pooled-population sequencing and haplotype-based data processing. *Plant Methods*. DOI: [10.1186/s13007-022-00852-8](https://doi.org/10.1186/s13007-022-00852-8)
    * **Strable, Josh** (2021). Developmental genetics of maize vegetative shoot architecture. *Molecular Breeding*. DOI: [10.1007/s11032-021-01208-1](https://doi.org/10.1007/s11032-021-01208-1)
    * **Reynolds, Matthew P and Lewis, Janet M and Ammar, Karim and Basnet, Bhoja R and Crespo-Herrera, Leonardo and Crossa, José and Dhugga, Kanwarpal S and Dreisigacker, Susanne and Juliana, Philomin and Karwat, Hannes and Kishii, Masahiro and Krause, Margaret R and Langridge, Peter and Lashkari, Azam and Mondal, Suchismita and Payne, Thomas and Pequeno, Diego and Pinto, Francisco and Sansaloni, Carolina and Schulthess, Urs and Singh, Ravi P and Sonder, Kai and Sukumaran, Sivakumar and Xiong, Wei and Braun, Hans J** (2021). Harnessing translational research in wheat for climate resilience. *Journal of Experimental Botany*. DOI: [10.1093/jxb/erab256](https://doi.org/10.1093/jxb/erab256)
    * **Pook, Torsten and Nemri, Adnane and Gonzalez Segovia, Eric Gerardo and Valle Torres, Daniel and Simianer, Henner and Schoen, Chris-Carolin** (2021). Increasing calling accuracy, coverage, and read-depth in sequence data by the use of haplotype blocks. *PLOS Genetics*. DOI: [10.1371/journal.pgen.1009944](https://doi.org/10.1371/journal.pgen.1009944)
    * **Li, Jeremiah H. and Mazur, Chase A. and Berisa, Tomaz and Pickrell, Joseph K.** (2021). Low-pass sequencing increases the power of GWAS and decreases measurement error of polygenic risk scores compared to genotyping arrays. *Genome Research*. DOI: [10.1101/gr.266486.120](https://doi.org/10.1101/gr.266486.120)
    * **Tao, Yongfu and Luo, Hong and Xu, Jiabao and Cruickshank, Alan and Zhao, Xianrong and Teng, Fei and Hathorn, Adrian and Wu, Xiaoyuan and Liu, Yuanming and Shatte, Tracey and Jordan, David and Jing, Haichun and Mace, Emma** (2021). Extensive variation within the pan-genome of cultivated and wild sorghum. *Nature Plants*. DOI: [10.1038/s41477-021-00925-x](https://doi.org/10.1038/s41477-021-00925-x)
    * **Wolfe, Marnin D and Chan, Ariel W and Kulakow, Peter and Rabbi, Ismail and Jannink, Jean-Luc** (2021). Genomic mating in outbred species: predicting cross usefulness with additive and total genetic covariance matrices. *Genetics*. DOI: [10.1093/genetics/iyab122](https://doi.org/10.1093/genetics/iyab122)
    * **Ruperao, Pradeep and Thirunavukkarasu, Nepolean and Gandham, Prasad and Selvanayagam, Sivasubramani and Govindaraj, Mahalingam and Nebie, Baloua and Manyasa, Eric and Gupta, Rajeev and Das, Roma Rani and Odeny, Damaris A. and Gandhi, Harish and Edwards, David and Deshpande, Santosh P. and Rathore, Abhishek** (2021). Sorghum Pan-Genome Explores the Functional Utility for Genomic-Assisted Breeding to Accelerate the Genetic Gain. *Frontiers in Plant Science*. DOI: [10.3389/fpls.2021.666342](https://doi.org/10.3389/fpls.2021.666342)
    * **Nyine, Moses and Adhikari, Elina and Clinesmith, Marshall and Aiken, Robert and Betzen, Bliss and Wang, Wei and Davidson, Dwight and Yu, Zitong and Guo, Yuanwen and He, Fei and Akhunova, Alina and Jordan, Katherine W. and Fritz, Allan K. and Akhunov, Eduard** (2021). The Haplotype-Based Analysis of Aegilops tauschii Introgression Into Hard Red Winter Wheat and Its Impact on Productivity Traits. *Frontiers in Plant Science*. DOI: [10.3389/fpls.2021.716955](https://doi.org/10.3389/fpls.2021.716955)
    * **Cerioli, Tommaso and Hernandez, Christopher and Angira, Brijesh and McCouch, Susan and Robbins, Kelly and Famoso, Adam** (2021). Development and validation of an optimized marker set for genomic selection in Southern U. S. rice breeding programs. *NA*. DOI: [10.1002/essoar.10508975.1](https://doi.org/10.1002/essoar.10508975.1)
    * **Cooper, M and Powell, O and Voss-Fels, K P and Messina, C D and Gho, C and Podlich, D W and Technow, F and Chapman, S C and Beveridge, C A and Ortiz-Barrientos, D and Hammer, G L** (2020). Modelling selection response in plant-breeding programs using crop models as mechanistic gene-to-phenotype (CGM-G2P) multi-trait link functions. *in silico Plants*. DOI: [10.1093/insilicoplants/diaa016](https://doi.org/10.1093/insilicoplants/diaa016)
    * **Torkamaneh, Davoud and Laroche, Jérôme and Valliyodan, Babu and O’Donoughue, Louise and Cober, Elroy and Rajcan, Istvan and Vilela Abdelnoor, Ricardo and Sreedasyam, Avinash and Schmutz, Jeremy and Nguyen, Henry T. and Belzile, François** (2020). Soybean (Glycine max) Haplotype Map (GmHapMap): a universal resource for soybean translational and functional genomics. *Plant Biotechnology Journal*. DOI: [10.1111/pbi.13466](https://doi.org/10.1111/pbi.13466)

=== "Book Chapters"

    * **Naik, Yogesh Dashrath and Zhao, Chuanzhi and Channale, Sonal and Nayak, Spurthi N. and Bhutia, Karma L. and Gautam, Ashish and Kumar, Rakesh and Niranjan, Vidya and Shah, Trushar M. and Mott, Richard and Punnuri, Somashekhar and Pandey, Manish K. and Wang, Xingjun and Varshney, Rajeev K. and Thudi, Mahendar** (2024). Bioinformatics for Plant Genetics and Breeding Research. *Frontiers Technologies for Crop Improvement*. DOI: [10.1007/978-981-99-4673-0_3](https://doi.org/10.1007/978-981-99-4673-0_3)
    * **Kaur, Ishveen and Relan, Ashima and Saini, Dinesh Kumar and Kaur, Gurleen and Biswas, Anju and Singh, Lovepreet and Kaur, Shivreet and Sandhu, Karansher Singh** (2023). Revisiting the Genomic Approaches in the Cereals and the Path Forward. *Smart Plant Breeding for Field Crops in Post-genomics Era*. DOI: [10.1007/978-981-19-8218-7_1](https://doi.org/10.1007/978-981-19-8218-7_1)
    * **Joukhadar, Reem and Daetwyler, Hans D.** (2022). Data Integration, Imputation, and Meta-analysis for Genome-Wide Association Studies. *Data Integration, Imputation, and Meta-analysis for Genome-Wide Association Studies*. DOI: [10.1007/978-1-0716-2237-7_11](https://doi.org/10.1007/978-1-0716-2237-7_11)
    * **Cooper, Mark and Messina, Carlos D. and Tang, Tom and Gho, Carla and Powell, Owen M. and Podlich, Dean W. and Technow, Frank and Hammer, Graeme L.** (2022). Predicting Genotype × Environment × Management (G × E × M) Interactions for the Design of Crop Improvement Strategies: Integrating Breeder, Agronomist, and Farmer Perspectives. *Plant Breeding Reviews*. DOI: [10.1002/9781119874157.ch8](https://doi.org/10.1002/9781119874157.ch8)
    * **Bajgain, Prabin and Crain, Jared L. and Cattani, Douglas J. and Larson, Steven R. and Altendorf, Kayla R. and Anderson, James A. and Crews, Timothy E. and Hu, Ying and Poland, Jesse A. and Turner, M. Kathryn and Westerbergh, Anna and DeHaan, Lee R.** (2022). Breeding Intermediate Wheatgrass for Grain Production. *Plant Breeding Reviews*. DOI: [10.1002/9781119874157.ch3](https://doi.org/10.1002/9781119874157.ch3)
    * **Ruperao, Pradeep and Gandham, Prasad and Rathore, Abhishek** (2022). Construction of Practical Haplotype Graph (PHG) with the Whole-Genome Sequence Data. *NA*. DOI: [10.1007/978-1-0716-2067-0_15](https://doi.org/10.1007/978-1-0716-2067-0_15)

=== "Reviews"

    * **Schreiber, Mona and Jayakodi, Murukarthick and Stein, Nils and Mascher, Martin** (2024). Plant pangenomes for crop improvement, biodiversity and evolution. *Nature Reviews Genetics*. DOI: [10.1038/s41576-024-00691-4](https://doi.org/10.1038/s41576-024-00691-4)
    * **Ruperao, Pradeep and Rangan, Parimalan and Shah, Trushar and Thakur, Vivek and Kalia, Sanjay and Mayes, Sean and Rathore, Abhishek** (2023). The Progression in Developing Genomic Resources for Crop Improvement. *Life*. DOI: [10.3390/life13081668](https://doi.org/10.3390/life13081668)
    * **Jha, Uday Chand and Nayyar, Harsh and von Wettberg, Eric J. B. and Naik, Yogesh Dashrath and Thudi, Mahendar and Siddique, Kadambot H. M.** (2022). Legume Pangenome: Status and Scope for Crop Improvement. *Plants*. DOI: [10.3390/plants11223041](https://doi.org/10.3390/plants11223041)
    * **Rai, Mayank and Tyagi, Wricha** (2022). Haplotype breeding for unlocking and utilizing plant genomics data. *Frontiers in Genetics*. DOI: [10.3389/fgene.2022.1006288](https://doi.org/10.3389/fgene.2022.1006288)
    * **Tay Fernandez, Cassandria G. and Nestor, Benjamin J. and Danilevicz, Monica F. and Marsh, Jacob I. and Petereit, Jakob and Bayer, Philipp E. and Batley, Jacqueline and Edwards, David** (2022). Expanding Gene-Editing Potential in Crop Improvement with Pangenomes. *International Journal of Molecular Sciences*. DOI: [10.3390/ijms23042276](https://doi.org/10.3390/ijms23042276)
    * **Tay Fernandez, Cassandria Geraldine and Nestor, Benjamin John and Danilevicz, Monica Furaste and Gill, Mitchell and Petereit, Jakob and Bayer, Philipp Emanuel and Finnegan, Patrick Michael and Batley, Jacqueline and Edwards, David** (2022). Pangenomes as a Resource to Accelerate Breeding of Under-Utilised Crop Species. *International Journal of Molecular Sciences*. DOI: [10.3390/ijms23052671](https://doi.org/10.3390/ijms23052671)
    * **Hübner, Sariel** (2022). Are we there yet? Driving the road to evolutionary graph-pangenomics. *Current Opinion in Plant Biology*. DOI: [10.1016/j.pbi.2022.102195](https://doi.org/10.1016/j.pbi.2022.102195)
    * **Jayakodi, Murukarthick and Schreiber, Mona and Stein, Nils and Mascher, Martin** (2021). Building pan-genome infrastructures for crop plants and their use in association genetics. *DNA Research*. DOI: [10.1093/dnares/dsaa030](https://doi.org/10.1093/dnares/dsaa030)
    * **Della Coletta, Rafael and Qiu, Yinjie and Ou, Shujun and Hufford, Matthew B. and Hirsch, Candice N.** (2021). How the pan-genome is changing crop genomics and improvement. *Genome Biology*. DOI: [10.1186/s13059-020-02224-8](https://doi.org/10.1186/s13059-020-02224-8)
    * **Bhat, Javaid Akhter and Yu, Deyue and Bohra, Abhishek and Ganie, Showkat Ahmad and Varshney, Rajeev K.** (2021). Features and applications of haplotypes in crop breeding. *Communications Biology*. DOI: [10.1038/s42003-021-02782-y](https://doi.org/10.1038/s42003-021-02782-y)
    * **Böndel, Katharina B. and Schmid, Karl J.** (2021). Quinoa Diversity and Its Implications for Breeding. *NA*. DOI: [10.1007/978-3-030-65237-1_7](https://doi.org/10.1007/978-3-030-65237-1_7)
    * **Richards, Chris** (2021). Genomics of Plant Gene Banks: Prospects for Managing and Delivering Diversity in the Digital Age. *NA*. DOI: [10.1007/13836_2021_95](https://doi.org/10.1007/13836_2021_95)



