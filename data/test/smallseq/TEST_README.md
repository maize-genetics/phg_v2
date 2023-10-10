# Small test data

The files here were created by lcj34 for testing the phg_v2 pipeline.  The methods used are below:

## creating the fasta files
I used the original PHG's CreateSmallGenomesPlugin to create the Ref.fa, LineA.fa and LineB.fa files.
This plugin will create a single chromosome's worth of data using the parameters supplied.  I ran it twice
so I would have 2 chromosomes worth of data. I then concatenated the two chromosomes together to create single fastas
with multiple chromosomes.

The first time it was run, this was the command and parameters used:
```
/Users/lcj34/git/tassel-5-standalone/run_pipeline.pl -debug -CreateSmallGenomesPlugin -geneLength 1000 -interGeneLength 5000 -refInterGeneDup 0.2 -refInterGeneDelete 0.1 -refGeneInsert 0.2 -numberOfGenes 10 -baseDir /Users/lcj34/notes_files/phg_v2/smallSeq_chrom1 -endPlugin
```

The fastas created were stored to a smallSeq_chrom1 directory.
I then ran the plugin again, this time changing the outputDir to smallSeq_chrom2, and the parameters to those shown below:
```
/Users/lcj34/git/tassel-5-standalone/run_pipeline.pl -debug -CreateSmallGenomesPlugin -geneLength 1000 -interGeneLength 5000 -refInterGeneDup 0.2 -refInterGeneDelete 0.1 -refGeneInsert 0.3 -numberOfGenes 10 -baseDir /Users/lcj34/notes_files/phg_v2/smallSeq_chrom2 -endPlugin
```

I edited the chrom2 version to have an id line with ">2".
In addition, I edited the chrom2 anchors.bed and anchors.gff files to have the correct chromosome number (2 instead of 1).

## final gff changes
In addition, in the anchors.gff file in the smallSeq_chrom2 folder, I changed the gene0 ... gene9 to be gene10..gene19.

## concat the files
I then concatenated the related smallseq_chrom1 and smallseq_chrom2 files together to create the Ref.fa, LineA.fa and LineB.fa files
that appear in a folder named smallseq_final.

## final gff changes
After the above steps, one more editing of the anchors.gff file.

## Additional assemblies
Though not yet stored here (update this README.md when they are!) we should have a LineC.fa
created that has a selection of genes removed.  This will be used to test the ability of
anchorwave to align when there is substantial missing data, and it will test our MAFToGVCF code
to verify it correctly add strings of N's to the output VCF file.

