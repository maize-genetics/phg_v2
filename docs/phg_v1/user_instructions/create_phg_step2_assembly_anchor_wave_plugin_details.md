!!! warning "Legacy Documentation - PHG Version 1"

    This section contains documentation for **PHG version 1**, which is
    no longer actively developed. It is preserved here for archival and
    historical reference only. If you are looking to use the Practical
    Haplotype Graph, please refer to the [PHG v2 documentation](../../index.md),
    which reflects the current version of the software.

# AssemblyMAFFromAnchorWavePlugin Details

The first step in loading assembly haplotypes to the database is to align the assemblies to the reference genome with the AssemblyMAFFromAnchorWavePlugin plugin.
It can be run separately using a shell script similar to the one detailed below:

```
WORKING_DIR=/workdir/lcj34/anchorwaveTesting/
DOCKER_CONFIG_FILE=/phg/config.txt

# correct parameters for this plugin should be in the config.txt file
# replace <version tag> with the specific PHG version you wish to run
docker1 run --name anchorwave_assemblies_maf --rm \
    -v ${WORKING_DIR}/:/phg/ \
    -t maizegenetics/phg:<version tag> \
    /tassel-5-standalone/run_pipeline.pl -Xmx100G -debug -configParameters ${DOCKER_CONFIG_FILE} \
    -AssemblyMAFFromAnchorWavePlugin -endPlugin
```

The config file must contain the entries below to run this step. The AssemblyMAFFromAnchorWavePlugin does not require database access.  Replace the entry values below with your own parameter values.  The entries listed here assume you have setup your directory structure using the MakeDefaultDirectoryPlugin.  It also assumes you are running via docker with the docker mount points as defined in the example script above.  Note the parameter values have paths that are docker specific.


```

AssemblyMAFFromAnchorWavePlugin.outputDir=/phg/outputDir
AssemblyMAFFromAnchorWavePlugin.keyFile=/phg/Ia453_K0326Y_keyfile.txt
AssemblyMAFFromAnchorWavePlugin.gffFile=/phg/inputDir/reference/Zm-B73-REFERENCE-NAM-5.0_Zm00001e.1.gff3
AssemblyMAFFromAnchorWavePlugin.refFasta=/phg/inputDir/reference/Zm-B73-REFERENCE-NAM-5.0.fa
AssemblyMAFFromAnchorWavePlugin.threadsPerRun=4
AssemblyMAFFromAnchorWavePlugin.numRuns=2
```


# Details

This class uses anchorwave to align full genome assembly fastas to a reference.
It expects a keyfile that contains information on the assembly fastas, including information
on the local and remote location of the fastas.  The former is used for aligning, the
latter (the remote location) is used in later steps to populate the genome_file_data table
with information on a public site where the reference and assembly fastas may be found.

For more information on the anchorwave alignment method, please check out [this link](https://github.com/baoxingsong/AnchorWave).

## Parameters

Run the following command to print a current list of parameters for this plugin:

```
docker run --rm  maizegenetics/phg /tassel-5-standalone/run_pipeline.pl -AssemblyMAFFromAnchorWavePlugin
```

## Parameter Descriptions

* outputDir: (required) Directory to which all output files (*.maf and others) will be written
* keyFile: (required) Name of the keyfile to process.  eetails on the keyfile are in a separate header below
* gffFile: (required) Full path to a  GFF file associated with the reference genome
* refFasta: (required) Full path to the reference fasta file to be used when aligning the assembly genomes
* threadsPerRun: (optional: default=2) The number of threads to use for each processed assembly.  This value plus the value for the numRUns parameter shoudl be determined based on available system threads and memory.
* numRuns: (optional: default=1) Number of simultaneous assemblies to process. The anchorwave application can take up to 50G per thread for each assembly processed, plus some overhead. Consider this memory factor when providing values for threadsPerRun and numRun
* minimap2Location: (optional: default="minimap2") Location of Minimap2 on file system.  This can be defaulted if it is on the PATH environment variable.  Anchorwave uses minimap2 for parts of the alignment.
* anchorwaveLocation: (optional: default="anchorwave") Location of anchorwave on file system.  This can be defaulted if it is on the PATH environment variable.
* refMaxAlignCov: (optional: default=1) anchorwave proali parameter R, indicating reference genome maximum alignment coverage.
* queryMaxAlignCov: (optional: default=1) anchorwave proali parameter Q, indicating query genome maximum alignment coverage

## Keyfile
The AssemblyMAFFromAnchorWavePlugin keyfile must have columns AssemblyServerDir, AssemblyGenomeFasta, AssemblyDir, AssemblyFasta, and AssemblyDBName.  The AssemblyFasta column should contain the name of the assembly fasta file for aligning.  This is a full genome fasta.  The AssemblyGenomeFasta column should contain the name of the full genome fasta from which the assembly fasta came (it may be the same name as the AssemblyGenomeFasta).

The keyfile may also have an optional "Description" column.  If present, the data from this column will be used by PopulatePHGDBPipelinePlugin when creating the keyfile for the LoadHaplotypesFromGVCFPlugin().  It becomes the description for the sample and should contain any assembly genome data the user would like stored in the genotypes table "description" field associated with the assembly taxon.

The keyfile should be of tab-deliminted format.  The headers should look as below:
`AssemblyServerDir       AssemblyGenomeFasta     AssemblyDir     AssemblyFasta   AssemblyDBName  Description`

An example keyfile for this plugin is below. Note all fields are tab-delimited and because this is run in a docker, the reference and assembly directories are docker paths.  The AssemblyServerDir should be the more permanent accessible location outside the docker where these assembly fastas are stored.


```
AssemblyServerDir       AssemblyGenomeFasta     AssemblyDir     AssemblyFasta   AssemblyDBName  Description
irods:AssemblyFilePath  LineA.fa        /phg/inputDir/assemblies/       LineA.fa        LineA_Assembly  smallSeq description for assembly LineA
irods:AssemblyFilePath  LineB.fa        /phg/inputDir/assemblies/       LineB.fa        LineB_Assembly  smallSeq description for assembly LineB
```

## Memory considerations
Running anchorwave can be very memory intensive. It can take up to 50G per thread when running anchorwave.  Consider both the number of threads available on your machine as well as the memory that will be used for each.
For example: If your machine has 512G, no more than 10 threads may be used when running anchorwave for a single run.  If you want to run 2 alignments in parallel and your machine has 512G, no more than 4 threadsPerRun.

```
memoryCost ~ <number of assemblies to process> * (80+(<numThreads-1>)*50)

e.g.: 2 assemblies with 4 threads each:
  memoryCost ~ 2 * (80+(4-1)*50) == 460GB
```

In the calculation above, the first 80G is for the thread that was subtracted in <numThreads-1>)  This
includes 50Gb for the alignment matrix, another 30Gb for genomes and longest path result etc.

When a user has 10 threads available, there are different options, e.g. 2 runs with 5 threads each,
or 5 runs with 2 threads each.  AnchorWave uses a thread for each collinear block, for some cases the last
single collinear block might take sometime during which it will be using a single thread.
So, setting of 5 runs with 2 threads each is faster than 2 runs with 5 threads each, but it incurs a higher memory cost.