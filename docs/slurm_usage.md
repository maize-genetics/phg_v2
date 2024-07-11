# PHGv2 - SLURM Usage Guidelines for `align-assemblies`

When running `align-assemblies` on a single machine, 
the `--total-threads` and `--in-parallel` parameters are used to 
determine how many threads to allocate for the alignment step and 
how many genomes to run in parallel. Users only need to call the 
`align-assemblies` command once and there is no pre-processing steps 
necessary.

However, for [HPC](https://en.wikipedia.org/wiki/High-performance_computing)
systems leveraging the [SLURM](https://en.wikipedia.org/wiki/Slurm_Workload_Manager)
job scheduler, parallel processes operate differently which will
require slight modifications to the alignment step. In this document, 
we will discuss how to use the `align-assemblies` command for usage 
with the SLURM work scheduler system.

## Workflow Overview

### 1 - Reference preparation
When running jobs in a SLURM array, the first step is to run the 
`align-assemblies` command with the `--just-ref-prep` parameter set 
to `true`. This will run the preliminary steps of aligning the 
reference genome to the GFF CDS, creating the `reference-sam` and 
`reference-cds-fasta` files needed for aligning the individual 
assemblies. The output of this step will be written to the user 
supplied output directory.

### 2 - SLURM data array file creation
The second step would be to create the SLURM data array file. This is 
a text file which has an `align-assemblies` command, one per line, 
for each assembly that will be run. It will use the `reference-sam` 
and `reference-cds-fasta` files created in the first step as input to 
the `align-assemblies` command. Users can create the file manually if 
they prefer, but for long lists of assemblies, manual curation can be 
tedious and prone to errors. PHGv2 provides a convenience command
(`prepare-slurm-align-file`) to automatically generate these files. See 
the ["Prepare SLURM File for Alignments"](#prepare-slurm-file-for-alignments)
section for further details.

### 3 - Job submission
The third step is to submit the SLURM job. When running a SLURM job, 
the `--in-parallel` parameter is not used as each `align-assemblies` 
command is assigned to separate computer nodes in the HPC system. 
Conversely, the number of threads to use is still specified by the 
`--total-threads` parameter.


## Prepare SLURM File for Alignments
The `prepare-slurm-align-file` command is a convenience method for 
creating a file that can be submitted to SLURM as a 
[data array job](https://slurm.schedmd.com/job_array.html):

```shell
phg prepare-slurm-align-file  \
    --phg-location /path/to/phg \
    --gff data/anchors.gff \
    --reference-file output/updated_assemblies/Ref.fa \
    --reference-sam output/alignment_files/Ref.sam \
    --reference-cds-fasta output/alignment_files/ref.cds.fasta \
    --asemblies data/assemblies_list.txt \
    --ref-max-align-cov 1 \
    --query-max-align-cov 1 \
    --total-threads 20 \
    --conda-env-prefix /path/to/conda/env \
    --output-dir /path/for/align-assemblies/output \
    --slurm-file output/slurm_align_file.txt \
    -o output/alignment_files
```


### Parameters
This command uses several parameters:
* `--phg-location` - The location of the phg executable.  The full path should be provided. This is needed to run the align-assemblies command.
  If it is not specified, the current directory, ie `./phg`, will be assumed.

* `--gff` - GFF file for the reference genome. This is used to
  identify full-length coding sequences to use as anchors

* `--reference-file` - The reference genome in
  [FASTA](https://en.wikipedia.org/wiki/FASTA_format) format.
    + > ℹ️ **Note**  
      The path to the reference genome should be the **updated version**
      that was created during the `prepare-assemblies` command.

* `--reference-sam` - Optional parameter. If this is specified, the 
  optional parameter `--reference-cds-fasta` must also be supplied. 
  When both are supplied, the software skips the creation of these 
  files and uses those supplied by the user. This is desirable when 
  the user is running multiple assembly alignments from a SLURM 
  data-array option and does not wish to realign the reference 
  multiple times. If specified, but `--reference-cds-fasta` is not, 
  **the software will throw an exception**.

* `--reference-cds-fasta` - Optional parameter. If this is specified, 
  the optional parameter `--reference-sam` must also be supplied. 
  When both are supplied, the software skips the creation of these 
  files and uses those supplied by the user. This is desirable when 
  the user is running multiple assembly alignments from a SLURM 
  data-array option and does not wish to realign the reference 
  multiple times. If specified, but `reference-sam` is not, the 
  software will throw an exception.

* `--assemblies` - A text file containing a list of **annotated**
  assembly genomes (_see the 
  ["Prepare Assembly FASTA files"](build_and_load.md#prepare-assembly-fasta-files) 
  section for further details_). The contents of the assembly list 
  file should be either full or relative paths to each uncompressed 
  assembly you would like to align. For example, since I am following
  the steps laid out in the 
  ["Build and Load Documentation"](build_and_load.md#align-assemblies),
  I can create a text file called `assemblies_list.txt` (placed in
  the `data/` subdirectory) and populate it with the following lines:

  ```
  output/updated_assemblies/LineA.fa
  output/updated_assemblies/LineB.fa
  ```
  
  Here, I am planning on aligning two genomes called `LineA` and
  `LineB`. Since these are created with the `prepare-assemblies` command
  and the output is located in a subdirectory called
  `output/updated_assemblies/` relative to my working directory, I 
  will also add that to the path.

  + > ⚠️ **Warning**  
    This text list **should not** contain the path to the reference
    genome since this is recognized in the `--reference-file` flag.

* `--ref-max-align-cov` - The maximum reference genome alignment 
  coverage. This is used in the `proali` 
  [command](build_and_load.md#internal-anchorwave-and-minimap2-commands).
  The default value is `1`.

* `--query-max-align-cov` - The maximum query genome alignment 
  coverage.  This is used in the `proali`
  [command](build_and_load.md#internal-anchorwave-and-minimap2-commands).
  The default value is `1`.

* `--total-threads` - How many threads would you like to allocate for
  the alignment step?

* `--conda-env-prefix` - Optional parameter that specifies the path
  to the Conda directory that contains the conda environment needed
  to run phg. If not set, conda env `phgv2-conda` in the default
  location will be used.

* `--slurm-file` - The name of the file that will be created that 
  contains the SLURM commands to run the `align-assemblies` command 
  for each assembly in the list of assemblies.

* `-o`, `--output-dir` - The name of the directory for the alignment outputs.


### Example output
Following along with the example data shown in the ["Build and Load"](build_and_load.md)
documentation, my example SLURM data-array file 
(`output/slurm_align_file.txt`) will look like the following since
I have two samples:

```
./phg align-assemblies --gff data/anchors.gff --output-dir output/alignment_files --reference-file output/update_assemblies/Ref.fa --reference-sam output/alignment_files/Ref.sam --reference-cds-fasta output/alignment_files/ref.cds.fasta --assembly-file data/test/smallseq/LineA.fa --total-threads 20 in-parallel 1  --ref-max-align-cov 1 --query-max-align-cov 1
./phg align-assemblies --gff data/anchors.gff --output-dir output/alignment_files --reference-file output/update_assemblies/Ref.fa --reference-sam output/alignment_files/Ref.sam --reference-cds-fasta output/alignment_files/ref.cds.fasta --assembly-file data/test/smallseq/LineB.fa --total-threads 20 in-parallel 1  --ref-max-align-cov 1 --query-max-align-cov 1
````

> [!NOTE]
> If the file specified by the `--assemblies` parameter contains _10_
> assemblies, the output file will contain _10_ lines, each with a 
> call to the `align-assemblies` command for a single assembly. If the 
> file specified by the `--assemblies` parameter contains _100_ 
> assemblies, the output file will contain _100_ lines.


### Integrating into SLURM jobs
Since the output from `prepare-slurm-align-file` is simply a list
of individual `align-assemblies` commands (each representing an
individual assembly), we must pass this along to an actual SLURM
array job. Below we have added an example SLURM script detailing
some example parameters and code setup you _may_ want to use for your
applications:

```shell
#!/bin/bash

#SBATCH --time=10:30:00  # walltime limit (HH:MM:SS)
#SBATCH --nodes=1  # number of nodes
#SBATCH --ntasks-per-node=40  # 40 processor core(s) per node X 2 threads per core
#SBATCH --mem=200G  # maximum memory per node
#SBATCH --partition=short  # standard node(s)
#SBATCH --job-name="10T_anchorwaveV2"
#SBATCH --mail-user=lcj34@cornell.edu  # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --array=0-4
# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
module load miniconda
module load java/17

eval "$(conda shell.bash hook)"

conda activate /project/buckler_lab_panand/lynn.johnson/phgv2-conda

echo "All jobs in this array have:"
echo "- SLURM array job id: ${SLURM_ARRAY_JOB_ID}"
echo "- SLURM array task count: ${SLURM_ARRAY_TASK_COUNT}"
echo "- SLURM array starting task: ${SLURM_ARRAY_TASK_MIN}"
echo "- SLURM array ending task: ${SLURM_ARRAY_TASK_MAX}"
echo "This job in the array has:"
echo "- SLURM job id: ${SLURM_JOB_ID}"
echo "- SLURM array task id: ${SLURM_ARRAY_TASK_ID}"

INPUTFILE=<your_align_command_input>

IFS=$'\n' read -d '' -r -a LINES < ${INPUTFILE}
LINE=${LINES[$SLURM_ARRAY_TASK_ID]}
eval ${LINE}
if [ $? -eq 0 ]
  then
    echo -e "$(date +"%D  %r")\tSuccess: ${LINE}"
    exit 0
  else
    echo -e "$(date +"%D  %r")\tFailed\t${LINE}"
    echo -e "$(date +"%D  %r")\tJobID\t${SLURM_JOB_ID}"
    echo -e "$(date +"%D  %r")\tTaskID\t${SLURM_ARRAY_TASK_ID}"
  exit 1
fi
```

...where the line that contains the `INPUTFILE` variable declaration
is the path to the `prepare-slurm-align-file` command output (in our
case, this would be `output/slurm_align_file.txt`). For more
information about how to get started with SLURM, please check out
the [official guides](https://slurm.schedmd.com/quickstart.html).
