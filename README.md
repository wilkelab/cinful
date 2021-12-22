# cinful



A fully automated pipeline to identify microcins along with their associated immunity proteins and export machinery



# Installation


First, make sure to clone this repository:

```bash
git clone https://github.com/tijeco/cinful.git
```
All software dependencies needed to run cinful are available through conda and are specified in `cinful_conda.yml`, the following helper script can be used to generate the cinful conda environment `scripts/build_conda_env.sh`, to run this script, you will need to have conda installed, as well as mamba (which helps speed up installation). To install mamba, use the following command:

```bash
conda install mamba -c conda-forge
```

Then simply run 
```bash
bash scripts/build_conda_env.sh
```

To set up the cinful environment, you can activate the environment with

```bash
conda activate cinful
```

# How to use

cinful takes a directory containing genome assemblies as input. All assemblies in the directory must end in `.fna`, if they end in a different extension, cinful will ignore them.

Snakemake is the core workflow management used by cinful, the main snakefile is located under `cinful/Snakefile`, which issues subroutines located in `cinful/rules`. To run cinful on your data set run the following command:

```bash
snakemake -d <assembly_directory> --threads <core_nums> --snakefile path/to/cinful/Snakefile
```

If installed properly, running `cinful -h` will produce the following output.

```
cinful

optional arguments:
  -h, --help            show this help message and exit
  -d DIRECTORY, --directory DIRECTORY
                        Must be a directory containing uncompressed FASTA formatted genome assemblies with
                        .fna extension. Files within nested directories are fine
  -o OUTDIR, --outDir OUTDIR
                        This directory will contain all output files. It will be nested under the input
                        directory.
  -t THREADS, --threads THREADS
                        This specifies how many threads to allow snakemake to have access to for
                        parallelization
```

## Workflow

The following workflow will be executed.
![cinful](figures/cinful_workflow.inkscape.svg)

Three output directories will be generated in your `assembly_directory` under a directory called `cinfulOut`.
* `00_dbs`
  * This is the initial location of the databases of verified microcins, CvaB, and immunity proteins.
* `01_orf_homology`
  * Prodigal will generate Open Reading Frame (ORF) predictions for the input assemblies
  * Those ORFs will be searched against the previously mentioned databases
* `02_homology_results`
  * The results from all the homology searches will be merged here
* `03_best_hits`
  * The top hits from the homology results will be placed here

# Example usage

There is a test dataset with an _E. coli_ genome assembly to test cinful on under `test/colcinV_Ecoli`, you can run cinful on this dataset by running the following:

```bash
snakemake -d test/colcinV_Ecoli --snakemake 
```
# Contributing

<!-- Currently, cinful is executed via directly issuing a snakemake command, what I will do in the future is create a python package that acts as a wrapper for snakemake to ease potential configuration of certain parameters within the workflow. -->

cinful currently exists as a wrapper to a series of snakemake subroutines, so adding functionality to it is as simple as adding additional subroutines. If there are any subroutines that you see are needed, feel free to raise an issue, and I will be glad to guide you through the process of making a pull request to add that feature.

Additionally, since cinful primarily works through snakemake, it can also be used by simply running the snakefiles separately, so if additional configuration is needed, in terms of the types of input files, this can probably be achieved that way.

<!-- Also, the pipeline currently runs end to end, though there may be cases where the user already has data for a certain part of the pipeline and would like to plug that in. Snakemake allows for that to be a possibility, so I will work to make a set of tutorials on how to do that through snakemake, and eventually the cinful python package will have options to do that as well. -->







<!-- ## Microcin homologs

* Filtering by length
   - Only need to search peptides that have < 150 AA
* Signal sequence
   - MSA of putative microcins can be used to evaluate the putative signal sequence based on what is known from the verifed dataset
 -->



<!-- ## Immunity protein homologs
* subcelluar localization and transmembrane helix will be predicted as a final filtering step
 -->

