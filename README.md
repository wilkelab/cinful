# cinful



A fully automated pipeline to identify microcins along with their associated immunity proteins and export machinery



# Installation


First, make sure to clone this repository:

```bash
git clone https://github.com/wilkelab/cinful.git
```
All software dependencies needed to run cinful are available through conda and are specified in `cinful_conda.yml`, the following helper script can be used to generate the cinful conda environment `scripts/build_conda_env.sh`, to run this script, you will need to have conda installed, as well as mamba (which helps speed up installation). To install mamba, use the following command:

```bash
conda install mamba -c conda-forge
```

To build the environment, run
```bash
bash scripts/build_conda_env.sh
```

Once setup is complete, you can activate the environment with
```bash
conda activate cinful
```

# How to use

cinful takes a directory containing genome assemblies as input. All assemblies in the directory must end in `.fna`, if they end in a different extension, cinful will ignore them.

Snakemake is the core workflow management used by cinful, the main snakefile is located under `cinful/Snakefile`, which issues subroutines located in `cinful/rules`.

If installed properly, running `python cinful.py -h` will produce the following output.

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

# Example usage

There is a test dataset with an _E. coli_ genome assembly to test cinful on under `test/colcinV_Ecoli`, you can run cinful on this dataset by running the following from the initial cinful directory:

```bash
python cinful/cinful.py -d test/colcinV_Ecoli -o <output_directory> -t <#_of_threads>
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

# Contributing

cinful currently exists as a wrapper to a series of snakemake subroutines, so adding functionality to it is as simple as adding additional subroutines. If there are any subroutines that you see are needed, feel free to raise an issue, and I will be glad to guide you through the process of making a pull request to add that feature.

Additionally, since cinful primarily works through snakemake, it can also be used by simply running the snakefiles separately, so if additional configuration is needed, in terms of the types of input files, this can probably be achieved that way.
