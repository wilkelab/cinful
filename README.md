# cinful

A fully automated pipeline to identify microcins along with their associated immunity proteins and export machinery

# Installation

There are two methods for installation, one uses pip and should be more user friendly. A backup method is to clone this repository and manually set up a conda environment.

## Recommended PyPI (pip) Installation

First you have to install anaconda (I would recommend miniconda), which can be found here: https://docs.conda.io/en/latest/miniconda.html

Once miniconda is installed, run the following commands in order. 

NOTE: Follow the instructions at each step and allow installations to complete before moving on to the next step. Do not paste all of the code at once into your shell.

```bash
conda create --name <your-env-name> python=3.8.13 pip

conda activate <your-env-name>

conda install mamba -c conda-forge

pip install cinful

cinful_init
```

Once installed, cinful can be called with `cinful` and can be used globally. I am working on a test to verify installation. As a workaround, you are able to download a test genome that contains microcin, MFP, PCAT, and immunity protein from https://github.com/wilkelab/cinful/blob/main/test/.

Once you've downloaded this test file, you can run cinful using instructions below.

## How to use

cinful takes a directory containing genome assemblies as input. All assemblies in the directory must end in `.fna`, if they end in a different extension, cinful will ignore them. 

Nested directories will explored recursively and all `.fna` files analyzed by cinful. Nested directories can be a good way to explore output, as the directory tree will be stored in as `cinful_id` in the output files.

Snakemake is the core workflow management used by cinful, the main snakefile is located under `cinful/Snakefile`, which issues subroutines located in `cinful/rules`.

If installed properly, running `cinful -h` will produce the following output:

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

Three output directories will be generated in your `assembly_directory` under a directory called `cinful_out`.
* `00_dbs`
  * This is the initial location of the databases of verified microcins, CvaB, and immunity proteins.
* `01_orf_homology`
  * Prodigal will generate Open Reading Frame (ORF) predictions for the input assemblies
  * Those ORFs will be searched against the previously mentioned databases
* `02_homology_results`
  * The results from all the homology searches will be merged here
* `03_best_hits`
  * The top hits from the homology results will be placed here

## Old Installation Method (not recommended)

Clone this repository:

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

There is a test dataset with an _E. coli_ genome assembly to test cinful on under `test/colcinV_Ecoli`, you can run cinful on this dataset by running the following from the initial cinful directory:

```bash
python cinful.py -d <genomes_directory> -o <output_directory> -t <threads>
```

# Contributing

cinful currently exists as a wrapper to a series of snakemake subroutines, so adding functionality to it is as simple as adding additional subroutines. If there are any subroutines that you see are needed, feel free to raise an issue, and I will be glad to guide you through the process of making a pull request to add that feature.

Additionally, since cinful primarily works through snakemake, it can also be used by simply running the snakefiles separately, so if additional configuration is needed, in terms of the types of input files, this can probably be achieved that way.
