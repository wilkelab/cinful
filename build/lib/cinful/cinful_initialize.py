import subprocess

def main():
    print('Running conda package installer for cinful..')
    subprocess.run(["mamba", "update", "python=3.8.13", "seqkit=0.15.0", "mafft=7.475", "hmmer=3.3.1", "blast=2.9.0", "diamond=2.0.11", "pandas=1.2.4", "numpy=1.19.2", "biopython=1.76", "snakemake=6.3.0", "prodigal=2.6.3", "pyhmmer=0.3.0", "--channel", "conda-forge", "--channel", "bioconda"])
    print('Installing pip packages..')
    subprocess.run(["pip", "install", "pyTMHMM==1.3.2"])
    subprocess.run(["pip", "install", "seqhash==1.0.0"])
    subprocess.run(["pip", "install", "blake3==0.2.0"])
    subprocess.run(["pip", "install", "cinful"])
    print('Environment setup complete for cinful')
    print('Please read log above to determine if setup succeeded')
    print('Run \'cinful -h\' to verify successful installation')

if __name__ == "__main__":
    main()
