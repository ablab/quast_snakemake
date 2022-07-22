# QUAST Snakemake pipeline

### Installation
You will need [Mamba](https://github.com/mamba-org/mamba)/[Mambaforge](https://github.com/conda-forge/miniforge#mambaforge) (**recommended**) or [Conda](https://conda.pydata.org/). 
If you already have Conda and want to switch to Mamba you can install it with
      
    $ conda install -n base -c conda-forge mamba

Once you have one of them, it is very easy to install [Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) (7.0.0 or higher) 
and other QUAST dependencies (specified in the [env YAML files](./envs/)) by typing  
     
    $ mamba create -c conda-forge -c bioconda -n quast_sm snakemake python=3.9  
    $ mamba env update -n quast_sm --file envs/basic.yaml  # all basic QUAST dependencies   
    $ mamba env update -n quast_sm --file envs/busco.yaml  # if conservative gene finding (BUSCO) is also needed
      
Note:
1. If you have Conda rather than Mamba, substitute `mamba` with `conda`. NB! We experienced a problem with installing BUSCO via `conda`.
2. Python 3.9 is needed since QUAST relies on pysam that is currently incompatible with Python 3.10 (as of 23.06.2022) 
3. You can use any other environment name instead of `quast_sm`

Now just activate the environment and you are ready to go!

    $ mamba activate quast_sm   
      
Note: you may need to run `mamba init` and restart the shell to enable `mamba activate` for the first time; or simply use `conda activate quast_sm`.

### Quick start:
To run QUAST on sample files from `./test_data/` type one of the following sample commands
    
    $ ./quast.py -o quast_test_output_with_ref test_data/contigs_1.fasta test_data/contigs_2.fasta -t 4 -r test_data/reference.fasta -g gene:test_data/genes.gff -g operon:test_data/operons.bed  # test with reference
    $ ./quast.py -o quast_test_output_no_ref test_data/contigs_1.fasta test_data/contigs_2.fasta