# Gene structure evolution
This repository contains code for research reproducibility purposes, related to the manuscript "The evolution of protein-coding gene structure in eukaryotes", currently in preparation.  
The code is not intended aa a standalone tool and it was not throughly tested. Feel free to use it at your own risk.
## To reproduce the analyses in the manuscript:
### 1. Environment setup
You'll need a Unix machine, with git available. You'll also need to setup conda or (preferebly) mamba for environment management. See [here](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html) for installation instructions and make sure you can create and activate environments from the command line before proceeding.
### 2. Clone the repository
E.g. `git clone git@github.com:MayroseLab/gene_structure_evolution.git`
### 3. Create the project environment
```
cd gene_structure_evolution
mamba env create -f env.yml
```
Once finished, activate with:
```
conda activate gene_structure_evolution
```
then run the following script to install additional R packages not available through conda:
```
Rscript install_R_packages.R
```
### 4. Option 1: Run the snakemake workflow
The snakemake workflow performs most of the analyses, starting from raw GFF3 files. This is rather heavy and will take some time and requires considerable computational resources. If you don't really need this, you can downolad the final results (see below) and skip this step. If you intend to run the snakemake workflow, it is recommended that you familiarize yourself with the [relevant documentation](https://snakemake.readthedocs.io/en/stable/tutorial/tutorial.html) first.
#### 4.1 Configure the compute system
Snakemake can run locally, on a compute cluster, or on the cloud. If you want to run on a cluster, you'll need to obtain an appropriate job submission script. Such a script that works for a PBS Pro cluster can be found under `snakemake/util/pbs_qsub_snakemake_wrapper.py`, but if you run on a different system, you might need to modify this script. See [here](https://snakemake.readthedocs.io/en/stable/tutorial/additional_features.html#cluster-execution) for more information on cluster execution and [here](https://snakemake.readthedocs.io/en/master/executing/cloud.html) for cloud execution.
#### 4.2 Configure the workflow
You'll need to configure the inputs and outputs of the pipeline. To reproduce the analysis described in the manuscript, just copy the two configuration files (yml and tsv) under `snakemake/config/`, and adjust the paths according to your file system.
#### 4.3 Run the workflow
The specific command you'll need for running the workflow depends on your working environment, but it will probably look something like this (assuming cluster execution):
```
conda activate gene_structure_evolution
snakefile="gene_structure_evolution/workflow/Snakefile"
qsub_script="gene_structure_evolution/util/pbs_qsub_snakemake_wrapper.py"
job_script="gene_structure_evolution/util/jobscript.sh"
snakemake -s $snakefile --configfile gene_structure_evolution/snakemake/config/conf.yml --cluster "python $qsub_script" --latency-wait 60 --use-conda -p -j 500 --jobscript "$job_script" --rerun-incomplete --keep-going --retries 1 >out 2>err
```
### 4. Option 2: Download the snakemake workflow result
If you don't want to run the workflow yourself, you can just download a pre-prepared result:
```
wget <DRYAD URL>
tar -zxvf RESULT.tar.gz
```
### 5. Run the jupyter notebooks
The notebooks contain post-analysis code which produces the final results, tables, and figures.  
Make sure you are in the right location and have the project environment activated:
```
cd gene_structure_evolution
conda activate gene_structure_evolution
```
Next, you need to start the jupyter server. The basic command is `jupyter lab`, howevr it may change depending on your computing environment. For example, if you are accessing the jupyter server on a remote machine, you may need to use some form of port-forwarding (see e.g. [here](https://ljvmiranda921.github.io/notebook/2018/01/31/running-a-jupyter-notebook/)).  
Once you have the server running, access it with a web browser. Navigate to `notebooks/`, and run first `notebook1.R.ipynb`, then `notebook2.py.ipynb`. In each of them, you'll need to set the value of the variable `snakemake_result_dir` to the path of your Snakemake result.
