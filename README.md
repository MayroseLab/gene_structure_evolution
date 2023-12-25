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

```
