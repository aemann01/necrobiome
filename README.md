# necrobiome

This repository describes the analysis performed in the paper *Tarone et al.* ///////.

## Setup

This repository assumes you are running in a Unix environment (e.g., Mac OSX or Linux) and you have conda installed. 

To get this repository:

- Install and set up up anaconda or miniconda as described at the [bioconda
  documentation](https://bioconda.github.io/user/install.html), including
  setting up channels.
- [You should also have QIIME2 installed as a conda environment.](https://docs.qiime2.org/2020.8/install/)
- Clone this repository to your machine and change into the directory with

```bash
git clone https://github.com/aemann01/necrobiome.git && cd necrobiome/
```

- Run the following command to install the one of the environments

```bash
conda env create -f environment.yml

```

- To load a given environment run

```bash
conda activate necrobiome
```

- To turn off the environment run

```bash
conda deactivate
```
