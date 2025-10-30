# MB_scRNAseq

## Code setup (`docker` directory)

**!!DOCKER WAS NOT USED AFTER JAN 2024; EVERYTHING WAS MOVED TO THE OICR HPC CLUSTER AND PACKAGES ARE MANAGED WITH CONDA ENVIRONMENTS!!**

The code env runs Rstudio server via Docker. These instructions are intended to be run in a terminal emulator such as iTerm2.

On the lab server run:

```
git clone git@github.com:RealPaiLab/MB_scRNAseq.git
cd MB_scRNAseq/
./startDocker.sh
```

This runs the Docker image on port 8877 on the lab server.

On your laptop run:

```
git clone git@github.com:RealPaiLab/MB_scRNAseq.git
cd MB_scRNAseq/
./mapSSH_onlaptop.sh
```

This creates an SSH tunnel from port 8877 on your laptop to port 8877 on the lab server.

You should now be able to run Rstudio on your web browser at: `localhost:8877`. Username is `rstudio`, Password is `pailab`.
Your software directory is mapped to `/software` and your data directory is at `/data`.

## `conda` directory

YAML files for the different conda environments are stored here. See the directory's README for more details about the installation of certain packages not available through conda.

## `figures` directory

This contains high quality figures (e.g. for presentations) and the code that was used to generate them.

## `HumanMouseUBC` directory

Scritps and figures associated with the human/mouse UBC paper. Output can be found in the `project/HumanMouseUBC` directory.

## `lab_notebook` directory
    
A 'virtual' lab notebook of code that was written and anlayses that were run.

## `MouseAdol-scRNAseq` directory

This contains code for the supplementary single-cell analysis for Shraddha's autism/GWAS paper.

## `software` directory

This contains all the source code that was used to analyze data and generate results.

The subdirectories are divided into code for human, mouse, integrated, and other general utilities/functions.

## `thesis` directory

This contains code to generate all tables and figures in the thesis.
