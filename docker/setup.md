# Docker and RStudio for scRNA-seq

This document outlines how to set up a Docker container with RStudio Server for single-cell RNA-seq analysis.

## Starting the Docker container

To start the Docker container, run the `startDocker.sh` script. Change the command arguments as needed. Note the port that is used (in the `-p` option). For example, to publish RStudio's default port (8787) to port 1234 on the server, you use:
```
docker run -p 1234:8787 ...
```

See the official documentation for more information about the `docker run` command.

## Installing additional Linux packages

Enter the container as root by running `enterContainer.sh`. Make sure to change `icheong_scRNAseq` to the name of your container.

Once inside, use upgrade all packages by running:
```
apt update
apt upgrade
```

Then install the following packages:
```
apt install zlib1g-dev libxml2 libglpk40
```

Exit the Docker container shell.

## Accessing RStudio

SSH into the lab server again, but this time, set up port forwarding from the server port (1234 in the previous example) to a port on your local machine. Access RStudio by typing in "localhost:abcd" in your browser, where "abcd" is your local port you chose when you SSH'd into the server.

## Installing R packages

(See the Dockerfile for an up-to-date list of Linux, R, and Bioconductor that need to be installed.))

```
# install tidyverse
install.packages("tidyverse")

# install Bioconductor
install.packages("BiocManager")

# install Seurat
install.packages("Seurat")

# install Monocle 3 (see https://cole-trapnell-lab.github.io/monocle3/docs/installation/)
BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'Matrix.utils'))
install.packages("devtools")
devtools::install_github('cole-trapnell-lab/leidenbase')
devtools::install_github('cole-trapnell-lab/monocle3')
```

If you get an error saying "Installation paths not writeable...", enter the container as root, start R using the command line, and run the install again.

## Troubleshooting

### `Setting LC_* failed` warning on start-up

When starting RStudio, you may encounter warning messages similar to this:

```
During startup - Warning messages:
1: Setting LC_CTYPE failed, using "C"
2: Setting LC_COLLATE failed, using "C"
3: Setting LC_TIME failed, using "C"
(etc.)
```

I suspect that the warning results from the locale not being set properly (not sure why this happens though).

The warnings can be eliminated as follows:

1. Enter the Docker container as root (you can run the `enterContainer.sh` script).

1. Open `/etc/locale.gen` and uncomment the line "en_US.UTF-8 UTF-8". Save and close the file. (To edit the file, you may need to install an editor like `nano` or `vim`.)

1. Back as the container root, run 

    ```
    root:/# locale-gen en_US.UTF-8
    root:/# update-locale LANG=en_US.UTF-8
    ```

1. Quit and restart RStudio (red power button in the top right corner).

The warnings should no longer appear when RStudio starts. You can check the locale by running `Sys.getlocale()` or `sessionInfo()`.

For a more permanent solution, setting the environment variable in the Dockerfile may be the way to go (see https://askubuntu.com/questions/581458/how-to-configure-locales-to-unicode-in-a-docker-ubuntu-14-04-container). I haven't tried it yet, so no clue if it'll work or not.
