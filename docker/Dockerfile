# Create an RStudio Docker image for scRNA-seq analysis
FROM rocker/rstudio:4.2.2

# install Linux packages
RUN apt-get update && apt-get upgrade -y && apt-get install -y \
    libbz2-dev \
    libgdal-dev \
    libglpk40 \
    libproj-dev \
    libudunits2-0 \
    libxml2 \
    libxt6 \
    libxt-dev \
    zlib1g-dev \
    libcairo2-dev \
    patch

# install R packages, including Seurat and monocle3
RUN R -e 'install.packages(c( \
    "tidyverse", \
    "argparse", \
    "BiocManager", \
    "bookdown", \
    "devtools", \
    "ggalluvial", \
    "gprofiler2", \
    "ggpubr", \
    "ggVennDiagram", \
    "harmony", \
    "hdf5r", \
    "pals", \
    "R.utils", \
    "Seurat", \
    "styler", \
    "svglite", \
    "quarto"))'
RUN R -e 'BiocManager::install(c( \
    "BiocGenerics", \
    "biomaRt", \
    "BSgenome.Hsapiens.UCSC.hg38", \
    "ComplexHeatmap", \
    "DelayedArray", \
    "DelayedMatrixStats", \
    "limma", \
    "memes", \
    "S4Vectors", \
    "SingleCellExperiment", \
    "SummarizedExperiment", \
    "batchelor", \
    "HDF5Array", \
    "terra", \
    "ggrastr", \
    "Matrix.utils", \
    "slingshot", \
    "scran", \
    "SingleR", \
    "clusterProfiler", \
    "org.Hs.eg.db"))'
RUN R -e 'devtools::install_github("cole-trapnell-lab/leidenbase"); \
    devtools::install_github("cole-trapnell-lab/monocle3"); \
    remotes::install_github("satijalab/seurat-wrappers"); \
    remotes::install_github("mojaveazure/seurat-disk"); \
    devtools::install_github("SONGDONGYUAN1994/PseudotimeDE")'

# install Quarto for command line
RUN /bin/bash -c 'wget https://github.com/quarto-dev/quarto-cli/releases/download/v1.2.329/quarto-1.2.329-linux-amd64.deb; \
    dpkg -i quarto-1.2.329-linux-amd64.deb; \
    rm -f quarto-1.2.329-linux-amd64.deb'

# install MEME Suite dependencies
RUN /bin/bash -c 'cpan \
    File::Which \
    HTML::PullParser \
    HTML::Template \
    HTML::TreeBuilder \
    JSON \
    XML::Simple \
    XML::Parser::Expat \
    Sys::Info'

# install MEME Suite command line tools
WORKDIR /home/rstudio/
RUN /bin/bash -c 'wget https://meme-suite.org/meme/meme-software/5.5.1/meme-5.5.1.tar.gz && \
    tar xzf meme-5.5.1.tar.gz && \
    cd meme-5.5.1/ && \
    ./configure --prefix=/home/rstudio/meme --enable-build-libxml2 --enable-build-libxslt && \
    make && \
    make test && \
    make install && \
    cd ../ && \
    rm -rf meme-5.5.1.tar.gz meme-5.5.1/'

# custom fonts and settings for RStudio
WORKDIR /home/rstudio/.config/rstudio
COPY rstudio-prefs.json ./
COPY ["JetBrains Mono/", "./fonts/JetBrains Mono"]
WORKDIR /
