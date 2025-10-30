#!/bin/bash

# start dnmt docker with mapped port for rstudio

docker run -d --rm \
        --name spai_bioc319 \
        -p 8880:8787 \
        -v /u/spai/software:/home/rstudio/software \
        -v /u/spai/data:/home/rstudio/data \
        -v /.mounts/labs/pailab/private:/home/rstudio/isilon \
        -v /.mounts/labs/pailab/src:/home/rstudio/src \
        bioc319_dnam_analysis:v4

docker exec -it spai_bioc319 /bin/bash
