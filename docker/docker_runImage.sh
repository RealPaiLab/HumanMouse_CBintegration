#!/bin/bash

# start dnmt docker with mapped port for rstudio

docker run -d --rm \
       	--name spai_epi \
       	-p 8880:8787 \
       	-v /u/spai/software:/home/rstudio/software \
       	-v /u/spai/data:/home/rstudio/data \
	-v /.mounts/labs/pailab:/home/rstudio/isilon \
	bioc319_dnam_analysis:v4

docker exec -it spai_epi /bin/bash
