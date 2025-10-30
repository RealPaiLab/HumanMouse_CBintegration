#!/bin/bash

# your code goes in software
# generated output files go in /data/username
# read input data from /isilon


docker run -d \
	-e PASSWORD=pailab \
	-p 8424:8787 \
	-v /.mounts/labs/pailab/private/icheong:/.mounts/labs/pailab/private/icheong \
	icheong/rstudio_scrnaseq:v4.2.2.20230129

	# -v /u/icheong/CBL_scRNAseq:/CBL_scRNAseq \
	# -v /data/icheong/scRNAseq/results:/CBL_scRNAseq/results \
	# --name icheong_rstudio_4.2.2 \