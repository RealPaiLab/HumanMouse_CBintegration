#!/bin/bash

# start jupyter on hpc
conda activate scenicplus
module load java
python -m ipykernel install --user --name=scenicplus --display-name "Python (scenicplus)"
jupyter-lab --no-browser --port=8888 --ip=0.0.0.0
