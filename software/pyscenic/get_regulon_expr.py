# Extract regulon expression values from a pySCENIC anndata object

import argparse
import os
import pandas as pd
import anndata as ad

parser = argparse.ArgumentParser()
parser.add_argument(
    # path to anndata file from pyscenic
    '--h5ad_file',
    required = True
)
parser.add_argument(
    # output directory
    '--out_dir',
    required = True
)
parser.add_argument(
    # output file name
    '--filename',
    required = True
)
args = parser.parse_args()

# load file
print('Loading anndata file from', args.h5ad_file)
adata = ad.read_h5ad(args.h5ad_file)

# save regulon expression (and other metadata) to CSV
fname = os.path.join(args.out_dir, args.filename)
print('Saving', fname)
adata.obs.to_csv(fname)
