# Extract genes from inferred regulons in a `*.regulons.dat` file

import argparse
import os
import pickle
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument(
    # inferred regulons in binary format from pyscenic
    '--regulon_file',
    required=True
)
parser.add_argument(
    # output directory
    '--out_dir',
    required=True
)
args = parser.parse_args()

# load file
print('Loading regulons from', args.regulon_file)
with open(args.regulon_file, 'rb') as f:
    regulons = pickle.load(f)

# get genes
for regulon in regulons:
    genes = pd.DataFrame({'gene': regulon.genes, 'weight': regulon.weights})

    fname = os.path.join(args.out_dir, (regulon.transcription_factor + '.csv'))
    print('Saving', fname)
    genes.to_csv(fname, index=False)

