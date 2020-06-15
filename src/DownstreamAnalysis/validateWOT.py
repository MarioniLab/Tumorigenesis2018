# This avoids using X 
import matplotlib as mpl
mpl.use('Agg')

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import wot

import argparse

# Read in Arguments
parser = argparse.ArgumentParser()
parser.add_argument('--pca',  
                    help='Path to the batch corrected PCA')
parser.add_argument('--meta', 
                    help='Path to meta data')
parser.add_argument('--growth', 
                    help='learned growth path g.txt')
parser.add_argument('--time', 
                    help='name of time',default="time")
parser.add_argument('--out', 
                    help='Path to figure')

args = parser.parse_args()
PCA_PATH = args.pca
OBS_PATH = args.meta
LEARNED_GROWH_PATH = args.growth 
TIME = args.time

#Load Data
adata = wot.io.read_dataset(path= PCA_PATH, obs=[OBS_PATH, LEARNED_GROWH_PATH], var=None)
adata.shape

otmodel = wot.ot.OTModel(adata, growth_rate_field='g2', day_field=TIME, growth_iters=1, local_pca=0)

all_triplets_summary = wot.ot.compute_validation_summary(otmodel)

# save results
all_triplets_summary.to_csv('validation_summary.txt')
all_triplets_stats = all_triplets_summary.groupby(['interval_mid', 'name'])['distance'].agg([np.mean, np.std])
all_triplets_stats.to_csv('validation_summary_stats.txt')

wot.graphics.plot_ot_validation_summary_stats(all_triplets_stats)
plt.savefig(args.out)
