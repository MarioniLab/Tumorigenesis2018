import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import argparse

import wot

# Read in Arguments
parser = argparse.ArgumentParser()
parser.add_argument('--pca',  
                    help='Path to the batch corrected PCA')
parser.add_argument('--meta', 
                    help='Path to meta data')
parser.add_argument('--out', 
                    help='Path to out folder')
parser.add_argument('--time', 
                    help='Path to out folder',default='time')

#Parser
args = parser.parse_args()
PCA_PATH = args.pca
OBS_PATH = args.meta
OUT_DIR = args.out
TIME=args.time

#Load Data
adata = wot.io.read_dataset(path= PCA_PATH, obs=OBS_PATH, var=None)
adata.shape

otmodel = wot.ot.OTModel(adata, epsilon = 0.05, lambda1 = 1, lambda2 = 50, day_field=TIME, growth_iters=3, local_pca=0)

otmodel.compute_all_transport_maps(tmap_out=OUT_DIR)
