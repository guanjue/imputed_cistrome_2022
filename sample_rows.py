#!/usr/bin/env python3
import gc
import os
import pickle

import fire
import h5py
import matplotlib.pyplot as plt
import seaborn as sns
from hyperopt.fmin import generate_trials_to_calculate
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import precision_recall_curve
from numpy import linalg as LA
import sklearn.metrics as metrics
import json
import lightgbm as lgb
import pandas as pd
import numpy as np



def sample_rows(inputfile, outputfile, seed, sample_nrows):
	d = np.array(pd.read_csv(inputfile, sep="\t", header=0))
	np.random.seed(seed)
	used_id_plot = np.random.choice(range(0,d.shape[0]), sample_nrows)
	ds = d[used_id_plot,:]
	np.savetxt(outputfile, ds, delimiter="\t", fmt="%s")



if __name__ == '__main__':
	fire.Fire(sample_rows)



