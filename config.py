import scanpy as sc
import anndata
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

import cell2location
import scvi

from matplotlib import rcParams
rcParams['pdf.fonttype'] = 42

results_folder = '.'
#relative path to where vis data lives
sp_data_folder = '/mnt/lustre/users/sd530/visium/SD22.5.1/'

# create paths and names to results folders for reference regression and cell2location models
ref_run_name = '/users/sd530/scratch/python/cell2location/SD22.5.1/reference_signatures'
run_name = f'{results_folder}/cell2location_map'

libraries = ["334-I4","345-N4","345-I3","334-N2-I4","334-I2","345-I5","334-N3","345-N1"]

sample_list = './data/sd22.5.1.csv'

exec(open("train_ref.py").read())