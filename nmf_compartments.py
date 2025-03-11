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

adata_file = f"{run_name}/sp.h5ad"
adata_vis = sc.read_h5ad(adata_file)
mod = cell2location.models.Cell2location.load(f"{run_name}", adata_vis)

##########################################################################
adata_vis1 = adata_vis[adata_vis.obs['sample'].isin(["334-I4","345-I3","334-I2","345-I5"])]

from cell2location import run_colocation
res_dict, adata_vis1 = run_colocation(
    adata_vis1,
    model_name='CoLocatedGroupsSklearnNMF',
    train_args={
      'n_fact': np.arange(5,30), # IMPORTANT: use a wider range of the number of factors (5-30)
      'sample_name_col': 'sample', # columns in adata_vis.obs that identifies sample
      'n_restarts': 3 # number of training restarts
    },
    export_args={'path': f'{run_name}/CoLocatedComb_d28/'}
)

##########################################################################
adata_vis2 = adata_vis[adata_vis.obs['sample'].isin(["345-N4","334-N3","345-N1"])]

from cell2location import run_colocation
res_dict, adata_vis2 = run_colocation(
    adata_vis2,
    model_name='CoLocatedGroupsSklearnNMF',
    train_args={
      'n_fact': np.arange(5, 30), # IMPORTANT: use a wider range of the number of factors (5-30)
      'sample_name_col': 'sample', # columns in adata_vis.obs that identifies sample
      'n_restarts': 3 # number of training restarts
    },
    export_args={'path': f'{run_name}/CoLocatedComb_naive/'}
)
