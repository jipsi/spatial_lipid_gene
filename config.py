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

# ──────────────────────────────────────────────────────────────
# USER-CONFIGURABLE PATHS
# Update these paths to match your local environment before running.
# ──────────────────────────────────────────────────────────────

results_folder = '.'

# Path to the directory containing Visium spatial data
sp_data_folder = '/mnt/lustre/users/sd530/visium/SD22.5.1/'

# Path to store / load reference regression model outputs
ref_run_name = '/users/sd530/scratch/python/cell2location/SD22.5.1/reference_signatures'

# Path to store cell2location spatial mapping outputs
run_name = f'{results_folder}/cell2location_map'

# Library IDs corresponding to Visium capture areas
libraries = [
    "334-I4", "345-N4", "345-I3", "334-N2-I4",
    "334-I2", "345-I5", "334-N3", "345-N1",
]

# CSV file listing sample names and directories
sample_list = './data/sd22.5.1.csv'


if __name__ == "__main__":
    from train_ref import train_reference_model
    from load_ref import load_reference_and_export
    from load_query import load_query_data
    from deconvolute import run_deconvolution

    # Step 1: Train the reference model
    adata_ref, inf_aver = train_reference_model(ref_run_name)

    # Step 2: Load reference and export signatures
    adata_ref, inf_aver = load_reference_and_export(ref_run_name)

    # Step 3: Load query (spatial) data
    adata_vis = load_query_data(
        sample_list, sp_data_folder, libraries
    )

    # Step 4: Run spatial deconvolution
    run_deconvolution(
        adata_vis, inf_aver, run_name, libraries
    )
