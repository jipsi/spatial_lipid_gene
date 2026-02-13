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
from cell2location import run_colocation


def run_nmf_analysis(run_name, sample_list):
    """Run NMF co-location analysis on cell2location results.

    Loads the spatial mapping output and runs co-location analysis
    separately for infected and naive sample groups.

    Parameters
    ----------
    run_name : str
        Path to cell2location mapping outputs (containing sp.h5ad).
    sample_list : str
        Path to the sample list CSV (unused here but kept for consistency).

    Returns
    -------
    adata_vis1 : anndata.AnnData
        AnnData for infected samples with NMF results.
    adata_vis2 : anndata.AnnData
        AnnData for naive samples with NMF results.
    """
    # Library IDs
    libraries = [
        "334-I4", "345-N4", "345-I3", "334-N2-I4",
        "334-I2", "345-I5", "334-N3", "345-N1",
    ]

    # Load cell2location output
    adata_file = f"{run_name}/sp.h5ad"
    adata_vis = sc.read_h5ad(adata_file)
    mod = cell2location.models.Cell2location.load(f"{run_name}", adata_vis)

    # ── NMF co-location: infected samples (day 28) ──
    infected_samples = ["334-I4", "345-I3", "334-I2", "345-I5"]
    adata_vis1 = adata_vis[adata_vis.obs['sample'].isin(infected_samples)]

    res_dict, adata_vis1 = run_colocation(
        adata_vis1,
        model_name='CoLocatedGroupsSklearnNMF',
        train_args={
            'n_fact': np.arange(5, 30),
            'sample_name_col': 'sample',
            'n_restarts': 3,
        },
        export_args={'path': f'{run_name}/CoLocatedComb_d28/'},
    )

    # ── NMF co-location: naive samples ──
    naive_samples = ["345-N4", "334-N3", "345-N1"]
    adata_vis2 = adata_vis[adata_vis.obs['sample'].isin(naive_samples)]

    res_dict, adata_vis2 = run_colocation(
        adata_vis2,
        model_name='CoLocatedGroupsSklearnNMF',
        train_args={
            'n_fact': np.arange(5, 30),
            'sample_name_col': 'sample',
            'n_restarts': 3,
        },
        export_args={'path': f'{run_name}/CoLocatedComb_naive/'},
    )

    return adata_vis1, adata_vis2


if __name__ == "__main__":
    # ──────────────────────────────────────────────────────────
    # USER-CONFIGURABLE PATHS
    # Update these paths to match your local environment.
    # ──────────────────────────────────────────────────────────
    results_folder = '.'
    run_name = f'{results_folder}/cell2location_map'
    sample_list = './data/sd22.5.1.csv'

    run_nmf_analysis(run_name, sample_list)
