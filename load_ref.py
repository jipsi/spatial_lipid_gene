import scanpy as sc
import pandas as pd
import numpy as np

import cell2location


def load_reference_and_export(ref_run_name):
    """Reload a trained reference model and export cell-type signatures.

    Parameters
    ----------
    ref_run_name : str
        Path where the reference model was saved.

    Returns
    -------
    adata_ref : anndata.AnnData
        Reference AnnData with model results.
    inf_aver : pd.DataFrame
        Inferred average expression per cell type.
    """
    adata_file = f"{ref_run_name}/sc.h5ad"
    adata_ref = sc.read_h5ad(adata_file)

    adata_ref.var.set_index(adata_ref.var.features, inplace=True)
    adata_ref.var.index.name = 'index'

    mod = cell2location.models.RegressionModel.load(f"{ref_run_name}", adata_ref)

    # Export inferred average expression per cell type
    factor_names = adata_ref.uns['mod']['factor_names']
    columns = [f'means_per_cluster_mu_fg_{i}' for i in factor_names]

    if 'means_per_cluster_mu_fg' in adata_ref.varm.keys():
        inf_aver = adata_ref.varm['means_per_cluster_mu_fg'][columns].copy()
    else:
        inf_aver = adata_ref.var[columns].copy()

    inf_aver.columns = factor_names
    inf_aver.iloc[0:5, 0:5]

    return adata_ref, inf_aver
