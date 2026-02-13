import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import cell2location
from cell2location.models import RegressionModel


def train_reference_model(ref_run_name):
    """Train the cell2location reference regression model.

    Loads the single-cell reference h5ad, renames reserved columns,
    trains a RegressionModel, exports the posterior, and saves results.

    Parameters
    ----------
    ref_run_name : str
        Path to store / load reference model outputs.

    Returns
    -------
    adata_ref : anndata.AnnData
        Reference AnnData with posterior estimates.
    inf_aver : pd.DataFrame
        Inferred average expression per cell type.
    """
    # Load the single-cell reference data
    adata_file = f"{ref_run_name}/sc_source.h5ad"
    print(adata_file)
    adata_ref = sc.read_h5ad(adata_file)

    # When converting from Seurat objects, rename '_index' column (reserved in Python)
    print(adata_ref.__dict__['_raw'].__dict__['_var'])
    adata_ref.__dict__['_raw'].__dict__['_var'] = (
        adata_ref.__dict__['_raw'].__dict__['_var']
        .rename(columns={'_index': 'features'})
    )
    print(adata_ref.__dict__['_raw'].__dict__['_var'])

    adata_ref.var.set_index(adata_ref.var.features, inplace=True)
    adata_ref.var.index.name = 'index'

    # Prepare anndata for the regression model
    RegressionModel.setup_anndata(
        adata=adata_ref,
        batch_key='orig.ident',
        labels_key='CellMeshFine',
    )

    # Create and train the regression model
    mod = RegressionModel(adata_ref)
    mod.view_anndata_setup()
    mod.train(max_epochs=500, use_gpu=True)

    mod.plot_history(20)
    plt.savefig('ref.mod.history.pdf', bbox_inches='tight')

    # Export posterior estimates
    adata_ref = mod.export_posterior(
        adata_ref,
        sample_kwargs={'num_samples': 1000, 'batch_size': 2500, 'use_gpu': True},
    )

    # QC plot
    mod.plot_QC()
    plt.savefig('ref.mod.qc.pdf', bbox_inches='tight')

    # Save model and anndata
    mod.save(f"{ref_run_name}", overwrite=True)
    adata_file = f"{ref_run_name}/sc.h5ad"
    adata_ref.write(adata_file)

    # Export inferred average expression per cell type
    inf_aver = _export_inf_aver(adata_ref)

    return adata_ref, inf_aver


def _export_inf_aver(adata_ref):
    """Extract inferred average expression from the reference model.

    Parameters
    ----------
    adata_ref : anndata.AnnData
        Reference AnnData after model training.

    Returns
    -------
    inf_aver : pd.DataFrame
        Inferred average expression per cell type.
    """
    factor_names = adata_ref.uns['mod']['factor_names']
    columns = [f'means_per_cluster_mu_fg_{i}' for i in factor_names]

    if 'means_per_cluster_mu_fg' in adata_ref.varm.keys():
        inf_aver = adata_ref.varm['means_per_cluster_mu_fg'][columns].copy()
    else:
        inf_aver = adata_ref.var[columns].copy()

    inf_aver.columns = factor_names
    return inf_aver
