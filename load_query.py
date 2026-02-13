import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

import warnings
warnings.filterwarnings('ignore')

from csv import DictReader
from scipy.sparse import csr_matrix


def read_and_qc(sample_name, path):
    """Read a 10X Visium spatial experiment and calculate QC metrics.

    Parameters
    ----------
    sample_name : str
        Name of the sample.
    path : str
        Path to the Visium data directory.

    Returns
    -------
    adata : anndata.AnnData
        Annotated data matrix with QC metrics.
    """
    adata = sc.read_visium(
        path,
        count_file=f"{sample_name}.h5",
        library_id=sample_name,
        load_images=True,
    )
    adata.obs['sample'] = sample_name
    adata.var['SYMBOL'] = adata.var_names
    print(adata.var['SYMBOL'])

    # Calculate QC metrics
    adata.X = adata.X.toarray()
    sc.pp.calculate_qc_metrics(adata, inplace=True)
    adata.X = csr_matrix(adata.X)

    # Mitochondrial gene fraction (may be absent in Visium FFPE)
    adata.var['MT'] = [gene.startswith('mt-') for gene in adata.var['SYMBOL']]
    adata.obs['mt_frac'] = (
        adata[:, adata.var['MT'].tolist()].X.sum(1).A.squeeze()
        / adata.obs['total_counts']
    )
    print(adata.obs['mt_frac'])

    # Add sample name to obs names
    adata.obs["sample"] = [str(i) for i in adata.obs['sample']]
    adata.obs_names = adata.obs["sample"] + '_' + adata.obs_names
    adata.obs.index.name = 'spot_id'
    adata.var_names_make_unique()

    return adata


def select_slide(adata, s, s_col='sample'):
    """Select data for one slide from a multi-sample spatial AnnData.

    Parameters
    ----------
    adata : anndata.AnnData
        AnnData with multiple spatial experiments.
    s : str
        Name of the selected experiment.
    s_col : str
        Column in adata.obs listing experiment name for each location.

    Returns
    -------
    slide : anndata.AnnData
        Subset AnnData for the selected slide.
    """
    slide = adata[adata.obs[s_col].isin([s]), :]
    s_keys = list(slide.uns['spatial'].keys())
    s_spatial = np.array(s_keys)[[s in k for k in s_keys]][0]
    slide.uns['spatial'] = {s_spatial: slide.uns['spatial'][s_spatial]}
    return slide


def load_query_data(sample_list, sp_data_folder, libraries):
    """Load, QC, and concatenate all spatial query samples.

    Parameters
    ----------
    sample_list : str
        Path to a CSV file listing sample names and directories.
    sp_data_folder : str
        Base path to the Visium data directories.
    libraries : list of str
        Library IDs corresponding to Visium capture areas.

    Returns
    -------
    adata_vis : anndata.AnnData
        Concatenated AnnData for all spatial samples.
    """
    sample_data = pd.read_csv(sample_list)

    slides = []
    with open(sample_list, 'r') as read_obj:
        csv_dict_reader = DictReader(read_obj)
        for row in csv_dict_reader:
            print(row['sample_name'], row['sample_dir'])
            sample = row['sample_name']
            sample_dir = sp_data_folder + row['sample_dir']
            slides.append(read_and_qc(sample, path=sample_dir))

    # Combine anndata objects
    adata_vis = slides[0].concatenate(
        slides[1:],
        batch_key="sample",
        uns_merge="unique",
        batch_categories=sample_data['sample_name'],
        index_unique=None,
    )

    # Plot QC for each sample
    fig, axs = plt.subplots(len(slides), 4, figsize=(15, 4 * len(slides) - 4))
    for i, s in enumerate(adata_vis.obs['sample'].unique()):
        slide = select_slide(adata_vis, s)

        sns.histplot(slide.obs['total_counts'], kde=False, ax=axs[i, 0])
        axs[i, 0].set_xlim(0, adata_vis.obs['total_counts'].max())
        axs[i, 0].set_xlabel(f'total_counts | {s}')

        sns.histplot(
            slide.obs['total_counts'][slide.obs['total_counts'] < 20000],
            kde=False, bins=40, ax=axs[i, 1],
        )
        axs[i, 1].set_xlim(0, 20000)
        axs[i, 1].set_xlabel(f'total_counts | {s}')

        sns.histplot(
            slide.obs['n_genes_by_counts'], kde=False, bins=60, ax=axs[i, 2],
        )
        axs[i, 2].set_xlim(0, adata_vis.obs['n_genes_by_counts'].max())
        axs[i, 2].set_xlabel(f'n_genes_by_counts | {s}')

        sns.histplot(
            slide.obs['n_genes_by_counts'][slide.obs['n_genes_by_counts'] < 6000],
            kde=False, bins=60, ax=axs[i, 3],
        )
        axs[i, 3].set_xlim(0, 6000)
        axs[i, 3].set_xlabel(f'n_genes_by_counts | {s}')

    plt.savefig("qc_vis_samples.pdf", bbox_inches='tight')

    return adata_vis
