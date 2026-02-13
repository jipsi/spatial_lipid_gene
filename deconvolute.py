import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

import cell2location


def run_deconvolution(adata_vis, inf_aver, run_name, libraries):
    """Run cell2location spatial mapping (deconvolution).

    Parameters
    ----------
    adata_vis : anndata.AnnData
        Spatial AnnData object (query data).
    inf_aver : pd.DataFrame
        Inferred average expression per cell type from the reference model.
    run_name : str
        Path to store cell2location mapping outputs.
    libraries : list of str
        Library IDs for spatial plotting.

    Returns
    -------
    adata_vis : anndata.AnnData
        AnnData with cell2location results.
    """
    # Find shared genes and subset both anndata and reference signatures
    intersect = np.intersect1d(adata_vis.var_names, inf_aver.index)
    adata_vis = adata_vis[:, intersect].copy()
    inf_aver = inf_aver.loc[intersect, :].copy()

    # Prepare anndata for cell2location model
    cell2location.models.Cell2location.setup_anndata(
        adata=adata_vis, batch_key="sample"
    )

    # Create and train the model
    mod = cell2location.models.Cell2location(
        adata_vis,
        cell_state_df=inf_aver,
        N_cells_per_location=30,
        detection_alpha=20,
    )
    mod.view_anndata_setup()

    mod.train(
        max_epochs=30000,
        batch_size=None,
        train_size=1,
        use_gpu=True,
    )

    # Plot ELBO loss history (removing first 1000 epochs)
    mod.plot_history(1000)
    plt.legend(labels=['full data training'])
    plt.savefig("cell2location_full_training.pdf", bbox_inches='tight')

    # Export posterior estimates
    adata_vis = mod.export_posterior(
        adata_vis,
        sample_kwargs={
            'num_samples': 1000,
            'batch_size': mod.adata.n_obs,
            'use_gpu': True,
        },
    )

    # QC plot
    mod.plot_QC()
    plt.savefig("cell2location_full_training_qc.pdf", bbox_inches='tight')

    # Save model and anndata
    mod.save(f"{run_name}", overwrite=True)
    adata_file = f"{run_name}/sp.h5ad"
    adata_vis.write(adata_file)

    # Spatial QC across batches
    mod.plot_spatial_QC_across_batches()
    plt.savefig("plot_spatial_QC_across_batches.pdf", bbox_inches='tight')

    # Add 5% quantile of cell abundance to adata.obs for plotting
    adata_vis.obs[adata_vis.uns['mod']['factor_names']] = (
        adata_vis.obsm['q05_cell_abundance_w_sf']
    )

    # ── KNN + UMAP ──
    sc.pp.neighbors(
        adata_vis, use_rep='q05_cell_abundance_w_sf', n_neighbors=15
    )
    sc.tl.leiden(adata_vis, resolution=1.1)
    adata_vis.obs["region_cluster"] = adata_vis.obs["leiden"].astype("category")
    sc.tl.umap(adata_vis, min_dist=0.3, spread=1)

    # Plot regions in UMAP coordinates
    with mpl.rc_context({'axes.facecolor': 'white', 'figure.figsize': [8, 8]}):
        sc.pl.umap(
            adata_vis, color=['region_cluster'], size=30,
            color_map='RdPu', ncols=2, legend_loc='on data',
            legend_fontsize=20,
        )
        plt.savefig("spots_umap_knn.pdf", bbox_inches='tight')

        sc.pl.umap(
            adata_vis, color=['sample'], size=30,
            color_map='RdPu', ncols=2, legend_fontsize=20,
        )
        plt.savefig("spots_umap_sample.pdf", bbox_inches='tight')

    # Plot in spatial coordinates
    for library in libraries:
        with mpl.rc_context({'axes.facecolor': 'black', 'figure.figsize': [4.5, 5]}):
            sc.pl.spatial(
                adata_vis, color=['region_cluster'],
                size=1.3, img_key='hires', alpha=0.5, library_id=library,
            )
            plt.savefig(f"spatial_umap_knn_{library}.pdf", bbox_inches='tight')

    return adata_vis
