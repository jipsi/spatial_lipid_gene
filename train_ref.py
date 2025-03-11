
#load the data
adata_file = f"{ref_run_name}/sc_source.h5ad"
print(adata_file)
adata_ref = sc.read_h5ad(adata_file)

################################################################################################################
#When converting from Seurat objects need to QC first by removing/renaming all columns named as "_index" as reserve word in python
print(adata_ref.__dict__['_raw'].__dict__['_var'])
adata_ref.__dict__['_raw'].__dict__['_var'] = adata_ref.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})
print(adata_ref.__dict__['_raw'].__dict__['_var'])
#When converting from Seurat objects need to QC first by removing/renaming all columns named as "_index" as reserve word in python
#print(adata_ref.var)
#del(adata_ref.var['_index'])
#print(adata_ref.var)
adata_ref.var.set_index(adata_ref.var.features, inplace=True)
adata_ref.var.index.name = 'index'
###############################################################################################################

# prepare anndata for the regression model
cell2location.models.RegressionModel.setup_anndata(adata=adata_ref,
                        # 10X reaction / sample / batch
                        batch_key='orig.ident',
                        # cell type, covariate used for constructing signatures
                        labels_key='CellMeshFine'
                       )


# create the regression model
from cell2location.models import RegressionModel
mod = RegressionModel(adata_ref)

# view anndata_setup as a sanity check
mod.view_anndata_setup()

mod.train(max_epochs=500, use_gpu=True)

mod.plot_history(20)
plt.savefig('ref.mod.history.pdf', bbox_inches='tight')

# In this section, we export the estimated cell abundance (summary of the posterior distribution).
adata_ref = mod.export_posterior(
    adata_ref, sample_kwargs={'num_samples': 1000, 'batch_size': 2500, 'use_gpu': True}
)


#mod qc
mod.plot_QC()
plt.savefig('ref.mod.qc.pdf', bbox_inches='tight')

# Save model
mod.save(f"{ref_run_name}", overwrite=True)

# Save anndata object with results
adata_file = f"{ref_run_name}/sc.h5ad"
adata_ref.write(adata_file)
adata_file

#Reload the data
adata_file = f"{ref_run_name}/sc.h5ad"
adata_ref = sc.read_h5ad(adata_file)

adata_ref.var.set_index(adata_ref.var.features, inplace=True)
adata_ref.var.index.name = 'index'

mod = cell2location.models.RegressionModel.load(f"{ref_run_name}", adata_ref)
    

# export estimated expression in each cluster
if 'means_per_cluster_mu_fg' in adata_ref.varm.keys():
    inf_aver = adata_ref.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_ref.uns['mod']['factor_names']]].copy()
else:
    inf_aver = adata_ref.var[[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_ref.uns['mod']['factor_names']]].copy()
inf_aver.columns = adata_ref.uns['mod']['factor_names']
inf_aver.iloc[0:5, 0:5]

exec(open("load_query.py").read())