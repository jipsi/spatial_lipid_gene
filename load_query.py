import seaborn as sns
# silence scanpy that prints a lot of warnings
import warnings
warnings.filterwarnings('ignore')

def read_and_qc(sample_name, path):
    r""" This function reads the data for one 10X spatial experiment into the anndata object.
    It also calculates QC metrics. Modify this function if required by your workflow.

    :param sample_name: Name of the sample
    :param path: path to data
    """
    #adding library_id to make sure key is consistent with sample name
    #for later methods like select_slide()
    adata = sc.read_visium(path,
                           count_file=f"{sample_name}.h5", library_id=sample_name, load_images=True)
    adata.obs['sample'] = sample_name
    adata.var['SYMBOL'] = adata.var_names
    print(adata.var['SYMBOL'])

    # Calculate QC metrics
    from scipy.sparse import csr_matrix
    adata.X = adata.X.toarray()
    sc.pp.calculate_qc_metrics(adata, inplace=True)
    adata.X = csr_matrix(adata.X)
    #The following lines are not applicable to Visium FFPE as ribosomal and mitochondrial
    #genes are not included
    adata.var['MT'] = [gene.startswith('mt-') for gene in adata.var['SYMBOL']]
    adata.obs['mt_frac'] = adata[:, adata.var['MT'].tolist()].X.sum(1).A.squeeze()/adata.obs['total_counts']
    print(adata.obs['mt_frac'])

    # add sample name to obs names
    adata.obs["sample"] = [str(i) for i in adata.obs['sample']]

    adata.obs_names = adata.obs["sample"] \
                          + '_' + adata.obs_names
                     
    adata.obs.index.name = 'spot_id'
    adata.var_names_make_unique()
 
    return adata

def select_slide(adata, s, s_col='sample'):
    r""" This function selects the data for one slide from the spatial anndata object.

    :param adata: Anndata object with multiple spatial experiments
    :param s: name of selected experiment
    :param s_col: column in adata.obs listing experiment name for each location
    """

    slide = adata[adata.obs[s_col].isin([s]), :]
    s_keys = list(slide.uns['spatial'].keys())
    
    s_spatial = np.array(s_keys)[[s in k for k in s_keys]][0]

    slide.uns['spatial'] = {s_spatial: slide.uns['spatial'][s_spatial]}

    return slide

#######################

# Read the list of spatial experiments
sample_data = pd.read_csv(sample_list)
    
from csv import DictReader
#read the data into anndata objects
slides = []
# open file in read mode
with open(sample_list, 'r') as read_obj:
    csv_dict_reader = DictReader(read_obj)
    for row in csv_dict_reader:
        print(row['sample_name'], row['sample_dir'])
        sample = row['sample_name']
        sample_dir = sp_data_folder + row['sample_dir']
        slides.append(read_and_qc(sample, path=sample_dir))
    
    

# Combine anndata objects together
adata_vis = slides[0].concatenate(
    slides[1:],
    batch_key="sample",
    uns_merge="unique",
    batch_categories=sample_data['sample_name'],
    index_unique=None
)

#create name for sample run as this is differerent from sample
#sample_run name = 30i_first_run
#slide = adata[adata.obs[s_col].isin([s]), :]
#s_keys = list(slide.uns['spatial'].keys())
#vis_slide_key=sample_data['sample_name'][0].split('_')[0]
#new_s_keys = list(f"{vis_slide_key}_{s_keys[0].split('_')[1]}" f"{vis_slide_key}_{s_keys[1].split('_')[1]}" f"{vis_slide_key}_{s_keys[2].split('_')[1]}")
#######################

# PLOT QC FOR EACH SAMPLE
fig, axs = plt.subplots(len(slides), 4, figsize=(15, 4*len(slides)-4))
for i, s in enumerate(adata_vis.obs['sample'].unique()):
    #fig.suptitle('Covariates for filtering')

    slide = select_slide(adata_vis, s)
    sns.distplot(slide.obs['total_counts'],
                 kde=False, ax = axs[i, 0])
    axs[i, 0].set_xlim(0, adata_vis.obs['total_counts'].max())
    axs[i, 0].set_xlabel(f'total_counts | {s}')

    sns.distplot(slide.obs['total_counts']\
                 [slide.obs['total_counts']<20000],
                 kde=False, bins=40, ax = axs[i, 1])
    axs[i, 1].set_xlim(0, 20000)
    axs[i, 1].set_xlabel(f'total_counts | {s}')

    sns.distplot(slide.obs['n_genes_by_counts'],
                 kde=False, bins=60, ax = axs[i, 2])
    axs[i, 2].set_xlim(0, adata_vis.obs['n_genes_by_counts'].max())
    axs[i, 2].set_xlabel(f'n_genes_by_counts | {s}')

    sns.distplot(slide.obs['n_genes_by_counts']\
                 [slide.obs['n_genes_by_counts']<6000],
                 kde=False, bins=60, ax = axs[i, 3])
    axs[i, 3].set_xlim(0, 6000)
    axs[i, 3].set_xlabel(f'n_genes_by_counts | {s}')

plt.savefig(f"qc_vis_samples.pdf", bbox_inches='tight')

exec(open("deconvolute.py").read())

