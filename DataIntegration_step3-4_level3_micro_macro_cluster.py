import anndata as ad
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import os
import pickle
from random import seed, sample
from collections import defaultdict
import glob
import pandas as pd


SEED = 1234
seed(SEED)

pd.set_option('display.max_columns', None)  # None means unlimited
pd.set_option('display.max_rows', None)     # None means unlimited


adata = ad.read_h5ad("../output_TACA2/immune_final2/integration_PROTEIN_CODING_reorganized_32bits_21datasets_micro_macro_batch_min10_HVG2000_epoch100_latent50_lr5e-4_all_nuclei_UMAP.h5ad")

sc.pp.neighbors(adata, use_rep = 'X_scVI')

sc.tl.louvain(adata, resolution = 1.0, random_state = SEED)

# sc.tl.leiden(adata, n_iterations = 2, resolution = 1.0, random_state = SEED)

adata.write_h5ad("../output_TACA2/immune_final2/integration_PROTEIN_CODING_reorganized_32bits_21datasets_micro_macro_batch_min10_HVG2000_epoch100_latent50_lr5e-4_all_nuclei_louvain_res10_UMAP.h5ad")


umap_coord = adata.obsm['X_umap']
# Plot the results
plt.figure(figsize=(9, 9))
sns.scatterplot(x=umap_coord[:,0], y=umap_coord[:, 1], hue=adata.obs['louvain'], s = 0.6, legend = False)
plt.xlabel('UMAP_1', fontsize = '16')
plt.ylabel('UMAP_2', fontsize = '16')
plt.xticks(fontsize = '16')
plt.yticks(fontsize = '16')
plt.savefig('../output_TACA2/immune_final2/micro_macro_l3_louvain_res10.png', bbox_inches = 'tight')


pre_adata = ad.read_h5ad("../output_TACA2/immune_final2/preintegration_PROTEIN_CODING_reorganized_32bits_21datasets_microglia_macrophage.h5ad")

assert adata.shape[0] == pre_adata.shape[0], 'Neuron nuclei # mismatch.'

pre_adata.obs['louvain'] = ''
pre_adata.obs.loc[adata.obs.index, 'louvain'] = adata.obs.loc[adata.obs.index, 'louvain']


sc.pp.normalize_total(pre_adata, target_sum=1e4)
sc.pp.log1p(pre_adata)


def downsample_nuclei(adata, meta_col, subsample_threshold):
    ct_count = adata.obs[meta_col].value_counts()
    ct_keep = set((ct_count[ct_count <= subsample_threshold]).index)
    ct_keep_nuclei = []
    for ct in ct_keep:
        ct_keep_nuclei += list(adata[adata.obs[meta_col] == ct].obs.index)
    ct_subsample =  set((ct_count[ct_count > subsample_threshold]).index)
    for ct in ct_subsample:
        ct_nuclei = list(adata[adata.obs[meta_col] == ct].obs.index)
        ct_sampled_nuclei = sample(ct_nuclei, subsample_threshold)
        ct_keep_nuclei += ct_sampled_nuclei
    sampled_adata = adata[ct_keep_nuclei]
    return sampled_adata


pre_adata = downsample_nuclei(pre_adata, meta_col = 'louvain', subsample_threshold = 20000)

sc.tl.rank_genes_groups(pre_adata, 'louvain', method='wilcoxon')

result = pre_adata.uns['rank_genes_groups']
deg_df = pd.DataFrame({
    group: result['names'][group]
    for group in result['names'].dtype.names
})
pvals_adj_df = pd.DataFrame({
    group: result['pvals_adj'][group]
    for group in result['pvals_adj'].dtype.names
})
logfc_df = pd.DataFrame({
    group: result['logfoldchanges'][group]
    for group in result['logfoldchanges'].dtype.names
})

collected_degs = defaultdict(set)
padj_threshold = 0.05
for cell_type in deg_df.columns:
    temp_df = pd.DataFrame({
        f'{cell_type}_gene': deg_df[cell_type],
        f'{cell_type}_padj': pvals_adj_df[cell_type],
        f'{cell_type}_logfc': logfc_df[cell_type]
    })
    filtered_temp_df = temp_df[(temp_df[f'{cell_type}_padj'] < padj_threshold)]
    top_temp_df = filtered_temp_df.sort_values(by=f'{cell_type}_logfc', ascending=False).head(500)
    top_temp_df.index = list(range(top_temp_df.shape[0]))
    collected_degs[f'{cell_type}'] = set(top_temp_df[f'{cell_type}_gene'])

dir_deg = '../output_TACA2/immune_final2/micro_macro_sampled_20000_top500_DEGs_epoch100_latent50_res10_louvain.pkl'
with open(dir_deg, 'wb') as file: 
    pickle.dump(collected_degs, file)


