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


adata = ad.read_h5ad("../output_TACA2/immune_final2/integration_PROTEIN_CODING_reorganized_32bits_21datasets_immune_batch_min10_HVG2000_epoch100_latent50_lr5e-4_all_nuclei_res05_UMAP.h5ad")
pre_adata = ad.read_h5ad("../output_TACA2/immune_final2/preintegration_PROTEIN_CODING_reorganized_32bits_21datasets_immune.h5ad")

assert adata.shape[0] == pre_adata.shape[0], 'Neuron nuclei # mismatch.'

pre_adata.obs['louvain'] = ''
pre_adata.obs.loc[adata.obs.index, 'louvain'] = adata.obs.loc[adata.obs.index, 'louvain']



adata.obs['immune_L2'] = 'Unknown'

adata.obs.loc[adata.obs['louvain'] == '0', 'immune_L2'] = 'Microglia_Macrophage'
adata.obs.loc[adata.obs['louvain'] == '1', 'immune_L2'] = 'Microglia_Macrophage'
adata.obs.loc[adata.obs['louvain'] == '2', 'immune_L2'] = 'Microglia_Macrophage'
adata.obs.loc[adata.obs['louvain'] == '3', 'immune_L2'] = 'Microglia_Macrophage'
adata.obs.loc[adata.obs['louvain'] == '4', 'immune_L2'] = 'Microglia_Macrophage'
adata.obs.loc[adata.obs['louvain'] == '5', 'immune_L2'] = 'Microglia_Macrophage'
adata.obs.loc[adata.obs['louvain'] == '6', 'immune_L2'] = 'Leukocyte_T cell'
adata.obs.loc[adata.obs['louvain'] == '7', 'immune_L2'] = 'Microglia_Macrophage'
adata.obs.loc[adata.obs['louvain'] == '8', 'immune_L2'] = 'Microglia_Macrophage'
adata.obs.loc[adata.obs['louvain'] == '9', 'immune_L2'] = 'Microglia_Macrophage'
adata.obs.loc[adata.obs['louvain'] == '10', 'immune_L2'] = 'Microglia_Macrophage'
adata.obs.loc[adata.obs['louvain'] == '11', 'immune_L2'] = 'Microglia_Macrophage'


adata.write_h5ad('../output_TACA2/immune_final2/integration_PROTEIN_CODING_reorganized_32bits_21datasets_immune_batch_min10_HVG2000_epoch100_latent50_lr5e-4_all_nuclei_res05_UMAP_L2_label.h5ad')


adata = ad.read_h5ad('../output_TACA2/immune_final2/integration_PROTEIN_CODING_reorganized_32bits_21datasets_immune_batch_min10_HVG2000_epoch100_latent50_lr5e-4_all_nuclei_res05_UMAP_L2_label.h5ad')

umap_coord = adata.obsm['X_umap']
plt.figure(figsize=(9, 9))
sns.scatterplot(x=umap_coord[:,0], y=umap_coord[:, 1], hue=adata.obs['immune_L2'], s = 0.6, legend = False)
plt.xlabel('UMAP_1', fontsize = '16')
plt.ylabel('UMAP_2', fontsize = '16')
plt.xticks(fontsize = '16')
plt.yticks(fontsize = '16')
plt.savefig('../output_TACA2/immune_final2/immune_umap_L2_level.png', bbox_inches = 'tight')


pre_adata.obs['immune_L2'] = ''
pre_adata.obs.loc[adata.obs.index, 'immune_L2'] = adata.obs.loc[adata.obs.index, 'immune_L2']
pre_adata.write_h5ad('../output_TACA2/immune_final2/preintegration_PROTEIN_CODING_reorganized_32bits_21datasets_immune_L2_label.h5ad')

micro_macro_adata = pre_adata[pre_adata.obs['immune_L2'] == 'Microglia_Macrophage']
micro_macro_adata.write_h5ad(os.path.join('../output_TACA2/immune_final2', 'preintegration_PROTEIN_CODING_reorganized_32bits_21datasets_microglia_macrophage.h5ad'))



