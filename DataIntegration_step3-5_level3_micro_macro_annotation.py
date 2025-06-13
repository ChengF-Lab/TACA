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


adata = ad.read_h5ad("../output_TACA2/immune_final2/integration_PROTEIN_CODING_reorganized_32bits_21datasets_micro_macro_batch_min10_HVG2000_epoch100_latent50_lr5e-4_all_nuclei_louvain_res10_UMAP.h5ad")
pre_adata = ad.read_h5ad("../output_TACA2/immune_final2/preintegration_PROTEIN_CODING_reorganized_32bits_21datasets_microglia_macrophage.h5ad")

assert adata.shape[0] == pre_adata.shape[0], 'Neuron nuclei # mismatch.'

pre_adata.obs['louvain'] = ''
pre_adata.obs.loc[adata.obs.index, 'louvain'] = adata.obs.loc[adata.obs.index, 'louvain']


adata.obs['immune_L3'] = ''
adata.obs['immune_L3'] = adata.obs['louvain']
adata.obs['immune_L3'] = adata.obs['immune_L3'].astype(str)


adata.obs.loc[adata.obs['louvain'] == '0', 'immune_L3'] = 'Tau Microglia'
adata.obs.loc[adata.obs['louvain'] == '1', 'immune_L3'] = 'Homeostasis Microglia1'
adata.obs.loc[adata.obs['louvain'] == '2', 'immune_L3'] = 'Homeostasis Microglia1'
adata.obs.loc[adata.obs['louvain'] == '3', 'immune_L3'] = 'DAM1'
adata.obs.loc[adata.obs['louvain'] == '4', 'immune_L3'] = 'Homeostasis Microglia2'
adata.obs.loc[adata.obs['louvain'] == '5', 'immune_L3'] = 'MHC Microglia1'
adata.obs.loc[adata.obs['louvain'] == '6', 'immune_L3'] = 'Unknown1'
adata.obs.loc[adata.obs['louvain'] == '7', 'immune_L3'] = 'DAM2'
adata.obs.loc[adata.obs['louvain'] == '8', 'immune_L3'] = 'Macrophage'
adata.obs.loc[adata.obs['louvain'] == '9', 'immune_L3'] = 'Unknown1'
adata.obs.loc[adata.obs['louvain'] == '10', 'immune_L3'] = 'MHC Microglia2'
adata.obs.loc[adata.obs['louvain'] == '11', 'immune_L3'] = 'Homeostasis Microglia1'
adata.obs.loc[adata.obs['louvain'] == '12', 'immune_L3'] = 'Inflammation Microglia'
adata.obs.loc[adata.obs['louvain'] == '13', 'immune_L3'] = 'Proliferation Microglia'
adata.obs.loc[adata.obs['louvain'] == '14', 'immune_L3'] = 'Monocyte'
adata.obs.loc[adata.obs['louvain'] == '15', 'immune_L3'] = 'Unknown2'
adata.obs.loc[adata.obs['louvain'] == '16', 'immune_L3'] = 'Unknown3'
adata.obs.loc[adata.obs['louvain'] == '17', 'immune_L3'] = 'MHC Microglia3'
adata.obs.loc[adata.obs['louvain'] == '18', 'immune_L3'] = 'Proliferation Microglia'


umap_coord = adata.obsm['X_umap']
# Plot the results
plt.figure(figsize=(9, 9))
sns.scatterplot(x=umap_coord[:,0], y=umap_coord[:, 1], hue=adata.obs['immune_L3'], s = 0.6, legend = False)
plt.xlabel('UMAP_1', fontsize = '16')
plt.ylabel('UMAP_2', fontsize = '16')
plt.xticks(fontsize = '16')
plt.yticks(fontsize = '16')
plt.savefig('../output_TACA2/immune_final2/micro_macro_L3_label.png', bbox_inches = 'tight')

adata.write_h5ad('../output_TACA2/immune_final2/integration_PROTEIN_CODING_reorganized_32bits_21datasets_micro_macro_batch_min10_HVG2000_epoch100_latent50_lr5e-4_all_nuclei_louvain_res10_UMAP_L3_label.h5ad')


pre_adata.obs['immune_L3'] = ''
pre_adata.obs.loc[adata.obs.index, 'immune_L3'] = adata.obs.loc[adata.obs.index, 'immune_L3']
pre_adata.write_h5ad('../output_TACA2/immune_final2/preintegration_PROTEIN_CODING_reorganized_32bits_21datasets_micro_macro_L3_label.h5ad')


pre_adata = ad.read_h5ad('../output_TACA2/immune_final2/preintegration_PROTEIN_CODING_reorganized_32bits_21datasets_micro_macro_L3_label.h5ad')
immune_L2 = ad.read_h5ad('../output_TACA2/immune_final2/preintegration_PROTEIN_CODING_reorganized_32bits_21datasets_immune_L2_label.h5ad')
immune_L2.obs['immune_L3'] = ''
immune_L2.obs['immune_L3'] = immune_L2.obs['immune_L2']
pre_adata.obs['immune_L3'] = pre_adata.obs['immune_L3'].astype(str)
immune_L2.obs['immune_L3'] = immune_L2.obs['immune_L3'].astype(str)
immune_L2.obs.loc[pre_adata.obs.index, 'immune_L3'] = pre_adata.obs.loc[pre_adata.obs.index, 'immune_L3']

immune_L2.write_h5ad('../output_TACA2/immune_final2/preintegration_PROTEIN_CODING_reorganized_32bits_21datasets_immune_L2_L3_label.h5ad')
