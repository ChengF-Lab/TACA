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


def load_celltype_markers(dir_folder, pattern, adata, header = False):
    cellmarkers = defaultdict(list)
    files = glob.glob(os.path.join(dir_folder, pattern))
    for dir_f in files:
        celltype = dir_f.split('_unique.txt')[0].split('/')[-1]
        with open(dir_f, mode='r') as f:
            if header:
                next(f)
            for line in f:
                bid = line.strip("\n").split("\t")[0]
                if bid in set(adata.var_names):
                    cellmarkers[celltype].append(bid)
    for ct in cellmarkers:
        if len(cellmarkers[ct]) > 30:
            cellmarkers[ct] = cellmarkers[ct][:30]
    return cellmarkers



adata = ad.read_h5ad("../output_TACA2/neuron_final/integration_PROTEIN_CODING_reorganized_32bits_21datasets_neuron_batch_min10_HVG2000_epoch100_latent50_lr5e-4_all_nuclei_UMAP.h5ad")

sc.pp.neighbors(adata, use_rep = 'X_scVI')

sc.tl.louvain(adata, resolution = 2.0, random_state = SEED)

# sc.tl.leiden(adata, n_iterations = 2, resolution = 2, random_state = SEED)

adata.write_h5ad("../output_TACA2/neuron_final/integration_PROTEIN_CODING_reorganized_32bits_21datasets_neuron_batch_min10_HVG2000_epoch100_latent50_lr5e-4_all_nuclei_res20_UMAP.h5ad")


adata = ad.read_h5ad("../output_TACA2/neuron_final/integration_PROTEIN_CODING_reorganized_32bits_21datasets_neuron_batch_min10_HVG2000_epoch100_latent50_lr5e-4_all_nuclei_res20_UMAP.h5ad")
pre_adata = ad.read_h5ad("../output_TACA2/neuron_final/preintegration_PROTEIN_CODING_reorganized_32bits_21datasets_neuron.h5ad")

assert adata.shape[0] == pre_adata.shape[0], 'Neuron nuclei # mismatch.'

pre_adata.obs['louvain'] = ''
pre_adata.obs.loc[adata.obs.index, 'louvain'] = adata.obs.loc[adata.obs.index, 'louvain']


sc.pp.normalize_total(pre_adata, target_sum=1e4)
sc.pp.log1p(pre_adata)

cell_markers = load_celltype_markers(dir_folder = '../data/neuron_markers', pattern = '*_unique.txt', adata = pre_adata)

for ct in cell_markers:
    print("ct = {}, len(cell_markers[ct] = {}".format(ct, len(cell_markers[ct])))

cell_markers['substantia_nigra'] = ['ENSG00000165646', 'ENSG00000142319', 'ENSG00000180176', 'ENSG00000167281']

adata_37824663 = pre_adata[pre_adata.obs['Dataset'] == '37824663']

sampled_nuclei_37824663 = sample(list(adata_37824663.obs.index), int(0.1*len(adata_37824663.obs.index)))

sampled_adata_37824663 = adata_37824663[sampled_nuclei_37824663, :]

plt.figure(figsize=(15, 9))
sc.pl.dotplot(sampled_adata_37824663, cell_markers, groupby='louvain', dendrogram=False)
plt.savefig('../output_TACA2/neuron_final/integration_PROTEIN_CODING_reorganized_32bits_21datasets_neuron_latent50_res20_dotplot_37824663_sampled_nuclei.png', bbox_inches = 'tight')

sampled_nuclei = sample(list(adata.obs.index), int(0.1*len(adata.obs.index)))
sampled_adata = adata[sampled_nuclei, :]

plt.figure(figsize=(15, 9))
sc.pl.dotplot(sampled_adata, cell_markers, groupby='louvain', dendrogram=False)
plt.savefig('../output_TACA2/neuron_final/integration_PROTEIN_CODING_reorganized_32bits_21datasets_neuron_latent50_res20_dotplot_all_sampled_nuclei.png', bbox_inches = 'tight')

