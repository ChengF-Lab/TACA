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
    files = glob.glob(os.path.join(dir_folder, '*'+pattern))
    for dir_f in files:
        celltype = dir_f.split(pattern)[0].split('/')[-1]
        with open(dir_f, mode='r') as f:
            if header:
                next(f)
            for line in f:
                s, bid = line.strip("\n").split("\t")
                if bid in set(adata.var_names):
                    cellmarkers[celltype].append(bid)
    return cellmarkers



adata = ad.read_h5ad("../output_TACA2/immune_final2/integration_PROTEIN_CODING_reorganized_32bits_21datasets_immune_batch_min10_HVG2000_epoch100_latent50_lr5e-4_all_nuclei_UMAP.h5ad")

sc.pp.neighbors(adata, use_rep = 'X_scVI')

sc.tl.louvain(adata, resolution = 0.5, random_state = SEED)

# sc.tl.leiden(adata, n_iterations = 2, resolution = 2, random_state = SEED)

adata.write_h5ad("../output_TACA2/immune_final2/integration_PROTEIN_CODING_reorganized_32bits_21datasets_immune_batch_min10_HVG2000_epoch100_latent50_lr5e-4_all_nuclei_res05_UMAP.h5ad")

pre_adata = ad.read_h5ad("../output_TACA2/immune_final2/preintegration_PROTEIN_CODING_reorganized_32bits_21datasets_immune.h5ad")

assert adata.shape[0] == pre_adata.shape[0], 'Neuron nuclei # mismatch.'

pre_adata.obs['louvain'] = ''
pre_adata.obs.loc[adata.obs.index, 'louvain'] = adata.obs.loc[adata.obs.index, 'louvain']


sc.pp.normalize_total(pre_adata, target_sum=1e4)
sc.pp.log1p(pre_adata)

cell_markers = load_celltype_markers(dir_folder = '../data/immune_markers', pattern = '_ensembleID.txt', adata = pre_adata)

for ct in cell_markers:
    print("ct = {}, len(cell_markers[ct] = {}".format(ct, len(cell_markers[ct])))

plt.figure(figsize=(15, 9))
sc.pl.dotplot(pre_adata, cell_markers, groupby='louvain', dendrogram=False)
plt.savefig('../output_TACA2/immune_final2/immune_dotplot_L2_level.png', bbox_inches = 'tight')

