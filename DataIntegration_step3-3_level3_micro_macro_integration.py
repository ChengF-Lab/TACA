import anndata as ad
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import os
import pickle
import numpy as np
import pandas as pd
import scvi
import umap

pd.set_option('display.max_columns', None)  # None means unlimited
pd.set_option('display.max_rows', None)     # None means unlimited

SEED = 1234
scvi.settings.seed = SEED

micro_macro_adata = ad.read_h5ad("../output_TACA2/immune_final2/preintegration_PROTEIN_CODING_reorganized_32bits_21datasets_microglia_macrophage.h5ad")

dir_hvg = "../output_TACA2/immune_final2/HVG_2000_PROTEIN_CODING_micro_macro_batch_min10.txt"

if os.path.isfile(dir_hvg):
    HVGs = []
    with open(dir_hvg, mode='r') as f:
        for line in f:
            hvg = line.strip("\n").split("\t")[0]
            HVGs.append(hvg)

    HVGs = sorted(HVGs)
    micro_macro_adata = micro_macro_adata[:, HVGs]

    micro_macro_adata = micro_macro_adata.copy()


    scvi.model.SCVI.setup_anndata(micro_macro_adata, batch_key = 'BatchID',
                                    continuous_covariate_keys=["n_genes_by_counts", "total_counts", "pct_counts_mt", "pct_counts_rb"])

    model = scvi.model.SCVI(micro_macro_adata, n_latent = 50)

    model.train(max_epochs = 100, plan_kwargs={'lr':5e-4})

    latent = model.get_latent_representation()

    micro_macro_adata.obsm["X_scVI"] = latent

    micro_macro_adata.write_h5ad(os.path.join('../output_TACA2/immune_final2', 'integration_PROTEIN_CODING_reorganized_32bits_21datasets_micro_macro_batch_min10_HVG2000_epoch100_latent50_lr5e-4_all_nuclei.h5ad'))

    umap_coord = umap.UMAP(random_state = SEED)

    X_scVI_umap = umap_coord.fit_transform(micro_macro_adata.obsm['X_scVI'])

    micro_macro_adata.obsm['X_umap'] = X_scVI_umap

    micro_macro_adata.write_h5ad(os.path.join('../output_TACA2/immune_final2', 'integration_PROTEIN_CODING_reorganized_32bits_21datasets_micro_macro_batch_min10_HVG2000_epoch100_latent50_lr5e-4_all_nuclei_UMAP.h5ad'))
else:
    if np.any(micro_macro_adata.X.sum(axis = 0) == 0):
        print("There is gene(s) whose expressions are zeros in all nuclei!")
        micro_macro_adata = micro_macro_adata[:, (micro_macro_adata.X.sum(axis = 0) != 0)]

    batch_micro_macro = micro_macro_adata.obs['BatchID'].value_counts()

    batch_micro_macro_keep = set((batch_micro_macro[batch_micro_macro >= 10]).index)

    if not len(batch_micro_macro_keep):
        raise Exception

    micro_macro_adata = micro_macro_adata[micro_macro_adata.obs['BatchID'].isin(batch_micro_macro_keep)]

    sc.pp.highly_variable_genes(
        micro_macro_adata,
        flavor="seurat_v3",
        n_top_genes=2000,
        batch_key= 'BatchID',
        subset=True,
        span = 1.0
    )

    with open(dir_hvg, "w") as f_out:
        for g in micro_macro_adata.var['highly_variable'].index:
            f_out.write(g + "\n")




