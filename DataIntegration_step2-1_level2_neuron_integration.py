import anndata as ad
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import scvi
from random import sample, seed
import os
import numpy as np
import re
import umap

SEED = 1234

seed(SEED)
scvi.settings.seed = SEED


pd.set_option('display.max_columns', None)  # None means unlimited
pd.set_option('display.max_rows', None)     # None means unlimited


# adata = ad.read_h5ad("../output_TACA2/integration_with_prediction_37824663_37824655_max_epoch_200_25_21_datasets.h5ad")
# pre_adata = ad.read_h5ad("../output_TACA2/preintegration_PROTEIN_CODING_reorganized_32bits_21datasets.h5ad")


# pre_adata.obs['C_scANVI'] = adata.obs['C_scANVI']
# neuron_adata = pre_adata[pre_adata.obs['C_scANVI'] == "Neuron"]

# neuron_adata.write_h5ad(os.path.join('../output_TACA2/neuron_final', 'preintegration_PROTEIN_CODING_reorganized_32bits_21datasets_neuron.h5ad'))


neuron_adata = ad.read_h5ad(os.path.join('../output_TACA2/neuron_final', 'preintegration_PROTEIN_CODING_reorganized_32bits_21datasets_neuron.h5ad'))


dir_hvg = '../output_TACA2/neuron_final/HVG_2000_reference_PROTEIN_CODING_reorganized_32bits_21datasets_neuron_batch_min10.txt'
if os.path.isfile(dir_hvg):
    HVGs = []
    with open(dir_hvg, mode='r') as f:
        for line in f:
            hvg = line.strip("\n").split("\t")[0]
            HVGs.append(hvg)

    HVGs = sorted(HVGs)
    neuron_adata = neuron_adata[:, HVGs]

    neuron_adata = neuron_adata.copy()

    scvi.model.SCVI.setup_anndata(neuron_adata, batch_key = 'BatchID',
                                    continuous_covariate_keys=["n_genes_by_counts", "total_counts", "pct_counts_mt", "pct_counts_rb"])

    model = scvi.model.SCVI(neuron_adata, n_latent = 50)

    model.train(max_epochs = 100, plan_kwargs={'lr':5e-4})

    latent = model.get_latent_representation()

    neuron_adata.obsm["X_scVI"] = latent

    neuron_adata.write_h5ad(os.path.join('../output_TACA2/neuron_final', 'integration_PROTEIN_CODING_reorganized_32bits_21datasets_neuron_batch_min10_HVG2000_epoch100_latent50_lr5e-4_all_nuclei.h5ad'))

    umap_coord = umap.UMAP(random_state = SEED)

    X_scVI_umap = umap_coord.fit_transform(neuron_adata.obsm['X_scVI'])

    neuron_adata.obsm['X_umap'] = X_scVI_umap

    neuron_adata.write_h5ad(os.path.join('../output_TACA2/neuron_final', 'integration_PROTEIN_CODING_reorganized_32bits_21datasets_neuron_batch_min10_HVG2000_epoch100_latent50_lr5e-4_all_nuclei_UMAP.h5ad'))
else:
    # # ====== add some preprcessing (sc.pp.filter_genes()) ====== #
    if np.any(neuron_adata.X.sum(axis = 0) == 0):
        print("There is gene(s) whose expressions are zeros in all nuclei!")
        neuron_adata = neuron_adata[:, (neuron_adata.X.sum(axis = 0) != 0)]

    batch_immune = neuron_adata.obs['BatchID'].value_counts()

    batch_immune_keep = set((batch_immune[batch_immune >= 10]).index)

    if not len(batch_immune_keep):
        raise Exception

    neuron_adata = neuron_adata[neuron_adata.obs['BatchID'].isin(batch_immune_keep)]

    sc.pp.highly_variable_genes(
        neuron_adata,
        flavor="seurat_v3",
        n_top_genes=2000,
        batch_key= 'BatchID',
        subset=True,
        span = 1.0
    )

    with open(dir_hvg, "w") as f_out:
        for g in neuron_adata.var['highly_variable'].index:
            f_out.write(g + "\n")







