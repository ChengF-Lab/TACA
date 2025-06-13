import anndata as ad
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import scvi
from random import sample, seed
import os
import numpy as np
import re
from collections import defaultdict
import pickle

SEED = 1234

seed(SEED)
scvi.settings.seed = SEED

pd.set_option('display.max_columns', None)  # None means unlimited
pd.set_option('display.max_rows', None)     # None means unlimited



inhib_adata = ad.read_h5ad("../output_TACA2/neuron_final/preintegration_PROTEIN_CODING_reorganized_32bits_21datasets_inhibitory_neuron_022425.h5ad")
inhib_adata.layers["counts"] = inhib_adata.X.copy()

dir_deg = '../output_TACA2/neuron_final/inhibneuron_37824655_top500_DEGs_1230.pkl'
dir_hvg = '../output_TACA2/neuron_final/inhibneuron_37824655_HVG_2000_reference_1230.txt' # check order of 2000 HVG

neuron_adata_subclass = ad.read_h5ad("/home/xuj2/isilon/Cheng-Qiu/TACA2/Data/37824655/processed/export/adata_subclass.h5ad")
common_index = set(inhib_adata.obs.index).intersection(set(neuron_adata_subclass.obs.index))
print(len(common_index))

inhib_adata.obs['subclass_label'] = ''
inhib_adata.obs.loc[list(common_index), 'subclass_label'] = neuron_adata_subclass.obs.loc[list(common_index), 'CrossArea_subclass']


reference_ds = inhib_adata[inhib_adata.obs['Dataset'] == '37824655']


reference_ds = reference_ds[reference_ds.obs['subclass_label'] != 'OPC']
reference_ds = reference_ds[reference_ds.obs['subclass_label'] != 'Oligo']
reference_ds = reference_ds[reference_ds.obs['subclass_label'] != 'Astro']
reference_ds = reference_ds[reference_ds.obs['subclass_label'] != 'Micro/PVM']
reference_ds = reference_ds[reference_ds.obs['subclass_label'] != 'Endo']
reference_ds = reference_ds[reference_ds.obs['subclass_label'] != 'VLMC']

reference_ds = reference_ds[reference_ds.obs['subclass_label'] != 'L2/3 IT']
reference_ds = reference_ds[reference_ds.obs['subclass_label'] != 'L5 IT']
reference_ds = reference_ds[reference_ds.obs['subclass_label'] != 'L4 IT']
reference_ds = reference_ds[reference_ds.obs['subclass_label'] != 'L6 IT']
reference_ds = reference_ds[reference_ds.obs['subclass_label'] != 'L6 CT']
reference_ds = reference_ds[reference_ds.obs['subclass_label'] != 'L6 IT Car3']
reference_ds = reference_ds[reference_ds.obs['subclass_label'] != 'L6b']
reference_ds = reference_ds[reference_ds.obs['subclass_label'] != 'L5/6 NP']
reference_ds = reference_ds[reference_ds.obs['subclass_label'] != 'L5 ET']



if os.path.isfile(dir_deg) and os.path.isfile(dir_hvg):

    inhib_adata.obs['seed_label'] = "Unknown"
    inhib_adata.obs.loc[reference_ds.obs.index, 'seed_label'] = inhib_adata.obs.loc[reference_ds.obs.index, 'subclass_label']
    print("inhib_adata.obs['seed_label'].value_counts() = ", inhib_adata.obs['seed_label'].value_counts())

    with open(dir_deg, "rb") as f:
        collected_degs = pickle.load(f)

    HVGs = []
    with open(dir_hvg, mode='r') as f:
        for line in f:
            hvg = line.strip("\n").split("\t")[0]
            HVGs.append(hvg)

    DEGs_HVGs = set()
    for ct in collected_degs:
        DEGs_HVGs |= collected_degs[ct]

    DEGs_HVGs |= set(HVGs)
    DEGs_HVGs = sorted(DEGs_HVGs)

    inhib_adata = inhib_adata[:, DEGs_HVGs]
    inhib_adata = inhib_adata.copy()

    # ==== each dataset should include batch id ==== #
    scvi.model.SCVI.setup_anndata(inhib_adata, batch_key = 'BatchID', labels_key = "seed_label", layer = "counts",
                                continuous_covariate_keys=["n_genes_by_counts", "total_counts", "pct_counts_mt", "pct_counts_rb"])

    model = scvi.model.SCVI(inhib_adata, n_latent = 50)

    model.module = model.module.float()

    model.train(max_epochs = 100, plan_kwargs={'lr':5e-4})

    scanvi_model = scvi.model.SCANVI.from_scvi_model(model, "Unknown")

    scanvi_model.module = scanvi_model.module.float()

    scanvi_model.train(25, plan_kwargs={'lr':5e-4})

    inhib_adata.obs["C_scANVI"] = scanvi_model.predict(inhib_adata)
    inhib_adata.obsm["X_scANVI"] = scanvi_model.get_latent_representation(inhib_adata)

    inhib_adata.write_h5ad(os.path.join('../output_TACA2/neuron_final/', 'integration_PROTEIN_CODING_reorganized_32bits_21datasets_inhib_neuron_raw_count_HVG2000_max_epoch100_epoch25_latent50_022625.h5ad'))
else:
    # ====== add some preprcessing ====== #
    if np.any(reference_ds.X.sum(axis = 0) == 0):
        print("There is gene(s) whose expressions are zeros in all nuclei!")
        reference_ds = reference_ds[:, (reference_ds.X.sum(axis = 0) != 0)]

    batch_neuron = reference_ds.obs['BatchID'].value_counts()
    batch_neuron_keep = set((batch_neuron[batch_neuron >= 10]).index)

    if not len(batch_neuron_keep):
        raise Exception

    reference_ds = reference_ds[reference_ds.obs['BatchID'].isin(batch_neuron_keep)]

    sc.pp.normalize_total(reference_ds, target_sum=1e4)
    sc.pp.log1p(reference_ds)

    # ====================== compute HVGs ====================== #
    sc.pp.highly_variable_genes(
        reference_ds,
        n_top_genes=2000,
        subset=True,
        layer="counts",
        flavor="seurat_v3",
        batch_key="BatchID",
        span = 1.0
    )

    with open(dir_hvg, "w") as f_out:
        for g in reference_ds.var['highly_variable'].index:
            f_out.write(g + "\n")

    # ====================== compute DEGs ====================== #
    sc.tl.rank_genes_groups(reference_ds, 'subclass_label', method='wilcoxon')
    result = reference_ds.uns['rank_genes_groups']
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
    with open(dir_deg, 'wb') as file:
        pickle.dump(collected_degs, file)











