import anndata as ad
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import scvi
from random import sample, seed
import os
import numpy as np
import pickle
from collections import defaultdict
import umap

SEED = 1234

seed(SEED)
scvi.settings.seed = SEED

adata = ad.read_h5ad("../output_TACA2/preintegration_PROTEIN_CODING_reorganized_32bits_21datasets.h5ad")

# ====== Create Seed Labels ====== #
'''
37824655 = TACA24CXG08
37824663 = TACA24CXG09
'''
adata.obs['seed_label'] = "Unknown"

adata.obs.loc[adata.obs.index.str.startswith(('TACA24CXG08', 'TACA24CXG09')), 'seed_label'] = adata.obs.loc[adata.obs.index.str.startswith(('TACA24CXG08', 'TACA24CXG09')), 'Super_Celltype']

# ======== Start integration ========= #
adata.layers["counts"] = adata.X.copy()

dir_deg = '../output_TACA2/top500_deg_reference_sampled_10percent_PROTEIN_CODING_reorganized_32bits_21datasets.pkl'
dir_hvg = '../output_TACA2/HVG_2000_reference_PROTEIN_CODING_reorganized_32bits_21datasets_123024.txt'



if os.path.isfile(dir_deg) and os.path.isfile(dir_hvg):

# if os.path.isfile(dir_hvg_deg):
#     HVG_DEGs_superCT = []
#     with open(dir_hvg_deg, mode='r') as f:
#         for line in f:
#             g = line.strip("\n").split("\t")[0]
#             HVG_DEGs_superCT.append(g)

    #both HVG and DEGs are from reference data only

    with open(dir_deg, 'rb') as f:  # 'rb' = read in binary mode
        degs_CT = pickle.load(f)

    DEGS_superCT = set()
    for ct in degs_CT:
        DEGS_superCT |= degs_CT[ct]

    HVGs = []
    with open(dir_hvg, mode='r') as f:
        for line in f:
            g = line.strip("\n").split("\t")[0]
            HVGs.append(g)

    HVG_DEGs_superCT = sorted(set(HVGs) | DEGS_superCT)

    adata = adata[:, HVG_DEGs_superCT]

    print(adata.shape)

    adata = adata.copy()

    # ==== each dataset should include batch id ==== #
    scvi.model.SCVI.setup_anndata(adata, batch_key = 'BatchID', labels_key = "seed_label", layer = "counts",
                                continuous_covariate_keys=["n_genes_by_counts", "total_counts", "pct_counts_mt", "pct_counts_rb"])

    model = scvi.model.SCVI(adata, n_latent = 30)

    model.module = model.module.float()

    model.train(max_epochs = 200, plan_kwargs={'lr':5e-4})

    scanvi_model = scvi.model.SCANVI.from_scvi_model(model, "Unknown")

    scanvi_model.module = scanvi_model.module.float()

    scanvi_model.train(25, plan_kwargs={'lr':5e-4})

    adata.obs["C_scANVI"] = scanvi_model.predict(adata)
    adata.obsm["X_scANVI"] = scanvi_model.get_latent_representation(adata)

    adata.write_h5ad(os.path.join('../output_TACA2/', 'integration_with_prediction_37824663_37824655_max_epoch_200_25_21_datasets.h5ad'))

    umap_coord = umap.UMAP(random_state = SEED)

    X_scANVI_umap = umap_coord.fit_transform(adata.obsm['X_scANVI'])

    adata.obsm['X_umap'] = X_scANVI_umap

    adata.write_h5ad(os.path.join('../output_TACA2/', 'integration_with_prediction_37824663_37824655_max_epoch_200_25_21_datasets_UMAP.h5ad'))
else:
    nuclei_37824655 = adata[adata.obs['Dataset'] == '37824655'].obs.index
    nuclei_37824663 = adata[adata.obs['Dataset'] == '37824663'].obs.index
    reference_nuclei = list(nuclei_37824655) + list(nuclei_37824663)
    reference_adata = adata[reference_nuclei, :]

    sc.pp.normalize_total(reference_adata, target_sum=1e4)
    sc.pp.log1p(reference_adata)

    sampled_nuclei_10_ref = sample(reference_nuclei, int(0.1*len(reference_nuclei)))
    sampled_10_ref_adata = reference_adata[sampled_nuclei_10_ref, :]

    sc.tl.rank_genes_groups(sampled_10_ref_adata, 'Super_Celltype', method='wilcoxon')

    result = sampled_10_ref_adata.uns['rank_genes_groups']

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


    # check repeated SampleID from CxG, datasetID_sampleID
    sc.pp.highly_variable_genes(
        reference_adata,
        n_top_genes=2000,
        subset=True,
        layer="counts",
        flavor="seurat_v3",
        batch_key="BatchID",
        span = 0.5
    ) 

    with open(dir_hvg, "w") as f_out:
        for g in reference_ds.var['highly_variable'].index:
            f_out.write(g + "\n")



