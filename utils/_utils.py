import os
from collections import defaultdict
import pickle
import anndata as ad
from anndata import AnnData
import scanpy as sc
import matplotlib.pyplot as plt
import pandas as pd
import re
import numpy as np
from typing import Optional, Dict, Any
from random import sample, seed
import scvi
from sklearn.model_selection import train_test_split
import optuna
import logging
logger = logging.getLogger("TACA")


def _DEGs(
    adata, 
    dir_deg, 
    deg_group,
    batch_key,
    top_deg_num: int = 500,
    use_sampling: bool = False,
    sampling_ratio: float = 0.0
):

    logger.info("Computing Differentially Expressed Genes...")

    if np.any(adata.X.sum(axis = 0) == 0):
        logger.info("There is gene(s) whose expressions are zeros in all nuclei!")
        adata = adata[:, (adata.X.sum(axis = 0) != 0)]

    batch_adata = adata.obs[batch_key].value_counts()

    batch_adata_keep = set((batch_adata[batch_adata >= 10]).index)

    if not len(batch_adata_keep):
        raise Exception

    adata = adata[adata.obs[batch_key].isin(batch_adata_keep)]

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    if use_sampling:
        assert sampling_ratio > 0 and sampling_ratio < 1, 'Sampling ratio should be between 0 and 1.'
        sampled_nuclei = sample(list(adata.obs.index), int(sampling_ratio * len(adata.obs.index)))
        sampled_adata = adata[sampled_nuclei, :]
        sc.tl.rank_genes_groups(sampled_adata, deg_group, method='wilcoxon')
        result = sampled_adata.uns['rank_genes_groups']
    else:
        sc.tl.rank_genes_groups(adata, deg_group, method='wilcoxon')
        result = adata.uns['rank_genes_groups']

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
        top_temp_df = filtered_temp_df.sort_values(by=f'{cell_type}_logfc', ascending=False).head(top_deg_num)
        top_temp_df.index = list(range(top_temp_df.shape[0]))
        collected_degs[f'{cell_type}'] = set(top_temp_df[f'{cell_type}_gene'])

    with open(dir_deg, 'wb') as file:  # 'wb' mode means "write binary"
        pickle.dump(collected_degs, file)


def _HVGs(
    adata,
    dir_hvg,
    top_hvg_genes: int = 2000,
    subset: bool = True,
    layer: str = "counts",
    batch_key: str | None = None,
    span: float = 0.5
):
    logger.info("Computing Highly Variable Genes...")

    if np.any(adata.X.sum(axis = 0) == 0):
        logger.info("There is gene(s) whose expressions are zeros in all nuclei!")
        adata = adata[:, (adata.X.sum(axis = 0) != 0)]

    batch_adata = adata.obs[batch_key].value_counts()

    batch_adata_keep = set((batch_adata[batch_adata >= 10]).index)

    if not len(batch_adata_keep):
        raise Exception

    adata = adata[adata.obs[batch_key].isin(batch_adata_keep)]

    sc.pp.highly_variable_genes(
        adata,
        n_top_genes = top_hvg_genes,
        subset = subset,
        layer = layer,
        flavor="seurat_v3",
        batch_key= batch_key,
        span = span
    ) 

    with open(dir_hvg, "w") as f_out:
        for g in adata.var['highly_variable'].index:
            f_out.write(g + "\n")


def _subset_adata_DEGs_HVGs(adata, dir_hvg, dir_deg):

    assert os.path.isfile(dir_hvg), 'Provide the current HVG directories.'
    HVGs = set() # just run 2 times, and see whether they are with same order
    with open(dir_hvg, mode='r') as f:
        for line in f:
            g = line.strip("\n").split("\t")[0]
            HVGs.add(g)
    
    if dir_deg is not None:
        assert os.path.isfile(dir_deg), 'Provide the current DEG directories.'
        DEGS_superCT = set()
        with open(dir_deg, 'rb') as f:
            degs_CT = pickle.load(f)
            for ct in degs_CT:
                DEGS_superCT |= degs_CT[ct]
        HVG_DEGs_superCT = sorted(HVGs | DEGS_superCT)
    else:
        HVG_DEGs_superCT = sorted(HVGs)

    adata = adata[:, HVG_DEGs_superCT].copy()

    return adata



def _train_valid_test_split(
    adata, 
    seed_value: int,
    train_ratio: float | None = 0.6, 
    test_ratio: float | None = 0.5,
    shuffle: bool = True,
    adjust_cols: list[Any] = ['BatchID']
):
    #TODO: add K-folds
    assert train_ratio > 0 and train_ratio < 1, 'train_ratio should be between 0 and 1.'
    assert test_ratio > 0 and test_ratio < 1, 'test_ratio should be between 0 and 1.'

    idx = np.arange(adata.n_obs)
    train_idx, temp_idx = train_test_split(idx, test_size= (1 - train_ratio), random_state=seed_value, shuffle = shuffle)
    temp_val_idx, _ = train_test_split(temp_idx, test_size= test_ratio, random_state=seed_value, shuffle = shuffle)
    train_idx = sorted(train_idx)
    temp_val_idx = sorted(temp_val_idx)

    adata_train = adata[train_idx]
    adata_train = adata_train.copy()

    adata_val_temp  = adata[temp_val_idx]

    for col in adjust_cols:
        train_cats = adata_train.obs[col].astype("category").cat.categories
        adata_val_temp.obs[col] = pd.Categorical(adata_val_temp.obs[col], categories=train_cats)
        adata_val_temp = adata_val_temp[~adata_val_temp.obs[col].isna()].copy()

    adata_val = adata_val_temp
    adata_val = adata_val.copy()
    del adata_val_temp

    test_idx = sorted(set(temp_idx) - set(adata_val.obs.index))

    adata_test = adata[test_idx]
    adata_test = adata_test.copy()

    return adata_train, adata_val, adata_test


def _parameter_helper(
    trial: optuna.Trial,
    name: str,
    default: Any,
    *,
    continuous_low: float | None = None,
    continuous_high: float | None = None,
    categorical_choices: list[Any] | None = None,
    log: bool = False,
):
    """Suggest from trial if provided, else fall back to default."""
    if trial is None:
        return default
    if categorical_choices is not None:
        return trial.suggest_categorical(name, categorical_choices)
    if continuous_low is not None and continuous_high is not None:
        return trial.suggest_float(name, continuous_low, continuous_high, log=log)
    raise ValueError("hp spec incomplete")



def _scIntegration(
    adata, 
    level: int,
    dir_scvi:str,
    dir_scanvi:str,
    best_params, 
    batch_key,
    labels_key,
    prediction_key,
    layer,
    categorical_covariate_keys,
    continuous_covariate_keys,
    accelerator,
    scanvi_epochs: int = 25,
    use_unsupervised: bool = True
):
    logger.info('Start Integration ... ')
    scvi.model.SCVI.setup_anndata(
        adata,
        batch_key=batch_key,
        labels_key = labels_key,
        layer = layer,
        categorical_covariate_keys = categorical_covariate_keys,
        continuous_covariate_keys = continuous_covariate_keys
    )

    model = scvi.model.SCVI(adata, 
                            n_latent = best_params["n_latent"],
                            n_hidden = best_params["n_hidden"], 
                            n_layers = best_params["n_layers"], 
                            dispersion = best_params["dispersion"]
                        )

    model.module = model.module.float()

    model.train(
        max_epochs = best_params['max_epochs'], 
        plan_kwargs={'lr': best_params['lr']}, 
        accelerator=accelerator)

    model.save(dir_scvi, overwrite = True)

    if use_unsupervised:
        adata.obsm[f"Level{level}_X_scVI"] = model.get_latent_representation(adata)
    else:
        scanvi_model = scvi.model.SCANVI.from_scvi_model(model, "Unknown")
        scanvi_model.module = scanvi_model.module.float()
        scanvi_model.train(scanvi_epochs, plan_kwargs={'lr':best_params['lr']}, accelerator=accelerator)
        scanvi_model.save(dir_scanvi, overwrite=True)
        adata.obs[f"Level{level}_{prediction_key}"] = scanvi_model.predict(adata)
        adata.obsm[f"Level{level}_X_scANVI"] = scanvi_model.get_latent_representation(adata)
        scanvi_model.save(dir_scanvi, overwrite = True)
    return adata



def _hyperparameter_tune(
    adata_train, 
    adata_val, 
    n_latent,
    n_hidden,
    n_layers,
    dispersion,
    lr,
    max_epochs,
    accelerator: str = 'cpu',
    labels_key: str | None = None,
    layer: str | None = None,
    batch_key: str | None = None,
    categorical_covariate_keys: list[Any] | None = None,
    continuous_covariate_keys: list[Any] | None = None,
    trial: Optional[optuna.Trial] = None,
    **kwargs
):

    n_latent = _parameter_helper(trial = trial, name = 'n_latent', **n_latent)
    n_hidden = _parameter_helper(trial = trial, name = 'n_hidden', **n_hidden)
    n_layers = _parameter_helper(trial = trial, name = 'n_layers', **n_layers)
    dispersion = _parameter_helper(trial = trial, name = 'dispersion', **dispersion)
    lr = _parameter_helper(trial = trial, name = 'lr', **lr)
    max_epochs = _parameter_helper(trial = trial, name = 'max_epochs', **max_epochs)


    scvi.model.SCVI.setup_anndata(
        adata_train,
        batch_key=batch_key,
        labels_key = labels_key,
        layer = layer,
        categorical_covariate_keys = categorical_covariate_keys,
        continuous_covariate_keys = continuous_covariate_keys
    )

    scvi.model.SCVI.setup_anndata(
        adata_val,
        batch_key=batch_key,
        labels_key = labels_key,
        layer = layer,
        categorical_covariate_keys = categorical_covariate_keys,
        continuous_covariate_keys = continuous_covariate_keys
    )

    model = scvi.model.SCVI(
        adata_train,
        n_latent = n_latent,
        n_hidden = n_hidden,
        n_layers = n_layers,
        dispersion = dispersion,
    )

    model.train(
        max_epochs = max_epochs,
        plan_kwargs = {"lr": lr},
        accelerator = accelerator, 
        check_val_every_n_epoch = None,
        early_stopping = False,
        train_size = 1.0
    )

    elbo_val = model.get_elbo(adata_val)
    return elbo_val  



def _meta2dict(df, prop):
    samples2meta = {}
    if prop not in df.columns:
        # Assign 'Unknown' to all SampleID if the property doesn't exist in the DataFrame
        return {sid: 'Unknown' for sid in set(df['SampleID'])}
    # Define classification mappings
    apoe_e4_set = {'E2E4', 'E3E4', 'E4E4'}
    apoe_non_e4_set = {'E2E2', 'E2E3', 'E3E3'}
    # Iterate over the DataFrame directly instead of using multiple loops
    for _, row in df.iterrows():
        sid = row['SampleID']
        value = row[prop]
        if prop == 'APOE':
            samples2meta[sid] = 'E4' if value in apoe_e4_set else ('non-E4' if value in apoe_non_e4_set else 'Unknown')
        elif prop == 'APOE4 status':
            samples2meta[sid] = 'E4' if value == 'Y' else ('non-E4' if value == 'N' else 'Unknown')
        else:
            samples2meta[sid] = value  # Keep original categorical values for other props
    return samples2meta