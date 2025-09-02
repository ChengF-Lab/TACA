import os
import sys
import anndata as ad
from anndata import AnnData
import scanpy as sc
import matplotlib.pyplot as plt
from sklearn.metrics import confusion_matrix, precision_score, recall_score
import re
import numpy as np
from typing import Optional, Dict, Any
from random import sample, seed
import scvi
from sklearn.model_selection import train_test_split
from ._utils import _DEGs, _HVGs, _subset_adata_DEGs_HVGs, _meta2dict, \
    _train_valid_test_split, _parameter_helper, _scIntegration, _hyperparameter_tune
from .constants import *
import optuna
import logging
logger = logging.getLogger("TACA")
from logging.handlers import WatchedFileHandler
import umap
import pandas as pd
pd.set_option('display.max_columns', None)  # None means unlimited
pd.set_option('display.max_rows', None)     # None means unlimited




def setup_logger(name="TACA", log_file=None, log_dir=None, level=logging.INFO):

    logger = logging.getLogger(name)
    logger.setLevel(level)

    # --- Clear old handlers (important in Jupyter) ---
    if logger.hasHandlers():
        logger.handlers.clear()

    # --- Console handler (always enabled) ---
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(level)
    formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s",
                                  datefmt="%H:%M:%S")
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)

    # --- Optional File handler ---
    if log_file is not None:
        if log_dir is not None:
            os.makedirs(log_dir, exist_ok=True)
            log_path = os.path.join(log_dir, log_file)
        else:
            log_path = log_file
        file_handler = WatchedFileHandler(log_path)
        file_handler.setLevel(level)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

    return logger


def create_directory(
    dir_job 
):
    if not os.path.exists(dir_job):
        logger.info("Create Directory '{}' .".format(os.path.abspath(dir_job)))
        os.makedirs(dir_job)
    else:
        logger.info("Directory '{}' Has Been Created Already.".format(os.path.abspath(dir_job)))


def setup_reference_label(
    adata,
    labels_key,
    ref_label,
    celltype,
    dataset_key: str = 'Dataset',
    reference_datasets: list[str] | None = None,
    dir_reference: str | None = None
):
    logger.info("Importing Seed Labels from Reference Datasets ...")
    adata.obs[labels_key] = "Unknown"
    if reference_datasets is not None:
        for ds in set(adata.obs['Dataset']):
            if ds in reference_datasets:
                ref_nuclei = adata[adata.obs['Dataset'] == ds].obs.index
                adata.obs.loc[ref_nuclei, labels_key] = adata.obs.loc[ref_nuclei, ref_label]
    elif dir_reference is not None:
        ref_adata = ad.read_h5ad(dir_reference)
        common_index = set(adata.obs.index).intersection(set(ref_adata.obs.index))
        if len(common_index) > 0:
            adata.obs.loc[list(common_index), labels_key] = ref_adata.obs.loc[list(common_index), ref_label]
    else:
        raise ValueError('Provide either reference_datasets or dir_reference.')

    if celltype in SUBCTs:
        allowed_ref_labels = SUBCTs[celltype]
        allowed_ref_labels.add("Unknown")
        if len(set(adata.obs[labels_key]) - allowed_ref_labels) > 0:
            for ct in set(adata.obs[labels_key]) - allowed_ref_labels:
                ct_nuclei = adata[adata.obs[labels_key] == ct].obs.index
                adata.obs.loc[ct_nuclei, labels_key] = 'Unknown'
    return adata



def DEGs_HVGs(
    adata, 
    dir_deg, 
    dir_hvg,
    labels_key: str,
    deg_group: str | None = None,
    top_deg_num: int | None = 500,
    use_sampling: bool | None = False,
    sampling_ratio: float| None = 1.0,
    use_reference: bool | None = False, 
    top_hvg_genes: int = 2000,
    subset: bool = True,
    layer: str = "counts",
    batch_key: str | None = None,
    span: float = 0.5
):
    if use_reference:

        ref_nuclei = list(adata[adata.obs[labels_key] != 'Unknown'].obs.index)
        ref_adata = adata[ref_nuclei,:]

        _DEGs(
            adata = ref_adata, 
            dir_deg = dir_deg, 
            deg_group = deg_group,
            batch_key = batch_key,
            top_deg_num = top_deg_num,
            use_sampling = use_sampling,
            sampling_ratio = sampling_ratio
        )
        _HVGs(
            adata = ref_adata,
            dir_hvg = dir_hvg,
            top_hvg_genes = top_hvg_genes,
            subset = subset,
            layer = layer,
            batch_key = batch_key,
            span = span
        )
    else:
        _HVGs(
            adata = adata,
            dir_hvg = dir_hvg,
            top_hvg_genes = top_hvg_genes,
            subset = subset,
            layer = layer,
            batch_key = batch_key,
            span = span
        )
    
    adata = _subset_adata_DEGs_HVGs(adata, dir_hvg, dir_deg)
    return adata



def correction_validity(
    correction: str | None, 
    batch_key: str | None, 
    continuous_covariate_keys: list[str] | None, 
    categorical_covariate_keys: list[str] | None
):

    if isinstance(correction, str):
        assert correction.lower() in ['no_batch_covariates', 'batch_only', 'covariates_only', 'batch_covariates'] , 'correction can be no_batch_covariates, batch_only, covariates_only, batch_covariates.'
    if correction == 'no_batch_covariates' or correction == None:
        assert batch_key is None and continuous_covariate_keys is None and categorical_covariate_keys is None
    elif correction ==  'batch_only':
        assert batch_key is not None and continuous_covariate_keys is None and categorical_covariate_keys is None
    elif correction ==  'covariates_only':
        assert batch_key is None and (continuous_covariate_keys is not None or categorical_covariate_keys is not None)
    else:
        assert batch_key is not None and (continuous_covariate_keys is not None or categorical_covariate_keys is not None)
    


def scIntegration(
    adata,
    nums_random_trials: int = 5,
    **kwargs
):

    required_hyperparameter_tune = ["enable_parameter_tuning", "n_latent", "n_hidden", "n_layers", "dispersion", "lr", "max_epochs"]
    assert all(k in kwargs for k in required_hyperparameter_tune), f"Missing required kwargs: {required_hyperparameter_tune}"

    required_integration = ["use_unsupervised", "dir_scvi"]
    assert all(k in kwargs for k in required_integration), f"Missing required kwargs: {required_integration}"
    if not kwargs.get("use_unsupervised"):
        required_integration = ["dir_scanvi"]
        assert all(k in kwargs for k in required_integration), f"Missing required kwargs: {required_integration}"

    if kwargs.get("enable_parameter_tuning"):

        logger.info("Enable Hyperparameter Tuning ...")

        required_train_val_test_split = ["seed_value", "adjust_cols"]
        assert all(k in kwargs for k in required_train_val_test_split), f"Missing required kwargs: {required_train_val_test_split}"

        split_params = {}
        split_params['adata'] = adata
        split_params['train_ratio'] = kwargs.get('train_ratio', 0.6)
        split_params['test_ratio'] = kwargs.get('test_ratio', 0.5)
        split_params['seed_value'] = kwargs.get('seed_value')
        split_params['shuffle'] = kwargs.get('shuffle', True)
        split_params['adjust_cols'] = kwargs.get("adjust_cols")

        adata_train, adata_val, adata_test =  _train_valid_test_split(**split_params)

        tune_params = {}
        tune_params['adata_train'] = adata_train
        tune_params['adata_val'] = adata_val
        tune_params['accelerator'] = kwargs.get("accelerator", 'cpu')
        tune_params['labels_key'] = kwargs.get("labels_key", None)
        tune_params['layer'] = kwargs.get("layer", None)
        tune_params['batch_key'] = kwargs.get("batch_key", None)
        tune_params['categorical_covariate_keys'] = kwargs.get("categorical_covariate_keys", None)
        tune_params['continuous_covariate_keys'] = kwargs.get("continuous_covariate_keys", None)
        tune_params['n_latent'] = kwargs.get("n_latent")
        tune_params['n_hidden'] = kwargs.get("n_hidden")
        tune_params['n_layers'] = kwargs.get("n_layers")
        tune_params['dispersion'] = kwargs.get("dispersion")
        tune_params['lr'] = kwargs.get("lr")
        tune_params['max_epochs'] = kwargs.get("max_epochs")

        study = optuna.create_study(
            direction="maximize",
            sampler=optuna.samplers.TPESampler(seed=kwargs.get('seed_value'))
        )
        study.optimize(
            lambda trial: _hyperparameter_tune(trial = trial, **tune_params),
            n_trials=nums_random_trials,
            show_progress_bar=False
        )
        logger.info("Using Best Hyperparameter Values for Data Integration ...")
        best_params = study.best_params
        for k, v in best_params.items():
            logger.info("Best Hyperparameters after Tuning: {} = {}".format(k, v))
    else:
        logger.info("Using Default Hyperparameter Values for Data Integration ...")
        best_params = dict()
        best_params["n_latent"] = kwargs.get("n_latent")
        best_params["n_hidden"] = kwargs.get("n_hidden")
        best_params["n_layers"] = kwargs.get("n_layers")
        best_params["dispersion"] = kwargs.get("dispersion")
        best_params["lr"] = kwargs.get("lr")
        best_params["max_epochs"] = kwargs.get("max_epochs")

        for k, v in best_params.items():
            logger.info("Default Hyperparameter Value: {} = {}".format(k, v))


    integration_params = {}
    integration_params['adata'] = adata
    integration_params['level'] = kwargs.get("level")
    integration_params['dir_scvi'] = kwargs.get("dir_scvi")
    integration_params['dir_scanvi'] = kwargs.get("dir_scanvi", None)
    integration_params['best_params'] = best_params
    integration_params['batch_key'] = kwargs.get("batch_key", None)
    integration_params['labels_key'] = kwargs.get("labels_key", None)
    integration_params['prediction_key'] = kwargs.get("prediction_key", None)
    integration_params['layer'] = kwargs.get("layer", None)
    integration_params['categorical_covariate_keys'] = kwargs.get("categorical_covariate_keys", None)
    integration_params['continuous_covariate_keys'] = kwargs.get("continuous_covariate_keys", None)
    integration_params['accelerator'] = kwargs.get("accelerator", 'cpu')
    integration_params['scanvi_epochs'] = kwargs.get("scanvi_epochs", 'cpu')
    integration_params['use_unsupervised'] = kwargs.get("use_unsupervised")
    
    adata = _scIntegration(**integration_params)

    return adata



def UMAP_visualization(
    adata,
    seed_value,
    level,
    use_unsupervised
):
    umap_coord = umap.UMAP(random_state = seed_value)

    if use_unsupervised:
        X_umap = umap_coord.fit_transform(adata.obsm[f'Level{level}_X_scVI'])
    else:
        X_umap = umap_coord.fit_transform(adata.obsm[f'Level{level}_X_scANVI'])

    adata.obsm[f'Level{level}_X_umap'] = X_umap

    return adata


def nuclei_clustering(
    adata,
    level,
    method,
    resolution,
    seed_value,
    use_unsupervised
):
    if use_unsupervised:
        sc.pp.neighbors(adata, use_rep = f'Level{level}_X_scVI')
    else:
        sc.pp.neighbors(adata, use_rep = f'Level{level}_X_scANVI')
    if method.lower() == 'leiden':
        sc.tl.leiden(adata, resolution = resolution, random_state = seed_value)
    elif method.lower() == 'louvain':
        sc.tl.louvain(adata, resolution = resolution, random_state = seed_value)
    else:
        raise ValueError(f"Provided clustering method {method} is not supported!")
    return adata




def load_protein_coding_genes(
    dir_file, 
    header = True
):
    protein_coding_genes = []
    with open(dir_file, mode='r') as f:
        if header:
            next(f)
        for line in f:
            ensid, is_pc = line.strip("\n").split("\t")
            if ensid != '' and is_pc == 'protein-coding gene':
                protein_coding_genes.append(ensid)
    return protein_coding_genes




def add_standard_meta(
    dir_meta: str, 
    adata: AnnData,
    dataset_key: str,
    sampleid_key: str,
    meta_cols: list[str]
):
    for ds in set(adata.obs[dataset_key]):
        curr_metadata = pd.read_csv(os.path.join(dir_meta, ds, 'metadata', 'sample_standard.tsv'), sep = "\t")
        curr_metadata = curr_metadata.fillna('Unknown')
        for meta in meta_cols:
            curr_samples2meta = _meta2dict(df = curr_metadata, prop = meta)
            for sid in set(curr_metadata[sampleid_key]):
                sel_nuclei = (adata.obs[sampleid_key] == str(sid)) & (adata.obs[dataset_key] == ds)
                if sum(sel_nuclei) > 0:
                    logger.info('Adding {} Information - Sample {} from Dataset {}'.format(meta, sid, ds))
                    adata.obs.loc[sel_nuclei, meta] = str(curr_samples2meta[sid])
    return adata



def merge_datasets(
    all_datasets: list[str],
    dir_datasets: str,
    ref_label: str,
    cell_type_key: str,
    dataset_key: str,
    sampleid_key: str,
    batch_key: str,
    kept_genes: str,
    protein_coding_genes: set | None,
    add_meta: bool,
    meta_cols : list[str] | str | None = None
):
    integrate_datasets = []

    assert len(all_datasets) > 0, 'No dataset ready for integration.'
    if kept_genes == GENE_SELECT_PROTEIN_CODING:
        assert protein_coding_genes is not None, 'protein_coding_genes can not be None, Missing protein coding genes information.'

    for ds in all_datasets:
        logger.info(f"========= {dataset_key}:{ds} =========")
        curr_adata = ad.read_h5ad(os.path.join(dir_datasets, ds, 'processed/export/adata.h5ad'))
        if not curr_adata.var_names.str.startswith('ENSG').all(): # might has transcriptomic id
            logger.info("Please check dataset {}, its gene IDs are not ensemble ID.".format(ds))
            break
        if cell_type_key in curr_adata.obs.columns:
            curr_adata.obs[ref_label] = curr_adata.obs[cell_type_key]
        else:
            curr_adata.obs[ref_label] = "Unknown"
        if ds in REFERENCE_CT2SUPERCT:
            if 'nonregex' in REFERENCE_CT2SUPERCT[ds]:
                curr_adata.obs[ref_label].replace(REFERENCE_CT2SUPERCT[ds]['nonregex'], inplace = True)
            for rex in REFERENCE_CT2SUPERCT[ds].get('regex', []):
                curr_adata.obs[ref_label] = curr_adata.obs[ref_label].str.replace(rex, REFERENCE_CT2SUPERCT[ds]['regex'][rex], regex = True, flags=re.IGNORECASE)
            if set(curr_adata.obs[ref_label]) - SUPERCTs:
                logger.info("There is something wrong with super cell type mapping in dataset {}!".format(ds))
                break

        curr_adata.obs[dataset_key] = ds
        if batch_key not in curr_adata.obs.columns:
            curr_adata.obs[batch_key] = curr_adata.obs[sampleid_key].astype(str)
        curr_adata.obs[batch_key] = ds + '_' + curr_adata.obs[batch_key].astype(str)
        curr_adata.var.drop(curr_adata.var.columns, axis=1, inplace=True)
        curr_adata.X = curr_adata.X.astype(np.float32)
        float64_cols = curr_adata.obs.select_dtypes(include=['float64']).columns
        if len(float64_cols) > 0:
            curr_adata.obs[float64_cols] = curr_adata.obs[float64_cols].astype('float32')
        int64_cols = curr_adata.obs.select_dtypes(include=['int64']).columns
        if len(int64_cols) > 0:
            curr_adata.obs[int64_cols] = curr_adata.obs[int64_cols].astype('int32')
        if kept_genes == GENE_SELECT_PROTEIN_CODING:
            logger.info("before filter with protein_coding genes, curr_adata.shape = (#nuclei:{:,}, #genes:{:,})".format(curr_adata.shape[0], curr_adata.shape[1]))
            pc_genes = set(curr_adata.var_names) & protein_coding_genes
            curr_adata = curr_adata[:,list(pc_genes)]
            logger.info("after filter with protein_coding genes, curr_adata.shape = (#nuclei:{:,}, #genes:{:,}) ".format(curr_adata.shape[0], curr_adata.shape[1]))
        integrate_datasets.append(curr_adata)
    

    if kept_genes == GENE_SELECT_INTERSECTION:
        adata = ad.concat(integrate_datasets, merge = "unique")
    elif kept_genes == GENE_SELECT_UNION or kept_genes == GENE_SELECT_PROTEIN_CODING:
        adata = ad.concat(integrate_datasets, merge = "unique", join = "outer")
    else:
        raise ValueError('kept_genes could be: intersection, union, protein_coding.')


    float64_cols = adata.obs.select_dtypes(include=['float64']).columns
    if len(float64_cols) > 0:
        adata.obs[float64_cols] = adata.obs[float64_cols].astype('float32')

    int64_cols = adata.obs.select_dtypes(include=['int64']).columns
    if len(int64_cols) > 0:
        adata.obs[int64_cols] = adata.obs[int64_cols].astype('int32')


    if adata.obsm is not None:
        del adata.obsm


    if add_meta:
        assert meta_cols is not None, 'meta_cols can not be None, missing meta information.'
        if isinstance(meta_cols, str):
            meta_cols = [meta_cols]
        adata = add_standard_meta(dir_meta = dir_datasets, 
                                adata = adata, 
                                dataset_key = dataset_key, 
                                sampleid_key = sampleid_key, 
                                meta_cols = meta_cols
                            )

    logger.info("Integrated data contains {:,} nuclei and {:,} genes.".format(adata.X.shape[0], adata.X.shape[1]))
    
    return adata





def evaluation(
    dir_adata: str,
    ref_label: str,
    prediction_key: str,
    level: int,
    meta_cols: list[str],
    dir_evaluation: str | None = None
):
    logger.info(f'Evaluate Level{level} Cell Type Annotation ... ')
    adata = ad.read_h5ad(dir_adata)
    if str(level) == '1':
        if 'regex' in NONREFERENCE_CT2SUPERCT:
            for rex in NONREFERENCE_CT2SUPERCT['regex']:
                adata.obs[ref_label] = adata.obs[ref_label].str.replace(rex, NONREFERENCE_CT2SUPERCT['regex'][rex], regex = True, flags=re.IGNORECASE)
        if 'nonregex' in NONREFERENCE_CT2SUPERCT:
            adata.obs[ref_label].replace(NONREFERENCE_CT2SUPERCT['nonregex'], inplace = True)

    logger.info(f'Evaluate Cell Type Annotation on Datasets with Original Cell Type Annotations ... ')
    if '' in set(adata.obs[ref_label]):
        adata = adata[adata.obs[ref_label] != '']
    elif 'Unknown' in set(adata.obs[ref_label]):
        adata = adata[adata.obs[ref_label] != 'Unknown']


    for meta in meta_cols:
        for ma in set(adata.obs[meta]):
            logger.info(f"========= {meta}:{ma} =========")
            curr_adata = adata[adata.obs[meta] == ma]
            curr_label = list(set(curr_adata.obs[ref_label]))
            curr_label.sort()
            curr_ps = precision_score(list(curr_adata.obs[ref_label]), list(curr_adata.obs[f'Level{level}_{prediction_key}']), labels=curr_label, average=None)
            curr_rs = recall_score(list(curr_adata.obs[ref_label]), list(curr_adata.obs[f'Level{level}_{prediction_key}']), labels=curr_label, average=None)
            curr_eval = pd.DataFrame({
                "CellType": curr_label,
                "PrecisionScore": list(curr_ps),
                "RecallScore": list(curr_rs)
            })
            if dir_evaluation is not None:
                curr_eval.to_csv(os.path.join(dir_evaluation, f"{meta}_{ma}_evaluation.txt"), sep = '\t', index = None)
            logger.info(f"{meta}:{ma} Precision & Recall = {curr_eval}.")
            curr_cf = pd.crosstab(curr_adata.obs[ref_label], curr_adata.obs[f'Level{level}_{prediction_key}'])
            logger.info(f"{meta}:{ma} Confusion Matrix: {curr_cf}.")





def save_anndata(
    dir_adata, 
    adata,
    level,
    prediction_key,
    use_unsupervised,
    dir_save,
):
    adata_preintegration = ad.read_h5ad(dir_adata)
    assert adata_preintegration.shape[0] == adata.shape[0], 'missing nuclei/cells during integration.'

    adata_preintegration.obsm[f'Level{level}_X_umap'] = adata.obsm.pop(f'Level{level}_X_umap')

    if not use_unsupervised:
        assert level is not None, 'provide level.'
        adata_preintegration.obsm[f'Level{level}_X_scANVI'] = adata_preintegration.obsm.pop(f'Level{level}_X_scANVI')
    else:
        adata_preintegration.obsm[f'Level{level}_X_scVI'] = adata_preintegration.obsm.pop(f'Level{level}_X_scVI')


    if (prediction_key is not None) and (f'Level{level}_{prediction_key}' in adata.obs.columns):
        adata_preintegration.obs.loc[adata.obs.index, f'Level{level}_{prediction_key}'] = adata.obs.loc[adata.obs.index, f'Level{level}_{prediction_key}']

    adata_preintegration.write_h5ad(os.path.join(dir_save, 'raw_count_after_integration.h5ad'))




