import os
import anndata as ad
import scanpy as sc
import matplotlib.pyplot as plt
import pandas as pd
import re
import numpy as np


def load_protein_coding_genes(dir_file, header = True):
    protein_coding_genes = []
    with open(dir_file, mode='r') as f:
        if header:
            next(f)
        for line in f:
            ensid, is_pc = line.strip("\n").split("\t")
            if ensid != '' and is_pc == 'protein-coding gene':
                protein_coding_genes.append(ensid)
    return protein_coding_genes


GENE_SELECT_INTERSECTION = "INTERSECTION"
GENE_SELECT_UNION = "UNION"
GENE_SELECT_PROTEIN_CODING = "PROTEIN_CODING"


#add normalized celltype information here

REFERENCE_CT2SUPERCT = {
    '37824663':{
        'nonregex':{
                    'oligodendrocyte': 'Oligodendrocyte', 
                    'neuron': 'Neuron', 
                    'astrocyte': 'Astrocyte', 
                    'oligodendrocyte precursor cell': 'OPC', 
                    'central nervous system macrophage': 'Immune', 
                    'fibroblast': 'Others', 
                    'Bergmann glial cell': 'Others', 
                    'choroid plexus epithelial cell': 'Others', 
                    'ependymal cell': 'Others', 
                    'endothelial cell':'Others', 
                    'pericyte':'Others', 
                    'leukocyte':'Immune', 
                    'vascular associated smooth muscle cell': 'Others'
        }          
    },
    '37824655':{
        'nonregex':{
                'oligodendrocyte': 'Oligodendrocyte', 
                'microglial cell': 'Immune', 
                'astrocyte of the cerebral cortex': 'Astrocyte', 
                'oligodendrocyte precursor cell': 'OPC', 
                'vascular leptomeningeal cell': 'Others', 
                'cerebral cortex endothelial cell': 'Others',
        },
        'regex':{
                r'.*neuron.*':'Neuron'
                # r'.*interneuron.*':'Neuron',
        }
    }
}


SUPERCTs = {'Astrocyte', 'Immune', 'Neuron', 'Oligodendrocyte', 'OPC', 'Others'}


KEPT_GENES = GENE_SELECT_INTERSECTION

# assert KEPT_GENES in ['intersection', 'union', 'protein_coding'], 'KEPT_GENES could be: intersection, union, protein_coding.'

dir_data = '/home/xuj2/isilon/Cheng-Qiu/TACA2/Data'
dir_integration = '/home/xuj2/isilon/Cheng-Qiu/TACA2/Integration_test'

all_ds = [name for name in os.listdir(dir_data) if os.path.isdir(os.path.join(dir_data, name))]
all_ds.remove('Template')

print('There are {:,} datasets used for integration.'.format(len(all_ds)))

protein_coding_genes = set(load_protein_coding_genes(dir_file = os.path.join(dir_integration, 'data', 'protein_coding_genes.tsv')))
print('There are totally {:,} protein coding genes.'.format(len(protein_coding_genes)))


integrate_datasets = []

for ds in all_ds:
    print("========= Dataset {} =========".format(ds))
    curr_adata = ad.read_h5ad(os.path.join(dir_data, ds, 'processed/export/adata.h5ad'))
    if not curr_adata.var_names.str.startswith('ENSG').all(): # might has transcriptomic id
        print("Please check dataset {}, its gene IDs are not ensemble ID.".format(ds))
        break
    if 'cell_type' in curr_adata.obs.columns:
        curr_adata.obs['Super_Celltype'] = curr_adata.obs['cell_type']
    else:
        curr_adata.obs['Super_Celltype'] = ""
    if ds in REFERENCE_CT2SUPERCT:
        if 'nonregex' in REFERENCE_CT2SUPERCT[ds]:
            curr_adata.obs['Super_Celltype'].replace(REFERENCE_CT2SUPERCT[ds]['nonregex'], inplace = True)
        for rex in REFERENCE_CT2SUPERCT[ds].get('regex', []):
            curr_adata.obs['Super_Celltype'] = curr_adata.obs['Super_Celltype'].str.replace(rex, REFERENCE_CT2SUPERCT[ds]['regex'][rex], regex = True, flags=re.IGNORECASE)
        if set(curr_adata.obs['Super_Celltype']) - SUPERCTs:
            print("There is something wrong with super cell type mapping in dataset {}!".format(ds))
            break

    curr_adata.obs['Dataset'] = ds
    if 'BatchID' not in curr_adata.obs.columns:
        curr_adata.obs['BatchID'] = curr_adata.obs['SampleID'].astype(str)
    curr_adata.obs['BatchID'] = ds + '_' + curr_adata.obs['BatchID'].astype(str)
    curr_adata.var.drop(curr_adata.var.columns, axis=1, inplace=True)
    curr_adata.X = curr_adata.X.astype(np.float32)
    float64_cols = curr_adata.obs.select_dtypes(include=['float64']).columns
    if len(float64_cols) > 0:
        curr_adata.obs[float64_cols] = curr_adata.obs[float64_cols].astype('float32')
    int64_cols = curr_adata.obs.select_dtypes(include=['int64']).columns
    if len(int64_cols) > 0:
        curr_adata.obs[int64_cols] = curr_adata.obs[int64_cols].astype('int32')
    if KEPT_GENES == GENE_SELECT_PROTEIN_CODING: # add check geneID type
        print("before filter with protein_coding genes, curr_adata.shape = (#nuclei:{:,}, #genes:{:,})".format(curr_adata.shape[0], curr_adata.shape[1]))
        pc_genes = set(curr_adata.var_names) & protein_coding_genes
        curr_adata = curr_adata[:,list(pc_genes)]
        print("after filter with protein_coding genes, curr_adata.shape = (#nuclei:{:,}, #genes:{:,}) ".format(curr_adata.shape[0], curr_adata.shape[1]))
    integrate_datasets.append(curr_adata)



if KEPT_GENES == GENE_SELECT_INTERSECTION:
    adata = ad.concat(integrate_datasets, merge = "unique")
elif KEPT_GENES == GENE_SELECT_UNION or KEPT_GENES == GENE_SELECT_PROTEIN_CODING:
    adata = ad.concat(integrate_datasets, merge = "unique", join = "outer")
else:
    raise ValueError('KEPT_GENES could be: intersection, union, protein_coding.')


float64_cols = adata.obs.select_dtypes(include=['float64']).columns
if len(float64_cols) > 0:
    adata.obs[float64_cols] = adata.obs[float64_cols].astype('float32')

int64_cols = adata.obs.select_dtypes(include=['int64']).columns
if len(int64_cols) > 0:
    adata.obs[int64_cols] = adata.obs[int64_cols].astype('int32')


print("Integrated data contains {:,} nuclei and {:,} genes.".format(adata.X.shape[0], adata.X.shape[1]))

adata.write_h5ad(os.path.join('../output_TACA2/', 'preintegration_{}_reorganized_32bits_{}datasets_BatchID.h5ad'.format(KEPT_GENES, len(all_ds))))



