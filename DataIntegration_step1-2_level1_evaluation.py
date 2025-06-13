import anndata as ad
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import scvi
from sklearn.metrics import confusion_matrix, precision_score, recall_score
import os
import re


pd.set_option('display.max_columns', None)  # None means unlimited
pd.set_option('display.max_rows', None)     # None means unlimited

adata = ad.read_h5ad("../output_TACA2/integration_with_prediction_37824663_37824655_max_epoch_200_25_21_datasets.h5ad")


NONREFERENCE_CT2SUPERCT = {
    'regex':{
            r'^Exc.*': 'Neuron',
            r'^Inh.*': 'Neuron',
            r'^Ast.*': 'Astrocyte',
            r'^Mic.*': 'Immune',
            r'^Fib.*': 'Others',
            # r'\b(neuron|inhibitory interneuron)\b':'Neuron',
            r'.*neuron.*':'Neuron',
    },
    'nonregex':{
                'Oli': 'Oligodendrocyte',
                'End': 'Others',
                'Per': 'Others',
                'CAMs': 'Immune',
                'T cells': 'Immune',
                'L2/3 IT':'Neuron',
                'Pvalb':'Neuron',
                'L5 IT':'Neuron',
                'Vip':'Neuron',
                'L4 IT':'Neuron',
                'Sst':'Neuron',
                'Lamp5':'Neuron',
                'L6 IT':'Neuron',
                'L6 CT':'Neuron',
                'Sncg':'Neuron',
                'Lamp5 Lhx6':'Neuron',
                'L5/6 NP':'Neuron',
                'L6b':'Neuron',
                'Chandelier':'Neuron',
                'L6 IT Car3':'Neuron',
                'Pax6':'Neuron',
                'L5 ET': 'Neuron',
                'Sst Chodl': 'Neuron',
                'Microglia-PVM':'Immune',
                'Oligodendrocytes': 'Oligodendrocyte',
                'Non-DA Neurons':'Neuron',
                'Astrocytes': 'Astrocyte',
                'DA Neurons':'Neuron',
                'Endothelial cells': 'Others',
                'Microglia': 'Immune',
                'oligodendrocyte': 'Oligodendrocyte',
                'microglial cell': 'Immune',
                'oligodendrocyte precursor cell':'OPC',
                'astrocyte': 'Astrocyte',
                'cerebral cortex neuron':'Neuron',
                'interneuron': 'Neuron',
                'pyramidal neuron': 'Neuron',
                'endothelial cell': 'Others',
                'pericyte': 'Others',
                'macrophage': 'Immune', 
                'smooth muscle cell': 'Others',
                'vascular leptomeningeal cell': 'Others', 
                'myeloid cell': 'Immune',
                'erythroid lineage cell': 'Others',
                'ODC': 'Oligodendrocyte',
                'Astro': 'Astrocyte',
                'MG': 'Immune',
                'Endo':'Others',
                'cerebellar granule cell':'Neuron',
                'mural cell': 'Others',
                'leukocyte': 'Immune',
                'differentiation-committed oligodendrocyte precursor': 'Oligodendrocyte',
                'vascular associated smooth muscle cell': 'Others',
                'central nervous system macrophage': 'Immune',
                'neuron': 'Neuron', 
                'progenitor cell': 'OPC',
                'VLMC': 'Others',
                'Endothelial': 'Others',
                'capillary endothelial cell': 'Others',
                'SMC': 'Others',
                'endothelial cell of artery': 'Others',
                'T cell': 'Immune',
                'T_Cell': 'Immune',
                'Mural': 'Others',
                'Oligo': 'Oligodendrocyte',
                'L2_L3': 'Neuron',
                'L3_L5': 'Neuron',
                'L4_L6': 'Neuron',
                'L4_L5': 'Neuron',
                'L5_L6': 'Neuron',
                'L6': 'Neuron',
                'L5': 'Neuron',
                'OPCs': 'OPC',
                'CPEC': 'Others',
                'Epd': 'Others',
                'CUX2+': 'Neuron',
                'Rosehip': 'Neuron',
                'SOM': 'Neuron',
                'PV': 'Neuron',
                '5HT3aR': 'Neuron',
                'B cell': 'Immune'
    }
}

# test order (TODO)

if 'regex' in NONREFERENCE_CT2SUPERCT:
    for rex in NONREFERENCE_CT2SUPERCT['regex']:
        adata.obs['Super_Celltype'] = adata.obs['Super_Celltype'].str.replace(rex, NONREFERENCE_CT2SUPERCT['regex'][rex], regex = True, flags=re.IGNORECASE)

if 'nonregex' in NONREFERENCE_CT2SUPERCT:
    adata.obs['Super_Celltype'].replace(NONREFERENCE_CT2SUPERCT['nonregex'], inplace = True)



dir_TACA = '/home/xuj2/isilon/Cheng-Qiu/TACA'
dir_data = '/home/xuj2/isilon/Cheng-Qiu/TACA2/Data'
all_ds = [name for name in os.listdir(dir_data) if os.path.isdir(os.path.join(dir_data, name))]
all_ds.remove('Template')

nonreference_ds = [ds for ds in all_ds]
nonreference_ds.remove('37824663')
nonreference_ds.remove('37824655')



for ds in nonreference_ds:
    print("========= Dataset {} =========".format(ds))
    curr_adata = adata[adata.obs['Dataset'] == ds]
    original_no_superct = (curr_adata.obs['Super_Celltype'] == '').all()
    if not original_no_superct:
        curr_label = list(set(curr_adata.obs['Super_Celltype']))
        curr_label.sort()
        curr_ps = precision_score(list(curr_adata.obs['Super_Celltype']), list(curr_adata.obs['C_scANVI']), labels=curr_label, average=None)
        curr_rs = recall_score(list(curr_adata.obs['Super_Celltype']), list(curr_adata.obs['C_scANVI']), labels=curr_label, average=None)
        curr_eval = pd.DataFrame({
            "CellType": curr_label,
            "PrecisionScore": list(curr_ps),
            "RecallScore": list(curr_rs)
        })
        print("curr_precision_recall = ", curr_eval)
        curr_compare = pd.crosstab(curr_adata.obs['Super_Celltype'], curr_adata.obs['C_scANVI'])
        print("curr_compare = ", curr_compare)
    else:
        print("{} has no original label.".format(ds))




















