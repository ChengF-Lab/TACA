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
pd.set_option('display.max_columns', None)  # None means unlimited
pd.set_option('display.max_rows', None)     # None means unlimited

pre_adata = ad.read_h5ad("../output_TACA2/neuron_final/preintegration_PROTEIN_CODING_reorganized_32bits_21datasets_neuron.h5ad")
adata = ad.read_h5ad("../output_TACA2/neuron_final/integration_PROTEIN_CODING_reorganized_32bits_21datasets_neuron_batch_min10_HVG2000_epoch100_latent50_lr5e-4_all_nuclei_res20_UMAP.h5ad")

assert adata.shape[0] == pre_adata.shape[0], 'Neuron nuclei # mismatch.'

pre_adata.obs['louvain'] = ''
pre_adata.obs.loc[adata.obs.index, 'louvain'] = adata.obs.loc[adata.obs.index, 'louvain']


pre_adata.obs['neuron_L2'] = 'Unknown'


pre_adata.obs.loc[pre_adata.obs['louvain'] == '0', 'neuron_L2'] = 'GABAergic inhibitory'
pre_adata.obs.loc[pre_adata.obs['louvain'] == '1', 'neuron_L2'] = 'Glutamatergic excitatory'
pre_adata.obs.loc[pre_adata.obs['louvain'] == '2', 'neuron_L2'] = 'Glutamatergic excitatory'
pre_adata.obs.loc[pre_adata.obs['louvain'] == '3', 'neuron_L2'] = 'Glutamatergic excitatory'
pre_adata.obs.loc[pre_adata.obs['louvain'] == '4', 'neuron_L2'] = 'Glutamatergic excitatory'
pre_adata.obs.loc[pre_adata.obs['louvain'] == '5', 'neuron_L2'] = 'Glutamatergic excitatory'
pre_adata.obs.loc[pre_adata.obs['louvain'] == '6', 'neuron_L2'] = 'Glutamatergic excitatory'
pre_adata.obs.loc[pre_adata.obs['louvain'] == '7', 'neuron_L2'] = 'Splatter'
pre_adata.obs.loc[pre_adata.obs['louvain'] == '8', 'neuron_L2'] = 'GABAergic inhibitory'
pre_adata.obs.loc[pre_adata.obs['louvain'] == '9', 'neuron_L2'] = 'Glutamatergic excitatory'

pre_adata.obs.loc[pre_adata.obs['louvain'] == '10', 'neuron_L2'] = 'Glutamatergic excitatory'
pre_adata.obs.loc[pre_adata.obs['louvain'] == '11', 'neuron_L2'] = 'GABAergic inhibitory'
pre_adata.obs.loc[pre_adata.obs['louvain'] == '12', 'neuron_L2'] = 'Glutamatergic excitatory'
pre_adata.obs.loc[pre_adata.obs['louvain'] == '13', 'neuron_L2'] = 'Glutamatergic excitatory'
pre_adata.obs.loc[pre_adata.obs['louvain'] == '14', 'neuron_L2'] = 'Glutamatergic excitatory'
pre_adata.obs.loc[pre_adata.obs['louvain'] == '15', 'neuron_L2'] = 'GABAergic inhibitory'
pre_adata.obs.loc[pre_adata.obs['louvain'] == '16', 'neuron_L2'] = 'GABAergic inhibitory'
pre_adata.obs.loc[pre_adata.obs['louvain'] == '17', 'neuron_L2'] = 'Glutamatergic excitatory'
pre_adata.obs.loc[pre_adata.obs['louvain'] == '18', 'neuron_L2'] = 'Glutamatergic excitatory'
pre_adata.obs.loc[pre_adata.obs['louvain'] == '19', 'neuron_L2'] = 'Glutamatergic excitatory'

pre_adata.obs.loc[pre_adata.obs['louvain'] == '20', 'neuron_L2'] = 'GABAergic inhibitory'
pre_adata.obs.loc[pre_adata.obs['louvain'] == '21', 'neuron_L2'] = 'Glutamatergic excitatory'
pre_adata.obs.loc[pre_adata.obs['louvain'] == '22', 'neuron_L2'] = 'Glutamatergic excitatory'
pre_adata.obs.loc[pre_adata.obs['louvain'] == '23', 'neuron_L2'] = 'Glutamatergic excitatory'
pre_adata.obs.loc[pre_adata.obs['louvain'] == '24', 'neuron_L2'] = 'Medium spiny neuron'
pre_adata.obs.loc[pre_adata.obs['louvain'] == '25', 'neuron_L2'] = 'Upper rhombic lip'
pre_adata.obs.loc[pre_adata.obs['louvain'] == '26', 'neuron_L2'] = 'Glutamatergic excitatory'
pre_adata.obs.loc[pre_adata.obs['louvain'] == '27', 'neuron_L2'] = 'GABAergic inhibitory'
pre_adata.obs.loc[pre_adata.obs['louvain'] == '28', 'neuron_L2'] = 'Glutamatergic excitatory'
pre_adata.obs.loc[pre_adata.obs['louvain'] == '29', 'neuron_L2'] = 'Midbrain-derived inhibitory'

pre_adata.obs.loc[pre_adata.obs['louvain'] == '30', 'neuron_L2'] = 'Glutamatergic excitatory'
pre_adata.obs.loc[pre_adata.obs['louvain'] == '31', 'neuron_L2'] = 'GABAergic inhibitory'
pre_adata.obs.loc[pre_adata.obs['louvain'] == '32', 'neuron_L2'] = 'GABAergic inhibitory'
pre_adata.obs.loc[pre_adata.obs['louvain'] == '33', 'neuron_L2'] = 'Glutamatergic excitatory'
pre_adata.obs.loc[pre_adata.obs['louvain'] == '34', 'neuron_L2'] = 'Glutamatergic excitatory'
pre_adata.obs.loc[pre_adata.obs['louvain'] == '35', 'neuron_L2'] = 'Glutamatergic excitatory'
pre_adata.obs.loc[pre_adata.obs['louvain'] == '36', 'neuron_L2'] = 'Glutamatergic excitatory'
pre_adata.obs.loc[pre_adata.obs['louvain'] == '37', 'neuron_L2'] = 'GABAergic inhibitory'
pre_adata.obs.loc[pre_adata.obs['louvain'] == '38', 'neuron_L2'] = 'Glutamatergic excitatory'
pre_adata.obs.loc[pre_adata.obs['louvain'] == '39', 'neuron_L2'] = 'Glutamatergic excitatory'

pre_adata.obs.loc[pre_adata.obs['louvain'] == '40', 'neuron_L2'] = 'GABAergic inhibitory'
pre_adata.obs.loc[pre_adata.obs['louvain'] == '41', 'neuron_L2'] = 'Glutamatergic excitatory'
pre_adata.obs.loc[pre_adata.obs['louvain'] == '42', 'neuron_L2'] = 'GABAergic inhibitory'
pre_adata.obs.loc[pre_adata.obs['louvain'] == '43', 'neuron_L2'] = 'Glutamatergic excitatory'
pre_adata.obs.loc[pre_adata.obs['louvain'] == '44', 'neuron_L2'] = 'Glutamatergic excitatory'
pre_adata.obs.loc[pre_adata.obs['louvain'] == '45', 'neuron_L2'] = 'GABAergic inhibitory'
pre_adata.obs.loc[pre_adata.obs['louvain'] == '46', 'neuron_L2'] = 'GABAergic inhibitory'
pre_adata.obs.loc[pre_adata.obs['louvain'] == '47', 'neuron_L2'] = 'Glutamatergic excitatory'
pre_adata.obs.loc[pre_adata.obs['louvain'] == '48', 'neuron_L2'] = 'Amygdala excitatory'
pre_adata.obs.loc[pre_adata.obs['louvain'] == '49', 'neuron_L2'] = 'Hippocampal CA1-3'

pre_adata.obs.loc[pre_adata.obs['louvain'] == '50', 'neuron_L2'] = 'GABAergic inhibitory'
pre_adata.obs.loc[pre_adata.obs['louvain'] == '51', 'neuron_L2'] = 'Thalamic excitatory'
pre_adata.obs.loc[pre_adata.obs['louvain'] == '52', 'neuron_L2'] = 'Unknown neuron'
pre_adata.obs.loc[pre_adata.obs['louvain'] == '53', 'neuron_L2'] = 'Hippocampal dentate gyrus'
pre_adata.obs.loc[pre_adata.obs['louvain'] == '54', 'neuron_L2'] = 'GABAergic inhibitory'
pre_adata.obs.loc[pre_adata.obs['louvain'] == '55', 'neuron_L2'] = 'Dopaminergic neuron'
pre_adata.obs.loc[pre_adata.obs['louvain'] == '56', 'neuron_L2'] = 'Lower rhombic lip'
pre_adata.obs.loc[pre_adata.obs['louvain'] == '57', 'neuron_L2'] = 'GABAergic inhibitory'
pre_adata.obs.loc[pre_adata.obs['louvain'] == '58', 'neuron_L2'] = 'GABAergic inhibitory'
pre_adata.obs.loc[pre_adata.obs['louvain'] == '59', 'neuron_L2'] = 'Eccentric medium spiny neuron'

pre_adata.obs.loc[pre_adata.obs['louvain'] == '60', 'neuron_L2'] = 'Glutamatergic excitatory'
pre_adata.obs.loc[pre_adata.obs['louvain'] == '61', 'neuron_L2'] = 'Amygdala excitatory'
pre_adata.obs.loc[pre_adata.obs['louvain'] == '62', 'neuron_L2'] = 'Glutamatergic excitatory'
pre_adata.obs.loc[pre_adata.obs['louvain'] == '63', 'neuron_L2'] = 'Glutamatergic excitatory'
pre_adata.obs.loc[pre_adata.obs['louvain'] == '64', 'neuron_L2'] = 'Glutamatergic excitatory'
pre_adata.obs.loc[pre_adata.obs['louvain'] == '65', 'neuron_L2'] = 'Glutamatergic excitatory'
pre_adata.obs.loc[pre_adata.obs['louvain'] == '66', 'neuron_L2'] = 'Glutamatergic excitatory'
pre_adata.obs.loc[pre_adata.obs['louvain'] == '67', 'neuron_L2'] = 'Splatter'
pre_adata.obs.loc[pre_adata.obs['louvain'] == '68', 'neuron_L2'] = 'Thalamic excitatory'
pre_adata.obs.loc[pre_adata.obs['louvain'] == '69', 'neuron_L2'] = 'Glutamatergic excitatory'

pre_adata.obs.loc[pre_adata.obs['louvain'] == '70', 'neuron_L2'] = 'Cerebellar inhibitory'
pre_adata.obs.loc[pre_adata.obs['louvain'] == '71', 'neuron_L2'] = 'Mammillary body'
pre_adata.obs.loc[pre_adata.obs['louvain'] == '72', 'neuron_L2'] = 'Glutamatergic excitatory'
pre_adata.obs.loc[pre_adata.obs['louvain'] == '73', 'neuron_L2'] = 'Splatter'
pre_adata.obs.loc[pre_adata.obs['louvain'] == '74', 'neuron_L2'] = 'Unknown neuron'
pre_adata.obs.loc[pre_adata.obs['louvain'] == '75', 'neuron_L2'] = 'Glutamatergic excitatory'
pre_adata.obs.loc[pre_adata.obs['louvain'] == '76', 'neuron_L2'] = 'Glutamatergic excitatory'
pre_adata.obs.loc[pre_adata.obs['louvain'] == '77', 'neuron_L2'] = 'Glutamatergic excitatory'
pre_adata.obs.loc[pre_adata.obs['louvain'] == '78', 'neuron_L2'] = 'Splatter'
pre_adata.obs.loc[pre_adata.obs['louvain'] == '79', 'neuron_L2'] = 'Hippocampal CA4' 

pre_adata.obs.loc[pre_adata.obs['louvain'] == '80', 'neuron_L2'] = 'Unknown neuron'
pre_adata.obs.loc[pre_adata.obs['louvain'] == '81', 'neuron_L2'] = 'Glutamatergic excitatory'
pre_adata.obs.loc[pre_adata.obs['louvain'] == '82', 'neuron_L2'] = 'Unknown neuron'
pre_adata.obs.loc[pre_adata.obs['louvain'] == '83', 'neuron_L2'] = 'Amygdala excitatory'
pre_adata.obs.loc[pre_adata.obs['louvain'] == '84', 'neuron_L2'] = 'Splatter'
pre_adata.obs.loc[pre_adata.obs['louvain'] == '85', 'neuron_L2'] = 'GABAergic inhibitory'
pre_adata.obs.loc[pre_adata.obs['louvain'] == '86', 'neuron_L2'] = 'Hippocampal CA1-3'
pre_adata.obs.loc[pre_adata.obs['louvain'] == '87', 'neuron_L2'] = 'Glutamatergic excitatory'
pre_adata.obs.loc[pre_adata.obs['louvain'] == '88', 'neuron_L2'] = 'Glutamatergic excitatory'
pre_adata.obs.loc[pre_adata.obs['louvain'] == '89', 'neuron_L2'] = 'Glutamatergic excitatory'

pre_adata.obs.loc[pre_adata.obs['louvain'] == '90', 'neuron_L2'] = 'Splatter'
pre_adata.obs.loc[pre_adata.obs['louvain'] == '91', 'neuron_L2'] = 'GABAergic inhibitory'
pre_adata.obs.loc[pre_adata.obs['louvain'] == '92', 'neuron_L2'] = 'Splatter'
pre_adata.obs.loc[pre_adata.obs['louvain'] == '93', 'neuron_L2'] = 'Splatter'
pre_adata.obs.loc[pre_adata.obs['louvain'] == '94', 'neuron_L2'] = 'Splatter'
pre_adata.obs.loc[pre_adata.obs['louvain'] == '95', 'neuron_L2'] = 'GABAergic inhibitory'
pre_adata.obs.loc[pre_adata.obs['louvain'] == '96', 'neuron_L2'] = 'Splatter'
pre_adata.obs.loc[pre_adata.obs['louvain'] == '97', 'neuron_L2'] = 'Upper rhombic lip'
pre_adata.obs.loc[pre_adata.obs['louvain'] == '98', 'neuron_L2'] = 'Splatter'
pre_adata.obs.loc[pre_adata.obs['louvain'] == '99', 'neuron_L2'] = 'Glutamatergic excitatory'

pre_adata.obs.loc[pre_adata.obs['louvain'] == '100', 'neuron_L2'] = 'GABAergic inhibitory'
pre_adata.obs.loc[pre_adata.obs['louvain'] == '101', 'neuron_L2'] = 'Splatter'
pre_adata.obs.loc[pre_adata.obs['louvain'] == '102', 'neuron_L2'] = 'Splatter'
pre_adata.obs.loc[pre_adata.obs['louvain'] == '103', 'neuron_L2'] = 'Splatter'
pre_adata.obs.loc[pre_adata.obs['louvain'] == '104', 'neuron_L2'] = 'GABAergic inhibitory'



# adata.obs['neuron_L2'] = ''
# adata.obs.loc[pre_adata.obs.index, 'neuron_L2'] = pre_adata.obs.loc[pre_adata.obs.index, 'neuron_L2']

# adata.write_h5ad('../output_TACA2/neuron_final/integration_PROTEIN_CODING_reorganized_32bits_21datasets_neuron_batch_min10_HVG2000_epoch100_latent50_lr5e-4_res20_UMAP_L2_label_022425.h5ad')



excit_adata = pre_adata[pre_adata.obs['neuron_L2'] == 'Glutamatergic excitatory']
excit_adata.write_h5ad('../output_TACA2/neuron_final/preintegration_PROTEIN_CODING_reorganized_32bits_21datasets_excitatory_neuron_022425.h5ad')

inhib_adata = pre_adata[pre_adata.obs['neuron_L2'] == 'GABAergic inhibitory']
inhib_adata.write_h5ad('../output_TACA2/neuron_final/preintegration_PROTEIN_CODING_reorganized_32bits_21datasets_inhibitory_neuron_022425.h5ad')






