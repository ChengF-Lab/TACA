GENE_SELECT_INTERSECTION = "INTERSECTION"
GENE_SELECT_UNION = "UNION"
GENE_SELECT_PROTEIN_CODING = "PROTEIN_CODING"

SUPERCTs = {'Astrocyte', 'Immune', 'Neuron', 'Oligodendrocyte', 'OPC', 'Others'}


SUBCTs = {
    'InhibNeuron': set(['Sst', 'Pvalb', 'Vip', 'Lamp5', 'Lamp5 Lhx6', 'Sncg', 'Chandelier', 'Pax6', 'Sst Chodl']),
    'ExcitNeuron': set(['L2/3 IT', 'L5 IT', 'L4 IT', 'L6 IT', 'L6 CT', 'L6 IT Car3', 'L6b', 'L5/6 NP', 'L5 ET']),
}


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




# CT2SUBCT = {
#     'nonregex':{

#     }
# }



BR_MAPPING = {
    'Prefrontal cortex': 'Frontal lobe',
    'Orbitofrontal cortex': 'Frontal lobe', 
    'dorsolateral prefrontal cortex': 'Frontal lobe',
    'occipital cortex': 'Occipital lobe',
    'middle temporal gyrus': 'Temporal lobe',
    'occipitotemporal cortex': 'Occipital lobe',
    'cerebral cortex': 'Cerebral cortex',
    'hippocampal formation': 'Limbic lobe',
    'somatosensory cortex': 'Parietal lobe',
    'BA9': 'Frontal lobe',
    'substantia nigra pars compacta': 'Basal ganglia',
    'BA4': 'Frontal lobe',
    'cerebral nuclei': 'Basal ganglia',
    'substantia nigra pars compacta (SNpc) of the midbrain': 'Basal ganglia',
    'prefrontal cortex': 'Frontal lobe',
    'entorhinal cortex': 'Temporal lobe',
    'thalamic complex': 'Thalamus',
    'midbrain': 'Midbrain',
    'pons': 'Pons',
    'Cerebellum': 'Cerebellum',
    'cerebellum': 'Cerebellum',
    'myelencephalon': 'Medulla',
    'ALS motor cortex': 'Frontal lobe',
    'hypothalamus': 'Hypothalamus',
    'Anterior cingulate cortex': 'Limbic lobe',
    'Primary visual cortex': 'Occipital lobe',
    'Dorsolateral prefrontal cortex': 'Frontal lobe',
    'Primary somatosensory cortex': 'Parietal lobe',
    'Primary auditory cortex': 'Temporal lobe',
    'Middle temporal gyrus': 'Temporal lobe',
    'spinal cord': 'Spinal cord',
    'Angular gyrus': 'Parietal lobe',
    'cervical spinal cord white matter': 'Spinal cord',
    'Primary motor cortex': 'Frontal lobe',
    'White matter': 'Unspecified',
    'Brodmann (1909) area 4': 'Frontal lobe',
    'white matter of cerebellum': 'Cerebellum'
}



