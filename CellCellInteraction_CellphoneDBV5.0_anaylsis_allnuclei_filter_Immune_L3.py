import os
import scanpy as sc
import datetime

path = "xxx"
path_out = os.path.join(path,"CCI_output_cellphonedbv5.0_250220")
path_out_data = os.path.join(path_out,"Data")
path_out_data_seq = os.path.join(path_out_data,"data_seq")
path_out_fig = os.path.join(path_out,"Figures")
start_time = datetime.datetime.now()
print("Script started at:", start_time)


# METHOD 2. Statistical inference of interaction specificity
from cellphonedb.src.core.methods import cpdb_statistical_analysis_method


threshold = 0.01
iterations = 10000
threads=1
result_precision = 4
cpdb_results = cpdb_statistical_analysis_method.call(
        cpdb_file_path = "xxx",
        meta_file_path = os.path.join(path_out_data_seq,"TACA2_data_meta_250220_AD_allnuclei_filter_Immune_L3.tsv"),
        counts_file_path = os.path.join(path_out_data_seq,"TACA2_data_normalized_log1p_matrix_250220_AD_AnndataV0.10.9_allnuclei_filter_Immune.h5ad"),
        counts_data = "ensembl",
        score_interactions = True,
        threshold = threshold,
        iterations = iterations,
        threads = threads,
        debug_seed = 42,                                 # debug randome seed. To disable >=0.
        result_precision = result_precision,                            # Sets the rounding for the mean values in significan_means.
        pvalue = 0.05,                                   # P-value threshold to employ for significance.
        subsampling = False,                             # To enable subsampling the data (geometri sketching).
        subsampling_log = False,                         # (mandatory) enable subsampling log1p for non log-transformed data inputs.
        subsampling_num_pc = 100,                        # Number of componets to subsample via geometric skectching (dafault: 100).
        subsampling_num_cells = 1000,                    # Number of cells to subsample (integer) (default: 1/3 of the dataset).
        separator = '|',                                 # Sets the string to employ to separate cells in the results dataframes "cellA|CellB".
        debug = False,                                   # Saves all intermediate tables employed during the analysis in pkl format.
        output_path = os.path.join(path_out_data,"TACA2_250220_AD_celltypeL3_threshold"+str(threshold)+"_iterations"+str(iterations)+"_threads"+str(threads)+"_allnuclei_filter_Immune"),
        output_suffix = None
)



end_time = datetime.datetime.now()
print("Script ended at:", end_time)

# Duration
print("")
duration = end_time - start_time
print("Duration of script execution:", duration)





# METHOD 3. Retrieval of differentially expressed interactions
from cellphonedb.src.core.methods import cpdb_degs_analysis_method

# cpdb_results = cpdb_degs_analysis_method.call(
#     cpdb_file_path="/home/doul2/beegfs/doul2/Work/database/CCI/cellphonedbv5.0.0/cellphonedb.zip",
#     meta_file_path=os.path.join(path_out_data_seq, "0_ALS_multiome_meta_Subcelltype_241220.tsv"),
#     counts_file_path=os.path.join(path_out_data_seq, "0_ALS_multiome_normalized_matrix_241220.h5ad"),
#          degs_file_path = os.path.join(path_out_data_seq, "0_ALS_multiome_DEGs_subcelltype_seurat_MAST_var_nCount_nFeature_241220.tsv"),
#          counts_data = 'hgnc_symbol',
#          threshold = 0.1,
#          output_path = path_out_data,
#         output_suffix = "2_ALS_multiome_cpdb_deg_MAST_var_nCount_nFeature"
# )



