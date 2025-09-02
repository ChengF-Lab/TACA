library(monocle3) 
library(ggplot2)
library(dplyr)
library(Seurat)
library(ggplot2)
library(patchwork) 

so <- readRDS('/home/houy2/isilon/Cheng-Qiu/TACA2/Integration_test/output_TACA2/immune_final2/micro_macro.rds')

dis<-"AD"  # Resilience PART AD
sex_type<-"all"

so <- so[, so$immune_L3 != 'Macrophage']
so <- so[, !(so$immune_L3 %in% c("Macrophage", "Monocyte", "DAM2"))]

# remove 2 reference datasets
so <- so[, so$Dataset != '37824655']
so <- so[, so$Dataset != '37824663']

sub_so <- so[,so$New_Group == dis]
#sub_so <- sub_so[, sub_so$Sex == sex_type]

# Create gene annotations
gene_annotation <- data.frame(gene_short_name = rownames(sub_so))
rownames(gene_annotation) <- rownames(sub_so)

# MODIFY HERE: Include New_Group and Sex in cell metadata
cell_metadata <- data.frame(
  barcode = sub_so@assays[["RNA"]]@counts@Dimnames[[2]], 
  batch = sub_so$BatchID, 
  celltype = sub_so$immune_L3,
  New_Group = sub_so$New_Group,  # Add New_Group
  Sex = sub_so$Sex             # Add Sex
)
rownames(cell_metadata) <- sub_so@assays[["RNA"]]@counts@Dimnames[[2]]

expression_matrix <- sub_so@assays[["RNA"]]@counts

cds_from_seurat <- new_cell_data_set(expression_matrix,
                                     cell_metadata = cell_metadata,
                                     gene_metadata = gene_annotation)


recreate.partition <- c(rep(1, length(cds_from_seurat@colData@rownames)))
names(recreate.partition) <- cds_from_seurat@colData@rownames
recreate.partition <- as.factor(recreate.partition)

print(str(recreate.partition))
cds_from_seurat@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partition


list_cluster <- sub_so@meta.data[[sprintf("immune_L3")]]
names(list_cluster) <- sub_so@assays[["RNA"]]@data@Dimnames[[2]]

cds_from_seurat@clusters@listData[["UMAP"]][["clusters"]] <- list_cluster

cds_from_seurat@clusters@listData[["UMAP"]][["louvain_res"]] <- "NA"

cds_from_seurat@int_colData$reducedDims@listData[["UMAP"]] <- sub_so@reductions[["umap"]]@cell.embeddings


get_earliest_principal_node <- function(cds_from_seurat, celltype){
  cell_ids <- which(colData(cds_from_seurat)[, "celltype"] == celltype)
  
  closest_vertex <- cds_from_seurat@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds_from_seurat), ])
  root_pr_nodes <- igraph::V(principal_graph(cds_from_seurat)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))]
  
  return(root_pr_nodes)
}

# find all possible partitions
all_partitions <- unique(cds_from_seurat@clusters$UMAP$partitions)
all_partitions <- all_partitions[all_partitions != "1"]
# set all partitions to 1
cds_from_seurat@clusters$UMAP$partitions[cds_from_seurat@clusters$UMAP$partitions %in% all_partitions] <- "1"

cds_from_seurat <- learn_graph(cds_from_seurat, close_loop = FALSE)
# cds_from_seurat <- learn_graph(cds_from_seurat)
cds_from_seurat <- order_cells(cds_from_seurat, reduction_method = "UMAP",  root_pr_nodes=get_earliest_principal_node(cds_from_seurat, celltype="Homeostasis Microglia1"))
celltype <- c(cds_from_seurat@colData@listData$celltype)
pseudotime <- c(cds_from_seurat@principal_graph_aux@listData$UMAP$pseudotime)
cells <- c(colnames(cds_from_seurat))

# VERIFY - Check that New_Group and Sex were added correctly
print("Checking New_Group and Sex in cds_from_seurat:")
print(table(cds_from_seurat@colData$New_Group))
print(table(cds_from_seurat@colData$Sex))

#=============================================
# Original input path
PATH <- '/home/houy2/isilon/Cheng-Hou/scRNA_AD/integ_TACA2/test/monocle3'

# Save the cds object for future use
saveRDS(cds_from_seurat, file = file.path(PATH, paste0(dis,"_",sex_type,"_cds_from_seurat.rds")))


# MODIFY HERE: Include Sex and New_Group in the saved trajectory data
df <- data.frame(cells, pseudotime, celltype, Sex=cds_from_seurat@colData$Sex, New_Group=cds_from_seurat@colData$New_Group)
write.table(df, 
            file = file.path(PATH, paste0(dis,"_",sex_type,"_micro_macro_trajectory_no_loop.tsv")), 
            sep = "\t", 
            row.names = FALSE, 
            col.names = TRUE,  # Changed to TRUE to include column names
            quote = FALSE)

# Plot trajectory
pdf(file.path(PATH, paste0(dis,"_",sex_type,"_micro_macro_trajectory_no_loop.pdf")), 
    width = 9, height = 6)
plot_cells(cds_from_seurat,
           color_cells_by = "pseudotime",
           norm_method = "size_only",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.0, 
           rasterize = TRUE)
dev.off()

# Create a more detailed plot showing the trajectory and cell types
pdf(file.path(PATH, paste0(dis,"_",sex_type,"_trajectory_with_celltypes.pdf")), 
    width = 10, height = 8)
plot_cells(cds_from_seurat,
           color_cells_by = "celltype",
           label_cell_groups = TRUE,
           label_leaves = TRUE,
           label_branch_points = TRUE,
           group_label_size = 4,
           cell_size = 1)
dev.off()

# NEW: Add trajectory plot colored by Sex
pdf(file.path(PATH, paste0(dis,"_",sex_type,"_trajectory_by_sex.pdf")), 
    width = 10, height = 8)
plot_cells(cds_from_seurat,
           color_cells_by = "Sex",
           label_cell_groups = FALSE,
           label_leaves = FALSE,
           label_branch_points = FALSE,
           group_label_size = 4,
           cell_size = 1)
dev.off()

# NEW: Add trajectory plot colored by New_Group (if multiple groups are present)
if(length(unique(cds_from_seurat@colData$New_Group)) > 1) {
  pdf(file.path(PATH, paste0(dis,"_",sex_type,"_trajectory_by_group.pdf")), 
      width = 10, height = 8)
  plot_cells(cds_from_seurat,
             color_cells_by = "New_Group",
             label_cell_groups = FALSE,
             label_leaves = FALSE,
             label_branch_points = FALSE,
             group_label_size = 4,
             cell_size = 1)
  dev.off()
}

#======================= Find driver genes================
# First check graph structure
print("Graph structure check:")
print(length(principal_graph(cds_from_seurat)[["UMAP"]]))
print(igraph::gsize(principal_graph(cds_from_seurat)[["UMAP"]]))

# Filter to expressed genes to reduce computational load
expression_matrix <- cds_from_seurat@assays@data$counts
expressed_genes <- rowSums(expression_matrix > 0) >= 10
cds_filtered <- cds_from_seurat[expressed_genes,]
print(paste("Testing", sum(expressed_genes), "genes out of", length(expressed_genes)))

# Run with single core first for debugging
ciliated_genes <- graph_test(cds_filtered, 
                             neighbor_graph="principal_graph", 
                             cores=4)

# Order genes by q-value
ciliated_genes_ordered <- ciliated_genes %>% 
  arrange(q_value)

# Save results
write.csv(ciliated_genes_ordered, 
          file = file.path(PATH, paste0(dis,"_",sex_type,"_driver_genes.csv")))
