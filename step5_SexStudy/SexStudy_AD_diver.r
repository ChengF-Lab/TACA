# --- Function to create correlation scatter plots for female-specific driver genes ---
create_female_driver_gene_plots <- function(cds_from_seurat, PATH, dis, gene_id_file, 
                                            correlation_threshold = 0.3, p_value_threshold = 0.05) {
  # --- 0. Setup ---
  gene_ids <- read.delim(gene_id_file)
  gene_ids <- gene_ids[!duplicated(gene_ids$gene_id), ]
  
  # Create mapping from gene_id to gene_name for output
  gene_id_to_name <- setNames(gene_ids$gene_name, gene_ids$gene_id)
  
  # Create output directory if it doesn't exist
  output_dir <- file.path(PATH, "driver_genes")
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  # --- 1. Extract data ---
  pseudotime <- cds_from_seurat@principal_graph_aux$UMAP$pseudotime
  celltype <- cds_from_seurat$celltype
  sex <- cds_from_seurat$Sex
  
  # Standardize sex labels
  sex <- tolower(sex)
  
  data <- data.frame(Pseudotime = pseudotime, CellType = celltype, Sex = sex)
  
  # --- 2. Create bins and calculate proportions ---
  num_bins <- 50
  breaks <- seq(min(data$Pseudotime), max(data$Pseudotime), length.out = num_bins + 1)
  data$PseudotimeBin <- cut(data$Pseudotime, breaks = breaks, labels = 1:num_bins, include.lowest = TRUE)
  data$PseudotimeBin <- as.numeric(as.character(data$PseudotimeBin))
  
  # --- 3. Filter data to DAM1 cells ---
  dam1_data <- data[data$CellType == "DAM1", ]
  
  # --- 4. Get expression data ---
  if("normalized" %in% names(cds_from_seurat@assays@data)){
    expression_matrix <- as.matrix(cds_from_seurat@assays@data$normalized)
  } else if ("counts" %in% names(cds_from_seurat@assays@data)){
    raw_counts <- as.matrix(cds_from_seurat@assays@data$counts)
    size_factors <- colSums(raw_counts)
    size_factors <- size_factors / median(size_factors)
    expression_matrix <- sweep(raw_counts, 2, size_factors, "/")
    expression_matrix <- log1p(expression_matrix)
  } else {
    stop("Neither 'normalized' nor 'counts' found in assay names.")
  }
  
  # --- 5. Read in driver genes from graph_test results ---
  all_driver_genes_file <- file.path(PATH, paste0(dis, "_all_driver_genes.csv"))
  driver_genes_df <- read.csv(all_driver_genes_file)
  
  # Use direct column names
  gene_col <- "gene_short_name"
  qval_col <- 'q_value'
  
  # Filter by q-value
  significant_driver_genes <- driver_genes_df[[gene_col]][driver_genes_df[[qval_col]] < p_value_threshold]
  significant_driver_genes <- significant_driver_genes[significant_driver_genes %in% rownames(expression_matrix)]
  
  # --- 6. Calculate DAM1 proportion along pseudotime ---
  # Sex-specific proportions
  sex_proportions <- data %>%
    group_by(PseudotimeBin, Sex, CellType) %>%
    summarise(Count = n(), MeanPseudotime = mean(Pseudotime), .groups = 'drop') %>%
    group_by(PseudotimeBin, Sex) %>%
    mutate(Total = sum(Count), Proportion = Count / Total) %>%
    ungroup()
  
  # Filter to DAM1 sex-specific proportions
  dam1_sex_proportions <- sex_proportions %>% filter(CellType == "DAM1")
  
  # Split by sex
  female_dam1_props <- dam1_sex_proportions %>% filter(Sex == "female")
  
  # --- 7. Calculate gene expression along pseudotime for female ---
  # For female AD
  female_gene_expr <- matrix(0, nrow = length(significant_driver_genes), ncol = num_bins)
  rownames(female_gene_expr) <- significant_driver_genes
  colnames(female_gene_expr) <- 1:num_bins
  
  for (bin in 1:num_bins) {
    bin_cells <- rownames(data[data$PseudotimeBin == bin & data$CellType == "DAM1" & data$Sex == "female", ])
    if (length(bin_cells) > 0) {
      female_gene_expr[, bin] <- rowMeans(expression_matrix[significant_driver_genes, bin_cells, drop = FALSE])
    }
  }
  
  # --- 8. Get proportion vector for correlation calculation ---
  # Female AD proportion
  female_prop_vector <- numeric(num_bins)
  for (bin in 1:num_bins) {
    bin_data <- female_dam1_props %>% filter(PseudotimeBin == bin)
    if (nrow(bin_data) > 0) female_prop_vector[bin] <- bin_data$Proportion
  }
  
  # --- 9. Calculate correlation for each gene for female ---
  # For female AD
  female_correlations <- numeric(length(significant_driver_genes))
  female_p_values <- numeric(length(significant_driver_genes))
  names(female_correlations) <- significant_driver_genes
  names(female_p_values) <- significant_driver_genes
  
  for (i in 1:length(significant_driver_genes)) {
    gene <- significant_driver_genes[i]
    female_expr <- female_gene_expr[gene, ]
    female_valid <- !is.na(female_expr) & !is.na(female_prop_vector) & !is.infinite(female_expr) & female_prop_vector > 0
    
    if (sum(female_valid) > 5) {
      corTest <- cor.test(female_expr[female_valid], female_prop_vector[female_valid], method = "pearson")
      female_correlations[i] <- corTest$estimate
      female_p_values[i] <- corTest$p.value
    } else {
      female_correlations[i] <- NA
      female_p_values[i] <- NA
    }
  }
  
  # --- 10. Filter significant female-specific genes ---
  # For female - ONLY POSITIVE CORRELATIONS
  significant_female_genes <- significant_driver_genes[
    !is.na(female_correlations) & 
      !is.na(female_p_values) & 
      female_correlations >= correlation_threshold &  # Only positive correlations
      female_p_values < p_value_threshold
  ]
  
  # If no female-specific genes found, use top female-positive genes
  if (length(significant_female_genes) == 0) {
    # Use genes with positive correlation
    positive_female_genes <- significant_driver_genes[female_correlations > 0 & !is.na(female_correlations)]
    sorted_female_genes <- positive_female_genes[order(female_correlations[positive_female_genes], decreasing = TRUE)]
    significant_female_genes <- head(sorted_female_genes, 10)
  }
  
  # Get top 10 female-specific genes
  top_female_genes <- significant_female_genes[order(female_correlations[significant_female_genes], decreasing = TRUE)]
  top_female_genes <- head(top_female_genes, min(10, length(top_female_genes)))
  
  # --- 11. Create correlation scatter plots for each of the top 10 female-specific genes ---
  # Convert gene IDs to gene names
  gene_names <- gene_id_to_name[top_female_genes]
  
  # Create a directory for individual plots
  scatter_dir <- file.path(output_dir, "female_gene_scatter_plots")
  if (!dir.exists(scatter_dir)) dir.create(scatter_dir, recursive = TRUE)
  
  # Set up for multi-plot PDF
  pdf(file.path(output_dir, paste0(dis, "_top10_female_genes_correlation_plots.pdf")), width = 12, height = 10)
  par(mfrow = c(3, 4))  # 3x4 grid for up to 12 plots
  
  # Individual scatter plots
  for (i in 1:length(top_female_genes)) {
    gene <- top_female_genes[i]
    gene_name <- gene_names[i]
    
    # Get expression values
    female_expr <- female_gene_expr[gene, ]
    
    # Filter out missing values
    valid <- !is.na(female_expr) & !is.na(female_prop_vector) & !is.infinite(female_expr) & female_prop_vector > 0
    x_values <- female_expr[valid]
    y_values <- female_prop_vector[valid]
    
    if (length(x_values) > 5) {  # Ensure enough data points
      # Create data frame for ggplot
      plot_data <- data.frame(
        Expression = x_values,
        Proportion = y_values
      )
      
      # Calculate correlation for title
      corr <- round(female_correlations[gene], 3)
      pval <- female_p_values[gene]
      pval_str <- if(pval < 0.001) "p < 0.001" else paste0("p = ", round(pval, 3))
      
      # Create individual scatter plot
      p <- ggplot(plot_data, aes(x = Expression, y = Proportion)) +
        geom_point(color = "#EC407A", alpha = 0.7) +
        geom_smooth(method = "lm", color = "#EC407A", fill = "#EC407A", alpha = 0.2) +
        labs(title = paste0(gene_name, " (", corr, ", ", pval_str, ")"),
             x = "Gene Expression", 
             y = "DAM1 Proportion") +
        theme_minimal() +
        theme(plot.title = element_text(size = 12, face = "bold"))
      
      # Save individual plot
      ggsave(file.path(scatter_dir, paste0(dis, "_", gene_name, "_scatter.pdf")), p, width = 5, height = 4)
      
      # For the combined PDF
      plot(x_values, y_values, 
           main = paste0(gene_name, " (r=", corr, ", ", pval_str, ")"),
           xlab = "Gene Expression", 
           ylab = "DAM1 Proportion",
           pch = 19, 
           col = "#EC407A")
      abline(lm(y_values ~ x_values), col = "#EC407A", lwd = 2)
    }
  }
  dev.off()
  
  # --- 12. Create a combined plot with small multiples ---
  # Set up data for combined plot
  combined_data <- data.frame()
  
  for (i in 1:length(top_female_genes)) {
    gene <- top_female_genes[i]
    gene_name <- gene_names[i]
    
    # Get expression values
    female_expr <- female_gene_expr[gene, ]
    
    # Filter out missing values
    valid <- !is.na(female_expr) & !is.na(female_prop_vector) & !is.infinite(female_expr) & female_prop_vector > 0
    x_values <- female_expr[valid]
    y_values <- female_prop_vector[valid]
    
    if (length(x_values) > 5) {  # Ensure enough data points
      # Create data for this gene and add to combined data
      gene_data <- data.frame(
        Expression = x_values,
        Proportion = y_values,
        Gene = gene_name
      )
      combined_data <- rbind(combined_data, gene_data)
    }
  }
  
  # Create the combined plot if we have data
  if (nrow(combined_data) > 0) {
    p_combined <- ggplot(combined_data, aes(x = Expression, y = Proportion)) +
      geom_point(color = "#EC407A", alpha = 0.7) +
      geom_smooth(method = "lm", color = "#EC407A", fill = "#EC407A", alpha = 0.2) +
      facet_wrap(~ Gene, scales = "free_x") +
      labs(title = paste0(dis, " - Top 10 Female-Specific Driver Genes"),
           x = "Gene Expression", 
           y = "DAM1 Proportion") +
      theme_minimal() +
      theme(plot.title = element_text(size = 14, face = "bold"))
    
    ggsave(file.path(output_dir, paste0(dis, "_top10_female_genes_combined.pdf")), 
           p_combined, width = 12, height = 10)
  }
  
  return(list(
    top_female_genes = top_female_genes,
    gene_names = gene_names
  ))
}

# --- Main Script ---
PATH <- '/home/houy2/isilon/Cheng-Hou/scRNA_AD/integ_TACA2/test/monocle3'
dis <- "AD"
gene_id_file <- "/home/houy2/isilon/Cheng-Hou/scRNA_AD/integ_TACA2/test/monocle3/EnsemblID2Symbol_39955.tsv" 

# Load the combined monocle object (with both sexes)
cds_from_seurat <- readRDS(file.path(PATH, paste0(dis, "_all_cds_from_seurat.rds")))

# Create the correlation plots
female_gene_results <- create_female_driver_gene_plots(
  cds_from_seurat, 
  PATH, 
  dis, 
  gene_id_file, 
  correlation_threshold = 0.3,
  p_value_threshold = 0.05
)

# Print the top female-specific genes that were used in the plots
cat("Created correlation plots for the following top female-specific driver genes:\n")
for (i in 1:length(female_gene_results$top_female_genes)) {
  cat(paste0(i, ". ", female_gene_results$gene_names[i], 
             " (", female_gene_results$top_female_genes[i], ")\n"))
}
