# classify_cardiac_immune.R
#
# Cardiac-specific k-nearest neighbor immune cell classifier
# Integrates query data with a curated cardiac single-cell atlas via Harmony,
# projects SMOTE-balanced synthetic training cells, and assigns cell type labels
# with confidence thresholding.
#
# Authors: Meite I, Pesce L, Feinstein M, Thorp EB, Lantz C
# Manuscript: Meite I et al., submitted to Am J Physiol Heart Circ Physiol (2026)
# Code DOI: 10.5281/zenodo.23112514
#
# Usage:
#   source("classify_cardiac_immune.R")
#   results <- classify_cardiac_immune(
#     query_seurat = your_seurat_object,
#     atlas_path   = "path/to/seurat_atlas.rds",
#     train_path   = "path/to/train_smote.rds"
#   )


classify_cardiac_immune <- function(
    query_seurat,
    atlas_path           = NULL,
    train_path           = NULL,
    gene_map             = NULL,
    best_k               = 5,
    prob_threshold       = 0.5,
    cell_type_col        = "lineage_r.08",
    return_probabilities = TRUE
) {
  
  if (is.null(atlas_path)) stop("Please provide a path to the atlas RDS file via atlas_path.")
  if (is.null(train_path)) stop("Please provide a path to the SMOTE-balanced training RDS file via train_path.")
  
  library(Seurat)
  library(harmony)
  library(class)
  library(FNN)
  
  detect_gene_scheme <- function(genes) {
    pct_ensembl <- mean(grepl("^ENSMUSG|^ENSG", genes))
    if (pct_ensembl > 0.5) return("ensembl") else return("symbol")
  }
  
  convert_ensembl_to_symbol <- function(genes, gene_map) {
    idx <- match(genes, gene_map$ensembl_gene_id)
    converted <- gene_map$mgi_symbol[idx]
    converted[is.na(converted)] <- genes[is.na(converted)]
    return(converted)
  }
  
  convert_symbol_to_ensembl <- function(genes, gene_map) {
    idx <- match(genes, gene_map$mgi_symbol)
    converted <- gene_map$ensembl_gene_id[idx]
    converted[is.na(converted)] <- genes[is.na(converted)]
    return(converted)
  }
  
  cat("Cardiac Immune Classifier\n")
  
  # [1/7] Load data
  cat("[1/7] Loading atlas and training data...\n")
  if (!file.exists(atlas_path)) stop("Atlas not found at: ", atlas_path)
  if (!file.exists(train_path)) stop("Training data not found: ", train_path)
  
  atlas       <- readRDS(atlas_path)
  train_atlas <- readRDS(train_path)
  DefaultAssay(atlas) <- "RNA"
  
  if (!cell_type_col %in% colnames(atlas@meta.data)) {
    stop("cell_type_col '", cell_type_col, "' not found. Available: ",
         paste(colnames(atlas@meta.data), collapse = ", "))
  }
  atlas$cell_type_use <- atlas@meta.data[[cell_type_col]]
  cat(sprintf("  Using label column: %s\n", cell_type_col))
  cat(sprintf("  Label preview: %s\n",
              paste(head(unique(atlas$cell_type_use), 5), collapse = ", ")))
  
  if (!is.null(gene_map)) {
    cat(sprintf("  Gene map loaded: %s entries\n", format(nrow(gene_map), big.mark = ",")))
  } else {
    cat("  No gene map provided; gene name conversion disabled\n")
  }
  
  cat(sprintf("  Atlas: %s cells\n",         format(ncol(atlas),        big.mark = ",")))
  cat(sprintf("  Training data: %s cells\n", format(ncol(train_atlas),  big.mark = ",")))
  cat(sprintf("  Query data: %s cells\n",    format(ncol(query_seurat), big.mark = ",")))
  
  # [2/7] Detect gene name schemes
  cat("\n[2/7] Detecting gene name schemes...\n")
  
  query_scheme <- detect_gene_scheme(rownames(query_seurat))
  atlas_scheme <- detect_gene_scheme(rownames(atlas))
  train_scheme <- detect_gene_scheme(rownames(train_atlas))
  
  cat(sprintf("  Query genes: %s (e.g. %s)\n", query_scheme, rownames(query_seurat)[1]))
  cat(sprintf("  Atlas genes: %s (e.g. %s)\n", atlas_scheme, rownames(atlas)[1]))
  cat(sprintf("  Train genes: %s (e.g. %s)\n", train_scheme, rownames(train_atlas)[1]))
  
  if (!is.null(gene_map) & atlas_scheme != query_scheme) {
    cat(sprintf("  Converting atlas: %s to %s\n", atlas_scheme, query_scheme))
    if (atlas_scheme == "ensembl") {
      new_genes <- convert_ensembl_to_symbol(rownames(atlas), gene_map)
    } else {
      new_genes <- convert_symbol_to_ensembl(rownames(atlas), gene_map)
    }
    atlas <- RenameGenesSeurat(atlas, newnames = new_genes)
    cat("  Atlas gene names converted\n")
  } else {
    cat("  Atlas and query schemes match\n")
  }
  
  real_train_cells  <- intersect(colnames(train_atlas), colnames(atlas))
  synth_train_cells <- setdiff(colnames(train_atlas),   colnames(atlas))
  cat(sprintf("  Real train cells:      %s\n", format(length(real_train_cells),  big.mark = ",")))
  cat(sprintf("  Synthetic train cells: %s\n", format(length(synth_train_cells), big.mark = ",")))
  
  # [3/7] Merge query + atlas
  cat("\n[3/7] Merging query with atlas (real cells only)...\n")
  query_seurat$dataset <- "query"
  atlas$dataset        <- "atlas"
  combined <- merge(query_seurat, atlas)
  cat(sprintf("  Combined: %s cells\n", format(ncol(combined), big.mark = ",")))
  
  cat("  Joining layers...\n")
  DefaultAssay(combined) <- "RNA"
  combined <- JoinLayers(combined)
  cat("  Layers joined\n")
  
  # [4/7] Harmony integration
  cat("\n[4/7] Running Harmony on real cells only...\n")
  combined <- NormalizeData(combined, verbose = FALSE)
  combined <- FindVariableFeatures(combined, nfeatures = 3000, verbose = FALSE)
  combined <- ScaleData(combined, verbose = FALSE)
  combined <- RunPCA(combined, npcs = 50, verbose = FALSE)
  combined <- RunHarmony(combined,
                         group.by.vars = "dataset",
                         theta    = 2,
                         dims.use = 1:50,
                         verbose  = FALSE)
  cat("  Harmony integration complete\n")
  
  all_harmony        <- Embeddings(combined, reduction = "harmony")
  query_harmony      <- all_harmony[colnames(query_seurat), 1:50]
  atlas_harmony      <- all_harmony[colnames(atlas),        1:50]
  real_train_harmony <- all_harmony[real_train_cells,       1:50]
  
  # [5/7] Project SMOTE synthetic cells
  cat("\n[5/7] Projecting SMOTE synthetic cells into Harmony space...\n")
  
  if (length(synth_train_cells) > 0) {
    
    synth_expr_full <- GetAssayData(train_atlas, assay = "RNA", layer = "counts")
    
    if (!is.null(gene_map) & train_scheme != query_scheme) {
      cat(sprintf("  Converting train genes: %s to %s\n", train_scheme, query_scheme))
      if (train_scheme == "ensembl") {
        new_train_genes <- convert_ensembl_to_symbol(rownames(synth_expr_full), gene_map)
      } else {
        new_train_genes <- convert_symbol_to_ensembl(rownames(synth_expr_full), gene_map)
      }
      valid_genes <- !is.na(new_train_genes) & !duplicated(new_train_genes)
      synth_expr_full <- synth_expr_full[valid_genes, ]
      rownames(synth_expr_full) <- new_train_genes[valid_genes]
      cat(sprintf("  Converted: %s valid, %s failed/duplicate removed\n",
                  sum(valid_genes), sum(!valid_genes)))
      cat(sprintf("  Sample converted genes: %s, %s, %s\n",
                  rownames(synth_expr_full)[1],
                  rownames(synth_expr_full)[2],
                  rownames(synth_expr_full)[3]))
    }
    
    shared_genes <- intersect(rownames(synth_expr_full), rownames(combined))
    cat(sprintf("  Shared genes for projection: %s\n", format(length(shared_genes), big.mark = ",")))
    
    if (length(shared_genes) > 0) {
      
      synth_expr <- synth_expr_full[shared_genes, synth_train_cells]
      
      cat(sprintf("  Synthetic cells to project: %s\n",
                  format(length(synth_train_cells), big.mark = ",")))
      
      col_sums <- Matrix::colSums(synth_expr)
      nonzero  <- col_sums > 0
      if (any(nonzero)) {
        synth_expr[, nonzero] <- log1p(
          t(t(synth_expr[, nonzero]) / col_sums[nonzero]) * 10000
        )
      }
      
      pca_rotation     <- Loadings(combined, reduction = "pca")
      shared_genes_pca <- intersect(shared_genes, rownames(pca_rotation))
      cat(sprintf("  Genes used for PCA projection: %s\n",
                  format(length(shared_genes_pca), big.mark = ",")))
      
      synth_pca <- t(as.matrix(synth_expr[shared_genes_pca, ])) %*%
        pca_rotation[shared_genes_pca, 1:50]
      
      synth_pca_nas <- rowSums(is.na(synth_pca)) > 0
      if (any(synth_pca_nas)) {
        cat(sprintf("  Removing %s synthetic cells with NA PCA coordinates\n",
                    sum(synth_pca_nas)))
        synth_pca         <- synth_pca[!synth_pca_nas, ]
        synth_train_cells <- synth_train_cells[!synth_pca_nas]
      }
      
    } else {
      cat("  0 shared genes - skipping synthetic cell projection\n")
      synth_train_cells <- character(0)
      synth_pca         <- matrix(0, nrow = 0, ncol = 50)
    }
    
    atlas_pca           <- Embeddings(combined, "pca")[colnames(atlas), 1:50]
    valid_atlas         <- rowSums(is.na(atlas_pca)) == 0
    atlas_pca_clean     <- atlas_pca[valid_atlas, ]
    atlas_harmony_clean <- atlas_harmony[valid_atlas, ]
    cat(sprintf("  Atlas cells with valid PCA: %s / %s\n",
                sum(valid_atlas), length(valid_atlas)))
    
    if (length(synth_train_cells) > 0) {
      
      k_proj <- 5
      nn_idx <- FNN::get.knnx(
        data  = atlas_pca_clean,
        query = synth_pca,
        k     = k_proj
      )$nn.index
      
      synth_harmony <- matrix(0, nrow = length(synth_train_cells), ncol = 50)
      for (i in seq_len(nrow(nn_idx))) {
        synth_harmony[i, ] <- colMeans(atlas_harmony_clean[nn_idx[i, ], ])
      }
      rownames(synth_harmony) <- synth_train_cells
      cat(sprintf("  Projected %s synthetic cells\n",
                  format(length(synth_train_cells), big.mark = ",")))
      
      real_labels  <- atlas$cell_type_use[match(real_train_cells,  colnames(atlas))]
      synth_labels <- train_atlas$cell_type[match(synth_train_cells, colnames(train_atlas))]
      
      cat("\nReal cell type distribution (from atlas):\n")
      print(sort(table(real_labels)))
      cat("\nSynthetic cell type distribution (from train_smote):\n")
      print(sort(table(synth_labels)))
      
      train_harmony <- rbind(real_train_harmony, synth_harmony)
      train_labels  <- c(real_labels, synth_labels)
      
    } else {
      cat("  No valid synthetic cells - using real train cells only\n")
      train_harmony <- real_train_harmony
      train_labels  <- atlas$cell_type_use[match(real_train_cells, colnames(atlas))]
    }
    
  } else {
    cat("  No synthetic cells - using real train cells only\n")
    train_harmony <- real_train_harmony
    train_labels  <- atlas$cell_type_use[match(real_train_cells, colnames(atlas))]
  }
  
  cat(sprintf("\n  Final training set: %s cells\n", format(nrow(train_harmony), big.mark = ",")))
  cat("Final training set cell type distribution:\n")
  print(sort(table(train_labels)))
  
  # [6/7] k-NN classification
  cat("\n[6/7] Applying k-NN classifier...\n")
  
  valid_train <- rowSums(is.na(train_harmony)) == 0
  valid_query <- rowSums(is.na(query_harmony)) == 0
  
  if (any(!valid_train) | any(!valid_query)) {
    cat(sprintf("  Removing %d train cells with NA harmony coords\n", sum(!valid_train)))
    cat(sprintf("  Removing %d query cells with NA harmony coords\n", sum(!valid_query)))
  }
  
  train_harmony <- train_harmony[valid_train, ]
  train_labels  <- train_labels[valid_train]
  query_harmony <- query_harmony[valid_query, ]
  
  predictions <- knn(
    train = train_harmony,
    test  = query_harmony,
    cl    = train_labels,
    k     = best_k,
    prob  = TRUE
  )
  
  max_prob <- attr(predictions, "prob")
  cat(sprintf("  Classification complete (k=%d)\n", best_k))
  
  # [7/7] Flag low confidence
  cat("\n[7/7] Flagging low-confidence predictions...\n")
  low_confidence      <- max_prob < prob_threshold
  predictions_flagged <- as.character(predictions)
  predictions_flagged[low_confidence] <- "Unassigned"
  
  cat(sprintf("  Probability threshold: %.2f\n", prob_threshold))
  cat(sprintf("  Assigned:   %s cells (%.1f%%)\n",
              format(sum(!low_confidence), big.mark = ","),
              100 * sum(!low_confidence) / length(low_confidence)))
  cat(sprintf("  Unassigned: %s cells (%.1f%%)\n",
              format(sum(low_confidence),  big.mark = ","),
              100 * sum(low_confidence)  / length(low_confidence)))
  
  cat("Classification complete\n")
  
  cat("Predicted cell type distribution:\n")
  type_counts <- sort(table(predictions_flagged), decreasing = TRUE)
  for (i in seq_len(min(10, length(type_counts)))) {
    cat(sprintf("  %s: %d\n", names(type_counts)[i], type_counts[i]))
  }
  if (length(type_counts) > 10) {
    cat(sprintf("  ... and %d more cell types\n", length(type_counts) - 10))
  }
  
  all_query_cells  <- colnames(query_seurat)
  predictions_full <- rep(NA_character_, length(all_query_cells))
  max_prob_full    <- rep(NA_real_,      length(all_query_cells))
  confidence_full  <- rep(NA_character_, length(all_query_cells))
  
  kept_idx <- match(rownames(query_harmony), all_query_cells)
  predictions_full[kept_idx] <- predictions_flagged
  max_prob_full[kept_idx]    <- max_prob
  confidence_full[kept_idx]  <- ifelse(low_confidence, "Low", "High")
  
  results <- data.frame(
    cell_id            = all_query_cells,
    predicted_celltype = predictions_full,
    max_probability    = max_prob_full,
    confidence         = confidence_full,
    stringsAsFactors   = FALSE
  )
  
  return(results)
}