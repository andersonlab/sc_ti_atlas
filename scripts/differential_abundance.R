suppressMessages({
  library(optparse)
  library(Matrix)
  library(readr)
  library(dplyr)
  library(variancePartition)
  library(edgeR)
  library(ggplot2)
  library(tidyr)
  library(forcats)
})


# This function calculates the abundance per sample and cluster from single_cell data and returns a DGEList object
# You need the metadata and pass the columns with the cluster and sample info
# Additionally, when creating the object, you can select columns (besides sample) that will be included in the sample information
# You can also filter samples based on their total number of cells
calculate_abundances <- function(df, cluster_col = NA, sample_col = NA, metadata_cols = c(), 
                                 filter_samples = F, n_cells = 0){
  # Filter samples with more than n_cells cells
  if(filter_samples){
    message('Filtering samples with less than ', n_cells , ' cells')
    n <- length(unique(data[[sample_col]]))
    df <- df %>%
      group_by(sanger_sample_id) %>%
      filter(n() > n_cells) 
    message(length(unique(data[[sample_col]])), '/', n, ' samples were kept')
  }
  
  print(nrow(df))
  # Calculate abundances per cluster and sample
  message('Calculating abundances')
  df %>%
    dplyr::count(.data[[sample_col]], .data[[cluster_col]]) %>%  # Fix column referencing
    pivot_wider(names_from = sample_col, values_from = 'n', values_fill = 0) %>%
    tibble::column_to_rownames(cluster_col)  -> abundances
  
  # Add sample if it was not there 
  metadata_cols <- unique(c(sample_col, metadata_cols))
  
  # Create metadata for samples
  df %>%
    ungroup() %>%
    distinct(across(all_of(metadata_cols))) %>% 
    tibble::column_to_rownames(sample_col) -> metadata
  metadata <- metadata[colnames(abundances),]
  
  # Create object
  gExpr <- DGEList(abundances, samples=metadata, lib.size = colSums(abundances))
  gExpr <- calcNormFactors(gExpr)
  
  return(gExpr)
}


run_one <- function(path) {
  # Read data
  data <- readr::read_csv(path, show_col_types = FALSE) %>% 
    mutate(fractioned = !is.na(f2_f1_ratio)) %>% 
    # filter(!(disease_status == 'CD' & inflammation_status_binary == 'uninflamed')) # Uncomment out for CD inflamed vs healthy analysis

  # Calculate abundances
  gExpr <- calculate_abundances(
    data,
    cluster_col  = "label",
    sample_col   = "sanger_sample_id",
    metadata_cols = c("sanger_sample_id", "age_binned", "sex",
                      "ses_cd", "disease_status", "patient_id", 'fractioned'),
    filter_samples = FALSE,
    n_cells = 500
  )
  
  # Fit model
  vobjDream <- voomWithDreamWeights(gExpr, ~disease_status + age_binned + sex + fractioned,
                                    gExpr$samples)
  fitmm     <- dream(vobjDream, ~disease_status + age_binned + sex + fractioned, gExpr$samples)
  fitmm     <- eBayes(fitmm)
  
  # Be a bit robust to the exact coefficient name (e.g., "disease_statusHealthy")
  coef_names  <- colnames(fitmm$coefficients)
  coef_target <- grep("^disease_status", coef_names, value = TRUE)
  if (length(coef_target) == 0) stop("No coefficient for disease_status found in model.")
  
  results <- variancePartition::topTable(
    fitmm, coef = coef_target[1], number = Inf, adjust.method = "fdr"
  ) %>% 
    as.data.frame() %>% 
    mutate(coef = "disease_status",
           formula = "~disease_status + age_binned + sex + fractioned") %>% 
    as_tibble(rownames = "Cell_types")
  
  results
}

# Paths for each cohort
files <- c(
  Discovery  = "data/anndata/obs/ti-Discovery_obs.csv",
  Replication = "data/anndata/obs/ti-Replication_obs.csv",
  Full = "data/anndata/obs/ti-Full_obs.csv"
)

# Run for both and return two data frames in a named list
results_list <- lapply(files, run_one)

disc <- results_list$Discovery %>%
  select(Cell_types, logFC, adj.P.Val) %>%
  dplyr::rename(logFC_Discovery = logFC, adj.P.Val_Discovery = adj.P.Val)

repl <- results_list$Replication %>%
  select(Cell_types, logFC, adj.P.Val) %>%
  dplyr::rename(logFC_Replication = logFC, adj.P.Val_Replication = adj.P.Val)

significant <- full_join(disc, repl, by = "Cell_types") %>% 
  mutate(cohorts_agree =  sign(logFC_Discovery) == sign(logFC_Replication) &
           adj.P.Val_Discovery < 0.05 &
           adj.P.Val_Replication < 0.05) %>% 
  filter(cohorts_agree == TRUE)

full <- results_list$Full %>%
  select(Cell_types, logFC, adj.P.Val) %>%
  mutate(manual_annotation = fct_reorder(manual_annotation, logFC, .desc = FALSE),
         cohorts_agree = Cell_types %in% significant$Cell_types)

# Save results
write_csv(full, './results/ti_differential_abundance_results_full.csv')


