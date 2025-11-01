library(glue)
library(data.table)
library(doParallel)
library(foreach)
library(highmean)


work_directory <- here::here()
residual_subset <- readRDS(glue(work_directory, 
                                '/yao_2023/data/intermediate_data/residual_matrix_all_in_paper.rds'))

residual_subset <- residual_subset[, -c("ID", "Cell_cycle_phase")]
residual_subset <- data.table(residual_subset)
colnames(residual_subset) <- gsub("\\.", "-", colnames(residual_subset))

clustering_file <- glue(work_directory, "/yao_2023/data/module_list_df_2000_genes.csv")
clustering <- data.table(read.csv(clustering_file))
clustering$gene_name <- gsub("\\.", "-", clustering$gene_name)

# Subset control samples
set.seed(1)
control_subsetting_sample_size <- 500
control <- residual_subset[Guides_collapsed_by_gene == "non-targeting", -"Guides_collapsed_by_gene"]
control <- control[sample(1:nrow(control), control_subsetting_sample_size), ]
control <- control[, match(clustering$gene_name, colnames(control)), with = FALSE]

dir.create(glue(work_directory, '/yao_2023/data/intermediate_data/A1_chen_method/'), 
           recursive = TRUE, showWarnings = FALSE)
dir.create(glue(work_directory, '/yao_2023/log/'), recursive = TRUE, showWarnings = FALSE)

preprocess_one_setting <- function(residual_subset, treatment_name, clustering) {
  X <- residual_subset[Guides_collapsed_by_gene == treatment_name, -"Guides_collapsed_by_gene"]

  X <- X[, match(clustering$gene_name, colnames(X)), with = FALSE]
  return(X)
}

all_treatment_names <- unique(residual_subset$Guides_collapsed_by_gene)
all_treatment_names <- all_treatment_names[all_treatment_names != 'non-targeting']

log_file <- glue(work_directory, '/yao_2023/log/logs.txt')

for(treatment_name in all_treatment_names){
  
  # set up parallel backend
  ncores <- max(1L, parallel::detectCores() - 1L)
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  
  res_list <- foreach(
    module_index = unique(clustering$cluster_index),
    .packages = c("data.table", "highmean", "glue"),
    .export   = c("preprocess_one_setting", "residual_subset", "treatment_name",
                  "clustering", "control", "work_directory")
  ) %dopar% {
    treatment <- preprocess_one_setting(residual_subset, treatment_name, clustering)
    gene_compared <- clustering[cluster_index == module_index, ]$gene_name
    
    control_one_module <- as.matrix(control[, ..gene_compared])
    treatment_one_module <- as.matrix(treatment[, ..gene_compared])
    
    test_result <- tryCatch(
      apval_Chen2010(control_one_module, treatment_one_module, eq.cov = FALSE),
      error = function(e) NA_real_
    )
    
    output_filename <- paste0(work_directory, 
                              '/yao_2023/data/intermediate_data/A1_chen_method/', treatment_name, '_module_', module_index,'.rds')
    saveRDS(test_result, file = output_filename)
  }
  stopCluster(cl)
}
