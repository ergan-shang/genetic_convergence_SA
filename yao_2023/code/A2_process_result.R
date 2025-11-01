
out_dir <- file.path(work_directory, "yao_2023", "data", "intermediate_data", "A1_chen_method")
stopifnot(dir.exists(out_dir))

# List all saved results
files <- list.files(out_dir, pattern = "^(.+)_module_([0-9]+)\\.rds$", full.names = TRUE)
if (!length(files)) stop(paste("No .rds files found under:", out_dir))

# Helper to parse filename: <treatment>_module_<number>.rds
parse_fname <- function(path) {
  bn <- basename(path)
  m <- regexec("^(.+)_module_([0-9]+)\\.rds$", bn)
  parts <- regmatches(bn, m)[[1]]
  if (length(parts) != 3L) return(list(treatment_name = NA_character_, module_index = NA_integer_))
  list(treatment_name = parts[2], module_index = as.integer(parts[3]))
}

# Read and collect
rows <- lapply(files, function(f) {
  meta <- parse_fname(f)
  result <- readRDS(f)
  pval <- result$pval
  data.table(
    treatment_name = meta$treatment_name,
    module_index   = meta$module_index,
    p_value        = pval
  )
})

summary_dt <- rbindlist(rows, use.names = TRUE, fill = TRUE)
setorder(summary_dt, treatment_name, module_index)

clustering_file <- glue(work_directory, "yao_2023/data/module_list_df_2000_genes.csv")
clustering <- data.table(read.csv(clustering_file))
clustering$gene_name <- gsub("\\.", "-", clustering$gene_name)
num_module <- max(clustering$cluster_index)

num_pert <- 50
correction_num <- (num_pert * (num_pert - 1) / 2) * num_module

summary_dt[, is.significant := p_value < 0.05 / correction_num]

all_treatment_names <- unique(summary_dt$treatment_name)
treatment_pairs <- expand.grid(treatment1 = all_treatment_names, treatment2 = all_treatment_names)
treatment_pairs <- as.data.table(treatment_pairs)
treatment_pairs <- treatment_pairs[treatment1 != treatment2, ]
treatment_pairs[, comparison := paste0(treatment1, '_vs_', treatment2)]
treatment_pairs[, p_value := 1]
treatment_pairs[, active_group := ""]

for (pair_index in seq_len(nrow(treatment_pairs))) {
  t1 <- treatment_pairs[pair_index, treatment1]
  t2 <- treatment_pairs[pair_index, treatment2]
  
  # rows for just these two treatments that are significant
  sum_pair <- summary_dt[
    treatment_name %in% c(t1, t2) & is.significant == TRUE,
    .N,
    by = module_index
  ]
  
  has_overlap <- any(sum_pair$N == 2)  # module appears for both t1 and t2
  # update by reference
  if(has_overlap){
    treatment_pairs[pair_index, p_value := 1e-10]
    sum_pair <- sum_pair[N == 2, ]
    sum_pair <- sum_pair[order(sum_pair, decreasing = FALSE), ]
    convergence_modules <- paste(sort(sum_pair$module_index), collapse = "/")
    treatment_pairs[pair_index, active_group := convergence_modules]
  }
}

# Save outputs
rds_path <- file.path(out_dir, "processed_data.rds")
saveRDS(treatment_pairs, file = rds_path)
hist(treatment_pairs$p_value)
