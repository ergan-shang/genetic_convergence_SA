
high_num_pair <- fread(file = here::here(work_directory, 
                                         "yao_2023/report/B1_interesting_pairs/greater_than_4_modules.csv"))

batch_name <- '51_pairwise'
processed_results <- fread(file = here::here(work_directory,
                                             'yao_2023/data/intermediate_data', batch_name, 'processed_results.csv'))

processed_results <- processed_results[, .(comparison, p_value, active_group)]
processed_results[, c("gene1", "gene2") := tstrsplit(comparison, "_vs_")]

merged_results <- merge(processed_results, processed_results, 
                        by.x = c("gene1", "gene2"), by.y = c("gene2", "gene1"))

merge(x = high_num_pair,
      y = merged_results[, .(gene1, gene2, active_group.x, active_group.y)],
      all.x = TRUE)

#Venn diagrams
genes_interest <- c("TRAF6", "TRIB1", "HSP90B1")
processed_results[
    gene1 %in% genes_interest &
    gene2 %in% genes_interest,
][order(gene1, gene2)] |>
  print()

genes_interest <- c("STAT1", "STAT2", "TYK2")
processed_results[
    gene1 %in% genes_interest &
    gene2 %in% genes_interest,
][order(gene1, gene2)] |>
  print()
