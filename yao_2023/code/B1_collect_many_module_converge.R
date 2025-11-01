library(glue)
library(data.table)
library(ggplot2)

work_directory <- here::here()
source(here::here(work_directory, 'R/collect_and_structure_results.R'))

batch_name <- '51_pairwise'
processed_results <- fread(glue(work_directory, 
                                '/yao_2023/data/intermediate_data/', batch_name, '/processed_results.csv'))

processed_results <- processed_results[, .(comparison, p_value, active_group)]
processed_results[, c("gene1", "gene2") := tstrsplit(comparison, "_vs_")]

merged_results <- merge(processed_results, processed_results, 
                        by.x = c("gene1", "gene2"), by.y = c("gene2", "gene1"))
# merged_results[, p_min := pmin(p_value.x, p_value.y)]
merged_results[, p_max := pmax(p_value.x, p_value.y)]

threshold <- 0.025
merged_results[, is_significant := 0]
merged_results[, is_significant := is_significant + (p_max < threshold)]
merged_results[, is_significant := is_significant + (p_max < threshold/2450)]
merged_results[, is_significant := factor(is_significant, levels = c(0, 1, 2))]

p_value_dt <- merged_results[gene1 < gene2, .(gene1, gene2, is_significant)]


parse_nums <- function(x) {
  if (is.na(x) || trimws(x) == "") return(integer(0))
  parts <- unlist(strsplit(x, "/", fixed = TRUE))
  nums  <- suppressWarnings(as.integer(trimws(parts)))
  nums  <- nums[!is.na(nums)]
  unique(nums)
}

count_number_strings <- function(s1, s2, mode = c("union", "intersection")) {
  mode <- match.arg(mode)
  a <- parse_nums(s1)
  b <- parse_nums(s2)
  if (mode == "union") {
    length(union(a, b))
  } else {
    length(intersect(a, b))
  }
}

merged_results[, relevant_module :=
                 mapply(count_number_strings, active_group.x, active_group.y,
                        MoreArgs = list(mode = "intersection"))]
module_num_dt <- merged_results[gene1 > gene2, .(gene1, gene2, relevant_module)]

p_value_n_module_num_dt <- merge(p_value_dt, module_num_dt, 
                                 by.x = c("gene1", "gene2"),
                                 by.y = c("gene2", "gene1"),)
dir.create(paste0(work_directory, 
                  "/yao_2023/report/B1_interesting_pairs/"))
high_num_pair <- p_value_n_module_num_dt[relevant_module >= 4 & is_significant == 2, ]
fwrite(high_num_pair, file = paste0(work_directory, 
                                    "/yao_2023/report/B1_interesting_pairs/greater_than_4_modules.csv"))

high_num_pair <- p_value_n_module_num_dt[relevant_module >= 3 & is_significant == 2, ]
fwrite(high_num_pair, file = paste0(work_directory, 
                                    "/yao_2023/report/B1_interesting_pairs/greater_than_3_modules.csv"))

high_num_pair <- p_value_n_module_num_dt[relevant_module >= 2 & is_significant == 2, ]
fwrite(high_num_pair, file = paste0(work_directory, 
                                    "/yao_2023/report/B1_interesting_pairs/greater_than_2_modules.csv"))

