source(file.path(work_directory, 
                 'R/collect_and_structure_results.R'))
batch_name <- '51_pairwise'
processed_results <- fread(file.path(work_directory, 
                                     'yao_2023/data/intermediate_data/', batch_name, '/processed_results.csv'))
# Assuming `prossed_results` is a data.table

parsed_results <- parse_download_results()

print(table(parsed_results$active_feature))

parsed_results[active_feature == 31,]

# Example usage
library(ggplot2)
plot_active_features(parsed_results, work_directory, batch_name, threshold = 1)



