load_results <- function(result_path) {
  if (file.exists(result_path)) {
    return(readRDS(result_path))
  } else {
    warning(glue("File {result_path} not found!"))
    return(NULL)
  }
}

process_result <- function(result, group = NULL) {
  if (is.null(result)) return(NULL)
  
  active_features <- collect_active_features(result, group = group)
  print(active_features)
  active_features[[1]] <- paste(active_features[[1]], collapse = "/")
  active_features[[2]] <- paste(active_features[[2]], collapse = "/")
  return(list(
    p_value = result$p_value,
    active_features = active_features[[1]],
    active_group = active_features[[2]]
  ))
  
}

parse_download_results <- function() {
  processed_results_copy <- copy(processed_results)[, -"active_features"]
  
  # Split comparison into two separate columns
  processed_results_copy[, c("gene_1", "gene_2") := tstrsplit(comparison, "_vs_", fixed = TRUE)]
  
  # Convert active_group to a list of numbers
  processed_results_copy[, active_group_list := strsplit(active_group, "/")]
  
  # Expand each row to have one active group per row
  processed_results_expanded <- processed_results_copy[, .(gene_1, gene_2, p_value, active_feature = unlist(active_group_list)), by = .(comparison)]
  
  # Convert active_feature to numeric
  processed_results_expanded[, active_feature := as.integer(active_feature)]
  
  return(processed_results_expanded)
}


plot_active_features <- function(parsed_results, work_directory, batch_name, threshold = 10) {
  # Define the output directory for plots
  plot_dir <- glue("{work_directory}/report/{batch_name}")

  # Create the directory if it doesn’t exist
  if (!dir.exists(plot_dir)) {
    dir.create(plot_dir, recursive = TRUE)
  }

  # Get the count of each active_feature
  feature_counts <- table(parsed_results$active_feature)

  # Filter features with large counts
  significant_features <- as.integer(names(feature_counts[feature_counts > threshold]))

  # Get unique gene pairs
  unique_gene_pairs <- unique(parsed_results[, .(gene_1, gene_2)])

  for (feature in significant_features) {
    # Extract rows where active_feature == feature
    active_feature_subset <- parsed_results[active_feature == feature,]

    # Merge with unique gene pairs, filling unmatched rows with NA
    gene_pair_merged <- merge(unique_gene_pairs,
                              active_feature_subset[, .(gene_1, gene_2, active_feature)],
                              by = c("gene_1", "gene_2"),
                              all.x = TRUE)

    # Create a new column indicating presence of active_feature
    gene_pair_merged[, has_active_feature := !is.na(active_feature)]

    # Generate the plot
    p <- ggplot(gene_pair_merged, aes(x = gene_1, y = gene_2, color = has_active_feature)) +
      geom_point(size = 3, alpha = 0.7) +
      scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "red")) +
      theme_bw() +
      labs(title = glue("Gene Pair Activation Status for Active Module {feature}"),
           x = "Exploratory perturbation",
           y = "Confirmative perturbation",
           color = "Active Feature") +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

    # Define the filename and save the plot
    plot_filename <- glue("{plot_dir}/active_feature_{feature}.png")
    ggsave(plot_filename, plot = p, width = 8, height = 6, dpi = 300)

    # Print message to confirm saving
    message(glue("Plot saved: {plot_filename}"))
  }
}
