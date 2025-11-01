
source(file.path(work_directory, 'R/collect_and_structure_results.R'))

# Define batch name and read in processed data
batch_name <- '51_pairwise'
processed_results <- fread(file.path(work_directory, 
                                     'yao_2023/data/intermediate_data/', batch_name, '/processed_results.csv'))

# Parse results from a custom function
parsed_results <- parse_download_results()

# Choose mode for selecting significant pairs ("both" requires both directions; "either" requires at least one)
mode <- "both"  

# Create output directory if it doesn't exist
output_dir <- file.path(work_directory, 
                        "yao_2023/report", paste0("54_", batch_name, "_igraph"))
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Remove all existing plots in the output directory
file.remove(list.files(output_dir, full.names = TRUE))

# Loop through all unique module indices
for (module_index in unique(parsed_results$active_feature)) {
  
  dt <- parsed_results[active_feature == module_index & p_value < 0.05]
  if (nrow(dt) == 0) next  # Skip if no significant pairs
  
  # Standardize gene pair ordering
  dt[, gene_low := pmin(gene_1, gene_2)]
  dt[, gene_high := pmax(gene_1, gene_2)]
  dt[, gene_pair := paste(gene_low, gene_high, sep = "_")]
  
  # Select pairs based on mode
  if (mode == "both") {
    selected_pairs <- dt[, .N, by = gene_pair][N == 2]
    dt_bidirectional <- dt[gene_pair %in% selected_pairs$gene_pair]
    
    dt_edges <- dt_bidirectional[, .(
      gene_1 = unique(gene_low),
      gene_2 = unique(gene_high),
      p_value = mean(p_value)
    ), by = gene_pair]
    
  } else if (mode == "either") {
    dt_plot <- unique(dt, by = "gene_pair")
    dt_edges <- dt_plot[, .(
      gene_1 = gene_low,
      gene_2 = gene_high,
      p_value = p_value
    ), by = gene_pair]
  } else {
    stop("mode must be either 'both' or 'either'")
  }
  
  # Skip if no edges
  if (nrow(dt_edges) == 0) next
  
  # Add edge weight
  dt_edges[, edge_width := -log10(p_value)]
  
  # Build igraph object
  g <- graph_from_data_frame(dt_edges[, .(gene_1, gene_2, edge_width)], directed = FALSE)
  set.seed(123)
  layout <- layout_with_fr(g)
  
  # Styling
  # scaled_width <- rescale(dt_edges$edge_width, to = c(1, 8))
  scaled_width <- pmin(dt_edges$edge_width, 3)

  
  edge_colors <- colorRampPalette(c("lightgray", "darkred"))(100)
  edge_col_indices <- as.integer(rescale(dt_edges$edge_width, to = c(1, 100)))
  edge_col <- edge_colors[edge_col_indices]
  
  node_degree <- degree(g)
  scaled_size <- rescale(node_degree, to = c(10, 25))
  scaled_label_cex <- rescale(node_degree, to = c(0.3, 0.8))
  
  vertex_labels <- V(g)$name
  vertex_labels[node_degree < 2] <- NA
  
  # Skip if no vertex names (e.g., empty graph)
  if (length(vertex_labels) == 0 || all(is.na(vertex_labels))) next
  
  node_count <- length(V(g))
  size_tag <- if (node_count >= 10) {
    "A"
  } else if (node_count >= 6) {
    "B"
  } else {
    "C"
  }
  
  # Save plot to file with size tag
  plot_path <- file.path(output_dir, paste0(size_tag, "_module_", module_index, "_mode_", mode, ".pdf"))
  
  # Load annotation for this module
  annotation_path <- file.path(work_directory, 
                               "yao_2023/data/final_data/module_GO/module_", module_index, ".csv")
  if (file.exists(annotation_path)) {
    annotation_dt <- fread(annotation_path)
    annotation_dt <- annotation_dt[p.adjust < 0.05,]
    annotation_text <- annotation_dt$Description[1:min(20, .N)]
    if(length(annotation_text) == 0){
      annotation_text <- "No annotation available."
    }
  } else {
    annotation_text <- "No annotation available."
  }
  
  pdf(file = plot_path, width = 8, height = 7, useDingbats = FALSE)  # adjust width/height as needed
  
  plot(
    g,
    layout = layout,
    vertex.label = vertex_labels,
    vertex.label.cex = scaled_label_cex,
    vertex.label.color = "black",
    vertex.size = scaled_size,
    vertex.frame.color = "black",
    vertex.frame.width = 1.5,
    vertex.color = "#fdd835",
    edge.width = scaled_width,
    edge.color = edge_col,
    main = glue("Perturbation convergence on module {module_index}"),
    margin = 0.2
  )
  
  # Position annotation in the top-left of the plotting area
  usr_coords <- par("usr")  # get current plot coordinates
  
  # Adjust starting position (top-left)
  x_pos <- usr_coords[1] + 0.01 * (usr_coords[2] - usr_coords[1])
  y_pos <- usr_coords[4] - 0.01 * (usr_coords[4] - usr_coords[3])
  
  # Line height offset
  line_spacing <- strheight("M", cex = 0.7) * 1.2
  
  # Draw each line
  for (i in seq_along(annotation_text)) {
    text(
      x = x_pos,
      y = y_pos - (i - 1) * line_spacing,
      labels = annotation_text[i],
      adj = c(0, 1),
      cex = 0.7,
      col = "gray30"
    )
  }
  
  dev.off()
}
