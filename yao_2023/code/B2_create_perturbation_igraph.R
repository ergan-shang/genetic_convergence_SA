
high_num_pair <- fread(file = paste0(work_directory, 
                                     "/yao_2023/report/B1_interesting_pairs/greater_than_3_modules.csv"))



edges_dt <- copy(high_num_pair)[
  , .(u = pmin(gene1, gene2), v = pmax(gene1, gene2))
]

# If the data can contain duplicate rows for the same pair, create an edge weight = count
edges_weighted <- edges_dt[, .(weight = .N), by = .(u, v)]

# Build the graph (undirected). Vertices are inferred from edge endpoints.
g <- graph_from_data_frame(edges_weighted, directed = FALSE)

# (Optional) remove self-loops if any slipped in
g <- simplify(g, remove.loops = TRUE, edge.attr.comb = "first")

# Quick checks
print(g)
cat("Number of vertices:", vcount(g), "\n")
cat("Number of edges:", ecount(g), "\n")

dir.create(paste0(work_directory, 
                  "/yao_2023/report/B2_perturbation_igraph/"), recursive = TRUE)
pdf(file = paste0(work_directory, 
                  "/yao_2023/report/B2_perturbation_igraph/greater_than_3_modules.pdf"), 
                  width = 10, height = 8)

set.seed(1)
plot(
  g,
  vertex.size = 15,
  vertex.color = "#f5b70a",
  vertex.frame.color = "#013e75",
  vertex.label.cex = 1,
  vertex.label.color = "#222222",
  edge.width = 2,
  edge.color = "#013e75"
)

dev.off()

# 
# batch_name <- '51_pairwise'
# processed_results <- fread(glue(work_directory, '/yao_2023/data/intermediate_data/', batch_name, '/processed_results.csv'))
# # Assuming `prossed_results` is a data.table
# 
# processed_results <- processed_results[, .(comparison, p_value, active_group)]
# processed_results[, c("gene1", "gene2") := tstrsplit(comparison, "_vs_")]
# 
# merged_results <- merge(processed_results, processed_results, 
#                         by.x = c("gene1", "gene2"), by.y = c("gene2", "gene1"))
# 
# merge(x = high_num_pair,
#       y = merged_results[, .(gene1, gene2, active_group.x, active_group.y)],
#       all.x = TRUE)
# 
# 
# processed_results[gene1 == 'TRAF6' & gene2 == "TRIB1", ]
# processed_results[gene1 == 'TRIB1' & gene2 == "TRAF6", ]
# 
# processed_results[gene1 == 'STAT1' & gene2 == "STAT2", ]
# processed_results[gene1 == 'STAT2' & gene2 == "STAT1", ]

