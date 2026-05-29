
install.packages("readr")
install.packages("dplyr")
install.packages("purrr")
install.packages("stringr")
library(readr)
library(dplyr)
library(purrr)
library(stringr)

# Mape, kur atrodas faili
in_dir <- "/home/evija.liva/Downloads/Data"
files <- list.files(in_dir, pattern = "\\.tsv$", full.names = TRUE)
basename(files)

#izvelk pacientu ID no failu nosaukumiem 
read_patient <- function(file) {
  sample_id <- stringr::str_extract(basename(file), "RS\\d+")
  
  df <- readr::read_tsv(file, show_col_types = FALSE) %>%
    dplyr::select(gene_name = 2, expression = 3) %>%
    dplyr::mutate(
      gene_name = as.character(gene_name),
      expression = as.numeric(expression)
    ) %>%
    dplyr::filter(!is.na(gene_name), gene_name != "") %>%
    
    dplyr::group_by(gene_name) %>%
    dplyr::summarise(expression = mean(expression, na.rm = TRUE), .groups = "drop") %>%
    
    dplyr::arrange(gene_name) %>%
    dplyr::rename(!!sample_id := expression)
  
  df
}

#Apvieno visus pacientu ID
patient_list <- lapply(files, read_patient)

#Gēnu ekspresijas matrica
exp_matrix <- reduce(patient_list, full_join, by = "gene_name")
any(is.na(exp_matrix))

exp_matrix <- as.data.frame(exp_matrix)
rownames(exp_matrix) <- exp_matrix$gene_name
exp_matrix <- as.matrix(exp_matrix[, -1])

#Matricas saglabāšana
write.csv(exp_matrix,
          file = "2026_BD_Gene_expression_matrix.csv",
          row.names = TRUE)

write.table(exp_matrix,
            file = "2026_DB_Gene_expression_matrix.tsv",
            sep = "\t",
            quote = FALSE,
            row.names = TRUE)

# Nolasa saglabāto matricu pirms filtrēšanas
mat_before <- read.csv("/home/evija.liva/2026_BD_Gene_expression_matrix.csv", 
                       row.names = 1)
nrow(mat_before)  # gēnu skaits pirms filtrēšanas
ncol(mat_before) 

#Gēnu filtrēšana
keep <- rowSums(exp_matrix <1 ) >= 9 #zemi ekspresēto gēnu izfiltrēšana
exp_matrix <- exp_matrix[!keep, ]
dim(exp_matrix)
install.packages('factoextra')
install.packages('ggplot2')
install.packages('ggfortify')

library(ggplot2)
library(ggfortify)
pca.plot <- autoplot(pca_exp_matrix, data = anno_match, colour = 'APAK_GRUPAS.ID')

#Visu pacientu anotācijas tabula
tabula <- read.csv("/home/evija.liva/Downloads/VISI_PACIENTI_LPCI.csv")
whole_annotation <- tabula[, c(1, 16)]
colnames(whole_annotation) <- c("PACIENTA.ID", "APAK_GRUPAS.ID")
whole_annotation <- whole_annotation[-(1:3), ]
whole_annotation <- whole_annotation[trimws(whole_annotation$APAK_GRUPAS.ID) != "", ]
rownames(whole_annotation) <- NULL
write.csv(whole_annotation,
          file = "2026_BD_Whole_annotation.csv",
          row.names = TRUE)

#Tumor transcriptome FCID pievienošana
library(readxl)
df <- read_excel("/home/evija.liva/Downloads/RS_pacientu_sekvenēšanas_flow_cells_lanes_barcodes.xls")
df_sub <- df %>%
  select(
    PACIENTA.ID = 1,
    Tumortranscriptome_FCID = 8
  ) %>%
  filter(!is.na(PACIENTA.ID))
df_sub <- df_sub[-((nrow(df_sub)-1):nrow(df_sub)), ]

df_sub_clean <- df_sub %>%
  mutate(
    PACIENTA.ID = trimws(as.character(PACIENTA.ID)),
    Tumortranscriptome_FCID = as.character(Tumortranscriptome_FCID)
  ) %>%
  filter(!is.na(PACIENTA.ID), PACIENTA.ID != "") %>%
  filter(PACIENTA.ID %in% whole_annotation$PACIENTA.ID) %>%
  group_by(PACIENTA.ID) %>%
  summarise(
    Tumortranscriptome_FCID =
      if (any(!is.na(Tumortranscriptome_FCID) & Tumortranscriptome_FCID != "")) {
        first(Tumortranscriptome_FCID[!is.na(Tumortranscriptome_FCID) & Tumortranscriptome_FCID != ""])
      } else {
        NA_character_
      },
    .groups = "drop"
  )
whole_annotation_FCID <- whole_annotation %>%
  mutate(PACIENTA.ID = trimws(as.character(PACIENTA.ID))) %>%
  left_join(df_sub_clean, by = "PACIENTA.ID")
nrow(whole_annotation)
nrow(whole_annotation_FCID)                      
sum(!is.na(whole_annotation_FCID$Tumortranscriptome_FCID))
sum(is.na(whole_annotation_FCID$Tumortranscriptome_FCID))

write.csv(whole_annotation_FCID,
          file = "2026_BD_Whole_annotation_FCID.csv",
          row.names = TRUE)

#PCA visiem pacientiem
whole_anno_match<-whole_annotation[match(colnames(exp_matrix),whole_annotation$PACIENTA.ID),]%>%na.omit()
whole_anno_match<-whole_annotation[whole_annotation$PACIENTA.ID%in%colnames(exp_matrix),]

exp_matrix_whole_anno<-exp_matrix[,match(whole_anno_match$PACIENTA.ID,colnames(exp_matrix))]
dim(exp_matrix_whole_anno)
identical(colnames(exp_matrix_whole_anno),whole_anno_match$PACIENTA.ID)


log_exp_matrix_whole <- log2(exp_matrix_whole_anno +1)
pca_exp_matrix_whole <- prcomp(t(log_exp_matrix_whole), scale=T) #PCA uz transpozētās matricas
pca_exp_matrix_whole$x
summary(pca_exp_matrix_whole)
pca.plot_whole <- autoplot(pca_exp_matrix_whole, data = whole_anno_match, colour = 'APAK_GRUPAS.ID')
pca.plot_whole
library(viridis)
pca.plot_whole +
  scale_color_viridis_d(option = "turbo", name = "Vēža tips") 

#PCA pēc FCID

whole_anno_match_FCID<-whole_annotation_FCID[match(colnames(exp_matrix),whole_annotation_FCID$PACIENTA.ID),]%>%na.omit()
whole_anno_match_FCID<-whole_annotation_FCID[whole_annotation_FCID$PACIENTA.ID%in%colnames(exp_matrix),]

exp_matrix_whole_anno_FCID<-exp_matrix[,match(whole_anno_match_FCID$PACIENTA.ID,colnames(exp_matrix))]
dim(exp_matrix_whole_anno_FCID)
identical(colnames(exp_matrix_whole_anno_FCID),whole_anno_match_FCID$PACIENTA.ID)


log_exp_matrix_whole_FCID <- log2(exp_matrix_whole_anno_FCID +1)
pca_exp_matrix_whole_FCID <- prcomp(t(log_exp_matrix_whole_FCID), scale=T) #PCA uz transponētās matricas
pca_exp_matrix_whole_FCID$x
summary(pca_exp_matrix_whole_FCID)

pca.plot_whole_FCID <- autoplot(pca_exp_matrix_whole_FCID, data = whole_anno_match_FCID, colour = 'Tumortranscriptome_FCID')
pca.plot_whole_FCID
pca.plot_whole_FCID +
  scale_color_viridis_d(option = "turbo", name = "Plūsmas kameras ID")

df_PC1 <- data.frame(
  PC1 = pca_exp_matrix_whole_FCID$x[, "PC1"],
  cancer = whole_anno_match_FCID$APAK_GRUPAS.ID,
  batch = whole_anno_match_FCID$Tumortranscriptome_FCID
)
residuals_PC1 <- df_PC1$PC1 - mean(df_PC1$PC1)
shapiro.test(residuals_PC1)
qqnorm(residuals_PC1)
qqline(residuals_PC1, col = "red")
kruskal.test(PC1 ~ cancer, data = df_PC1)


#LIONESS 
variance <- matrixStats::rowVars(as.matrix(log_exp_matrix_whole_FCID))
variance_2 <- data.frame(Seq = seq(1:nrow(log_exp_matrix_whole_FCID)), rowVars = variance [order (variance, decreasing = TRUE)])
ggplot(variance_2, aes(x = Seq, y = rowVars)) +geom_line() + scale_y_log10()

exp_top <- log_exp_matrix_whole_FCID[which(variance >= quantile(variance, c(0.75))), ]
lioness_res <-lioness(exp_top, motif = NULL, ppi = NULL, network.inference.method = "pearson", nscores = 1, mode = "union")
BALL_ind<-which(whole_anno_match_FCID$APAK_GRUPAS.ID=='BALL')
names(BALL)<-whole_anno_match_FCID[which(whole_anno_match_FCID$APAK_GRUPAS.ID=='BALL'),]$PACIENTA.ID

#NMV - permutācijas un empīrisko kvantiļu metodes variācija
flatten_lioness <- function(lioness_list) {
  # Get gene names
  genes <- rownames(lioness_list[[1]])
  
  # Create edge names (upper triangle only to avoid duplicates)
  edge_index <- which(upper.tri(matrix(1, length(genes), length(genes))), arr.ind = TRUE)
  edge_names <- paste(genes[edge_index[,1]], genes[edge_index[,2]], sep = "|")
  
  # Initialize output matrix
  E <- matrix(NA, nrow = length(lioness_list), ncol = length(edge_names))
  colnames(E) <- edge_names
  
  # Fill with edge weights
  for (i in seq_along(lioness_list)) {
    mat <- lioness_list[[i]]
    E[i, ] <- mat[edge_index]
  }
  
  return(E)
}


flat_lioness <- flatten_lioness(lioness_res)
dim(flat_lioness)

set.seed(2)
shapiro.test(sample(as.vector(flat_lioness), size=4000,replace=F))
hist(sample(as.vector(flat_lioness), size=4000,replace=F),breaks=100)
save.image()

huffle_genes<-function(gene_sample){
  print(dim(gene_sample))
  sam_perm<-apply(gene_sample,1,function(x){x<-x[sample(c(1:length(x)),size=length(x), replace=F)]
  })
  rownames(sam_perm)=colnames(gene_sample)
  colnames(sam_perm)=rownames(gene_sample)
  print(dim(sam_perm))
  lion_perm<-lioness(t(sam_perm), motif = NULL, ppi = NULL, network.inference.method = "pearson", ncores = 1)
  print(dim(lion_perm[[1]]))
  diff_cor_vecs<-lapply(lion_perm,function(x){x<-x[upper.tri(x)]
  return(x)})
  diff_cor_vecs<-unlist(diff_cor_vecs)
  print(paste0('overall we generated ',length(diff_cor_vecs),' random differential correlations'))
  
  return(diff_cor_vecs)}

gene_picks<-list()

for(i in c(1:50)){sam<-sample(c(1:nrow(exp_top)), size=300,replace=F)
sam<-exp_top[sam,]
cor_vecs<-replicate(100, shuffle_genes(sam))

gene_picks[[i]]<-cor_vecs}

gene_picks<-unlist(gene_picks)

hist(gene_picks)
length(gene_picks)
lim_pos<-quantile(gene_picks,0.975)  
lim_neg<-quantile(gene_picks,0.025)

test_ssn<-lioness_res[[1]]
test_ssn<-graph_from_adjacency_matrix(test_ssn,weighted=T,diag=F,mode='undirected')
test_ssn<-delete_edges(test_ssn,E(test_ssn)[E(test_ssn)$weight<lim_pos & E(test_ssn)$weight>lim_neg])

summary(E(test_ssn)$weight)

#transform network to unweighted
E(test_ssn)$weight<-1

cl<-cluster_louvain(test_ssn)
membership(cl)%>%table()

degree(test_graph)

lioness_list<-lapply(lioness_res,function(x){x<-x[upper.tri(x)]
return(x)})

length(unlist(lioness_list))

#Atlasa statistiski nozīmīgās šķautnes, Statistiski būtiskās diferenciālās koekspresijas robežvērtības -2.21 un 2.2
lim_pos <- 2.21
lim_neg <- -2.21

#Tīkli igraph ar nozīmīgajām taisnēm

lioness_graphs <- lapply(lioness_res, function(ssn) {
  graph <- graph_from_adjacency_matrix(ssn,weighted = TRUE,diag = FALSE,mode = "undirected")
  graph <- delete_edges(graph,E(graph)[E(graph)$weight < lim_pos & E(graph)$weight > lim_neg]) #atstāj būtiskās šķautnes
  graph <- delete_vertices(graph, degree(graph) == 0)  # noņem singletone nodes
  E(graph)$weight <- 1
  graph
})

saveRDS(graph, "graph.rds")

#Tikai BALL tīklu igraph
BALL_res <- lioness_res[BALL_ind]
names(BALL_res) <- whole_anno_match_FCID$PACIENTA.ID[BALL_ind]

BALL_graphs <- lioness_graphs[BALL_ind]
names(BALL_graphs) <- whole_anno_match_FCID$PACIENTA.ID[BALL_ind]

# virsotņu skaits katrā tīklā
num_vertices <- sapply(BALL_graphs, vcount)

# šķautņu skaits katrā tīklā
num_edges <- sapply(BALL_graphs, ecount)

# Virsotņu statistika
mean(num_vertices)
sd(num_vertices)
min(num_vertices)
max(num_vertices)

# Šķautņu statistika
mean(num_edges)
sd(num_edges)
min(num_edges)
max(num_edges)

# kopā vienā tabulā
network_summary <- data.frame(
  sample = names(BALL_graphs),
  vertices = num_vertices,
  edges = num_edges
)

network_stats <- data.frame(
  Patient = names(BALL_graphs),
  Vertices = sapply(BALL_graphs, vcount),
  Edges = sapply(BALL_graphs, ecount),
  Avg_clustering_coef = round(sapply(BALL_graphs, transitivity, type = "average"), 4),
  Avg_shortest_path = round(sapply(BALL_graphs, mean_distance, directed = FALSE, unconnected = TRUE), 4),
  Centralization = round(sapply(BALL_graphs, function(g) centr_degree(g)$centralization), 4),
  Diameter = sapply(BALL_graphs, diameter, unconnected = TRUE),
  Radius = sapply(BALL_graphs, radius),
  Density = round(sapply(BALL_graphs, edge_density), 4)
  )

write.csv(network_stats, 
          file = "BALL_network_stats.csv", 
          row.names = FALSE)

#Tīklu attēli

plot_list <- imap(BALL_graphs, function(graph, sample_id) {
  row <- network_summary[network_summary$sample == sample_id, ]
  
  ggraph(graph, layout = "fr") +
    geom_edge_link(alpha = 0.2, colour = "grey60") +
    geom_node_point(size = 1.5, colour = "#2C7FB8") +
    labs(
      title = sample_id,
      subtitle = paste0(
        "Virsotnes: ", row$vertices,
        " | Šķautnes: ", row$edges
      )
    ) +
    theme_void() +
    theme(
      plot.title = element_text(size = 9, hjust = 0.5),
      plot.subtitle = element_text(size = 8, hjust = 0.5)
    )
})

wrap_plots(plot_list)

id_key <- data.frame(
  sample = names(BALL_graphs),
  anon_id = sprintf("PBALL%02d", seq_along(BALL_graphs))
)
# Anonimizē pacientu ID
names(BALL_graphs) <- id_key$anon_id

# pievieno anonīmos ID summary tabulai
network_summary_anon <- merge(
  network_summary,
  id_key,
  by = "sample"
)

plot_list <- Map(function(graph, sample_id) {
  
  row <- network_summary_anon[
    network_summary_anon$anon_id == sample_id,
  ]
  
  ggraph(graph, layout = "kk") +
    geom_edge_link(alpha = 0.15, colour = "grey70") +
    geom_node_point(size = 1.5, colour = "#2C7FB8") +
    ggtitle(sample_id) +
    labs(
      subtitle = paste0(
        "V = ", row$vertices,
        " | E = ", row$edges
      )
    ) +
    theme_void() +
    theme(
      plot.title = element_text(size = 9, hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(size = 8, hjust = 0.5),
      
      
      panel.background = element_rect(
        fill = "white",
        colour = "grey60",
        linewidth = 0.5
      ),
      plot.background = element_rect(fill = "white", colour = NA),
      plot.margin = margin(5, 5, 5, 5)
    )
  
}, BALL_graphs, names(BALL_graphs))

grid.arrange(
  grobs = plot_list,
  ncol = 4
)

#Klasterizēšana

BALL_clusters <- lapply(BALL_graphs, function(graph) {
  cluster_louvain(graph, weights = E(graph)$weight, resolution = 1.5)
})

# BALL_clusters_test <- lapply(BALL_graphs, function(graph) {
#   cluster_louvain(graph, weights = E(graph)$weight, resolution = 1)
# })

# Klasteru skaits katram grafikam/tīklam
n_clusters_res15 <- sapply(BALL_clusters, function(cl) length(unique(membership(cl))))
n_clusters_res1 <- sapply(BALL_clusters_test, function(cl) length(unique(membership(cl))))

n_clusters_res15
n_clusters_res1
names(BALL_clusters) <- names(BALL_graphs)
names(BALL_clusters_test) <- names(BALL_graphs)

BALL_membership <- lapply(BALL_clusters, membership)
sizes(BALL_clusters[[1]])


BALL_graphs <- mapply(function(g, memb) {
  V(g)$membership <- memb[V(g)$name]
  g
}, BALL_graphs, BALL_membership, SIMPLIFY = FALSE) #pievieno atribūtu membership

#TOP2 klasteri
BALL_top2_cl <- lapply(BALL_clusters, function(cl) {
  cluster_sizes <- sizes(cl)
  top2 <- names(sort(cluster_sizes, decreasing = TRUE))[1:2]
  
  memb <- membership(cl)
  cluster_genes <- split(names(memb), memb)
  return(cluster_genes[top2])
})

#Tīkls lielākajiem klasteriem
genes_to_plot <- BALL_top2_cl[[15]][[2]]
sub_net <- induced_subgraph(BALL_graphs[[15]], vids = genes_to_plot)
plot(
  sub_net,
  vertex.label = V(sub_net)$name,
  vertex.size = 6,
  vertex.label.cex = 0.6,
  edge.width = 0.5,
  layout = layout_with_fr(sub_net),
)
sub_tbl <- as_tbl_graph(sub_net)
sub_tbl <- sub_tbl %>%
  activate(nodes) %>%
  mutate(
    degree = centrality_degree(),
    label = name
  )
p <- ggraph(sub_tbl, layout = "fr") +
  
  geom_edge_link(
    alpha = 0.25,
    linewidth = 0.3,
    colour = "grey60"
  ) +
  
  # Mezgli (nodes) — krāsoti pēc  pakāpes centralitātes
  geom_node_point(
    aes(size = degree, fill = degree),
    shape = 21,
    colour = "white",
    stroke = 0.5,
    alpha = 0.95
  ) +
  
  # Gēnu nosaukumi
  geom_node_text(
    aes(label = label),
    repel = TRUE,
    size = 2.5,
    colour = "grey10",
    fontface = "italic",
    family = "sans",
    max.overlaps = Inf,
    point.padding = unit(0.4, "lines"),
    box.padding = unit(0.3, "lines"),
    segment.colour = "grey50",
    segment.size = 0.25,
    segment.alpha = 0.5
  ) +
  
  # Krāsu skala mezgliem (zils → sarkans pēc centralitātes)
  scale_fill_gradientn(
    colours = c("#2166AC", "#4DAC26", "#D7191C"),
    name = "Pakāpes centralitāte",
    guide = guide_colorbar(
      barwidth = 6,
      barheight = 0.5,
      title.position = "top",
      direction = "horizontal"
    )
  ) +
  
  # Izmēra skala
  scale_size_continuous(
    range = c(2.5, 10),
    name = "Pakāpes centralitāte",
    guide = "none" 
  ) +
  
  theme_void(base_family = "sans") +
  theme(
    legend.position = c(0.9, 0.9),
    legend.direction = "horizontal",
    legend.title = element_text(
      face = "bold",
      size = 8,
      colour = "grey20", 
      hjust = 0.5
    ),
    legend.text = element_text(
      size = 7,
      colour = "grey20"
    ),
    plot.margin = margin(25, 25, 25, 25),
    plot.background = element_rect(
      fill = "white",
      colour = NA
    )
  )

p
ggsave(
  "cluster_network_15_.png",
  plot = p,
  width = 180,
  height = 160,
  units = "mm",
  dpi = 600,
  bg = "white"
)

sub_net_check <- induced_subgraph(BALL_graphs[[15]], BALL_top2_cl[[15]][[2]])
vcount(sub_net_check)
ecount(sub_net_check)
degree_centr <- degree(sub_net_check, mode = "all")
hub_df <- data.frame(
  gene = V(sub_net_check)$name,
  degree = degree_centr
)
hub_df <- hub_df[order(-hub_df$degree), ]

# Iegūt izmērus katram no top 2 klasteriem visos paraugos
lapply(BALL_top2_cl, function(x) {
  lapply(x, length)})

#Gene enrichement ar GO
# test_genes <- BALL_top2_cl[[1]]$genes[[1]] #TESTS
# 
# test_enrich <- enrichGO(
#   gene = test_genes,
#   keyType='SYMBOL',
#   OrgDb = org.Hs.eg.db,
#   ont = "BP",
#   universe = rownames(exp_top),
#   pAdjustMethod = "BH",
#   pvalueCutoff = 0.05,
#   readable = TRUE
# )
# dotplot(test_enrich, showCategory = 10) #TESTS

GO_genes <- unlist(BALL_top2_cl, recursive = FALSE)

Enrich_results <- lapply(GO_genes, function(genes) {
  enrichGO(
    gene = genes,
    keyType = "SYMBOL",
    OrgDb = org.Hs.eg.db,
    ont = "BP",
    universe = rownames(exp_top),
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05, 
    readable = TRUE
  )
})

GO_sig_res <- lapply(Enrich_results, \(x) {
  x2 <- x
  x2@result <- x2@result |>
    dplyr::filter(p.adjust < 0.05)
  x2
})
dir.create("GO_enrich_results", showWarnings = FALSE)

for (i in seq_along(GO_sig_res)) {
  file_name <- names(GO_sig_res)[i]
  
  if (is.null(file_name) || is.na(file_name) || file_name == "") {
    file_name <- paste0("GO_result_", i)
  }
  
  write.csv(
    as.data.frame(GO_sig_res[[i]]),
    file = file.path("GO_enrich_results", paste0(file_name, ".csv")),
    row.names = FALSE
  )
}

go_all_tbl <- purrr::imap_dfr(GO_sig_res, \(x, nm) {
  as.data.frame(x) |>
    as_tibble() |>
    mutate(source = nm)
})

go_all_tbl <- purrr::imap_dfr(GO_sig_res, \(x, nm) {
  as.data.frame(x) |>
    as_tibble() |>
    mutate(
      source = nm,
      log10_p_adjust = -log10(p.adjust)
    )
})

write.csv(
  go_all_tbl,
  file = file.path("GO_enrich_results", "GO_all_results.csv"),
  row.names = FALSE
)

go_top_tbl <- go_all_tbl |>
  group_by(Description) |>
  summarise(
    mean_log10_p_adjust = mean(log10_p_adjust, na.rm = TRUE),
    mean_p_adjust = mean(p.adjust, na.rm = TRUE),
    n_sources = n(),
    .groups = "drop"
  ) |>
  arrange(desc(mean_log10_p_adjust)) |>
  slice_head(n = 20) |>
  mutate(
    Description = forcats::fct_reorder(Description, mean_log10_p_adjust)
  )

#attēls ar vidējām p.adjust vērtībām
go_plot_p.adj <- ggplot(
  go_top_tbl,
  aes(
    x = mean_log10_p_adjust,
    y = forcats::fct_reorder(Description, mean_log10_p_adjust),
    fill = mean_log10_p_adjust
  )
) +
  geom_col(width = 0.75, color = "grey20", linewidth = 0.2) +
  scale_fill_gradient(
    low = "#CFE8F3",
    high = "#0B6E8A",
    guide = "none"
  ) +
  labs(
    x = expression("Mean " * -log[10] * "(p.adjust)"),
    y = "GO termini (Bioloģiskie procesi)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.y = element_text(size = 10, color = "black"),
    axis.text.x = element_text(color = "black"),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "grey85"),
    plot.margin = margin(10, 15, 10, 10)
  )

ggsave("Figure_9.png", plot = go_plot_p.adj, width = 9, height = 7, dpi = 300)

GO_plot_15 <- dotplot(GO_sig_res[[29]], showCategory = 10) + #15.Pacienta otrs lielākais klasteris Nr.6
  scale_y_discrete(labels = \(x) stringr::str_wrap(x, width = 35)) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    plot.margin = margin(5.5, 5.5, 5.5, 20)
  )

common_terms <- go_all_tbl |>
  dplyr::distinct(source, ID, Description) |>
  dplyr::count(ID, Description, sort = TRUE)

common_terms_2plus <- go_all_tbl |>
  dplyr::distinct(source, ID, Description) |>
  dplyr::count(ID, Description, sort = TRUE) |>
  dplyr::filter(n >= 2)

#TOP 3 kopīgie GO terms
top_terms <- common_terms |>
  dplyr::slice_head(n = 3)

#Grafiks ar biežāk sastopamajiem GO terminiem
GO_plot_freq <- common_terms |>
  dplyr::slice_head(n = 20) |>
  dplyr::mutate(
    Description = stringr::str_wrap(Description, width = 50)
  ) |>
  ggplot(aes(x = n, y = reorder(Description, n))) +
  geom_col(
    width = 0.75,
    fill = "#0B6E8A",
    color = "grey20",
    linewidth = 0.2
  ) +
  labs(
    x = "Biežums",
    y = "GO termini (Bioloģiskie procesi)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.y = element_text(size = 10, color = "black"),
    axis.text.x = element_text(color = "black"),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "grey85"),
    plot.margin = margin(10, 15, 10, 10)
  )

ggsave(filename = "Figure_7.png",plot = GO_plot_freq,width = 9,height = 7,dpi = 300)

#Sadala pa gēniem
go_genes_long <- go_all_tbl |>
  dplyr::select(source, ID, Description, geneID, Count, p.adjust) |>
  tidyr::separate_rows(geneID, sep = "/") |>
  dplyr::rename(gene = geneID)

head(go_genes_long)

common_genes <- go_genes_long |>
  dplyr::distinct(source, gene) |>
  dplyr::count(gene, sort = TRUE)

head(common_genes, 20)

common_genes_in_top_terms <- go_genes_long |>
  dplyr::distinct(source, ID, Description, gene) |>
  dplyr::count(ID, Description, gene, sort = TRUE)

#Hub Virsotnes pēc degree un TOP30
BALL_network_hubs <- lapply(BALL_graphs, degree)
BALL_network_top30 <- lapply(BALL_network_hubs, function(x) head(sort(x, decreasing = TRUE), 30))
BALL_network_top30_df <- do.call(rbind, lapply(names(BALL_network_top30), function(patient) {
  data.frame(
    patient = patient,
    gene = names(BALL_network_top30[[patient]]),
    degree = as.numeric(BALL_network_top30[[patient]])
  )
}))

#HUB virsotnes top2 klasteriem un TOP20
BALL_top2_hubs <- mapply(function(graph, top2_clusters) {
  lapply(top2_clusters, function(cluster_genes) {
    subg <- igraph::induced_subgraph(graph, vids = cluster_genes)
    igraph::degree(subg)
  })
}, BALL_graphs, BALL_top2_cl, SIMPLIFY = FALSE)

BALL_top2_top20 <- lapply(BALL_top2_hubs, function(sample_clusters) {
  lapply(sample_clusters, function(deg_vec) {
    head(sort(deg_vec, decreasing = TRUE), 20)
  })
})

#Kopīgās HUB virsotnes tīkliem starp top30 virsotnēm
hub_freq_networks <- sort(
  table(unlist(lapply(BALL_network_top30, names))),
  decreasing = TRUE
)

hub_freq_networks_df <- data.frame(
  gene = names(hub_freq_networks),
  n_patients = as.integer(hub_freq_networks),
  row.names = NULL
)

patient_ids <- names(BALL_network_top30)
if (is.null(patient_ids) || any(patient_ids == "")) {
  patient_ids <- paste0("patient_", seq_along(BALL_network_top30))
}

hub_patient_membership <- lapply(names(hub_freq_networks), function(gene) {
  present_in <- patient_ids[
    vapply(
      BALL_network_top30,
      function(x) gene %in% names(x),
      logical(1)
    )
  ]
  
  data.frame(
    gene = gene,
    n_patients = length(present_in),
    patients = paste(present_in, collapse = ", "),
    row.names = NULL
  )
})

hub_patient_membership <- do.call(rbind, hub_patient_membership)

hub_patient_membership
write.csv(hub_patient_membership, "hub_patient_membership.csv", row.names = FALSE)

#Atšķirīgās HUB virsotnes
gene_freq <- table(unlist(lapply(BALL_network_top30, names)))
unique_genes <- gene_freq[gene_freq == 1]

unique_genes
length(unique_genes)

#Kopīgās HUB virsotnes starp lielākajiem klasteriem:
#starp 1.klasteriem

BALL_cluster1_top20 <- lapply(BALL_top2_top20, function(x) x[[1]])
hub_freq_cluster1 <- sort(
  table(unlist(lapply(BALL_cluster1_top20, names))),
  decreasing = TRUE
)
hub_freq_cluster1_df <- data.frame(
  gene = names(hub_freq_cluster1),
  n_patients = as.integer(hub_freq_cluster1),
  row.names = NULL
)

BALL_cluster2_top20 <- lapply(BALL_top2_top20, function(x) x[[2]])
patient_ids <- names(BALL_cluster1_top20)

if (is.null(patient_ids) || any(patient_ids == "")) {
  patient_ids <- paste0("patient_", seq_along(BALL_cluster1_top20))
}

hub_patient_membership_cluster1 <- lapply(names(hub_freq_cluster1), function(gene) {
  present_in <- patient_ids[
    vapply(
      BALL_cluster1_top20,
      function(x) gene %in% names(x),
      logical(1)
    )
  ]
  
  data.frame(
    gene = gene,
    n_patients = length(present_in),
    patients = paste(present_in, collapse = ", "),
    row.names = NULL
  )
})
hub_patient_membership_cluster1 <- do.call(rbind, hub_patient_membership_cluster1)

#starp 2.klasteriem
BALL_cluster2_top20 <- lapply(BALL_top2_top20, function(x) x[[2]])

hub_freq_cluster2 <- sort(
  table(unlist(lapply(BALL_cluster2_top20, names))),
  decreasing = TRUE
)

hub_freq_cluster2_df <- data.frame(
  gene = names(hub_freq_cluster2),
  n_patients = as.integer(hub_freq_cluster2),
  row.names = NULL
)

patient_ids <- names(BALL_cluster2_top20)

if (is.null(patient_ids) || any(patient_ids == "")) {
  patient_ids <- paste0("patient_", seq_along(BALL_cluster2_top20))
}

hub_patient_membership_cluster2 <- lapply(names(hub_freq_cluster2), function(gene) {
  present_in <- patient_ids[
    vapply(
      BALL_cluster2_top20,
      function(x) gene %in% names(x),
      logical(1)
    )
  ]
  
  data.frame(
    gene = gene,
    n_patients = length(present_in),
    patients = paste(present_in, collapse = ", "),
    row.names = NULL
  )
})

hub_patient_membership_cluster2 <- do.call(rbind, hub_patient_membership_cluster2)

#Mutāciju tabulas: SNV Germline, SNV somatic, SV germline, SV somatic
SNV_germline <- read.delim("HM_RakstaKohorta_Collected_VEP_SNVgermline.tsv")
SNV_somatic <- read.delim("HM_RakstaKohorta_Collected_VEP_SNVsomatic.tsv")
SV_germline <- read.delim("HM_RakstaKohorta_Collected_VEP_SVgermline.tsv")
SV_somatic <- read.delim("HM_RakstaKohorta_Collected_VEP_SVsomatic.tsv")

hub_genes_by_patient <- lapply(BALL_network_top30, names)

#SNV germline TOP30 hub virsotnēm
SNV_germline_clean <- SNV_germline |>
  dplyr::select(RS, SYMBOL) |>
  dplyr::filter(!is.na(RS), !is.na(SYMBOL)) |>
  dplyr::distinct()

snv_germline_genes_by_patient <- split(
  SNV_germline_clean$SYMBOL,
  SNV_germline_clean$RS
)

common_patients_snv_germline <- intersect(
  names(hub_genes_by_patient),
  names(snv_germline_genes_by_patient)
)

snv_germline_hub_overlap <- lapply(common_patients_snv_germline, function(patient_id) {
  intersect(
    hub_genes_by_patient[[patient_id]],
    snv_germline_genes_by_patient[[patient_id]]
  )
})

names(snv_germline_hub_overlap) <- common_patients_snv_germline

#SNV somatic TOP30 hub virsotnēm
SNV_somatic_clean <- SNV_somatic |>
  dplyr::select(RS, SYMBOL) |>
  dplyr::filter(!is.na(RS), !is.na(SYMBOL)) |>
  dplyr::distinct()

snv_somatic_genes_by_patient <- split(
  SNV_somatic_clean$SYMBOL,
  SNV_somatic_clean$RS
)

common_patients_snv_somatic <- intersect(
  names(hub_genes_by_patient),
  names(snv_somatic_genes_by_patient)
)

snv_somatic_hub_overlap <- lapply(common_patients_snv_somatic, function(patient_id) {
  intersect(
    hub_genes_by_patient[[patient_id]],
    snv_somatic_genes_by_patient[[patient_id]]
  )
})

names(snv_somatic_hub_overlap) <- common_patients_snv_somatic

#SV germline TOP30 hub virsotnēm
SV_germline_clean <- SV_germline |>
  dplyr::select(RS, SYMBOL) |>
  dplyr::filter(!is.na(RS), !is.na(SYMBOL)) |>
  dplyr::distinct()

sv_germline_genes_by_patient <- split(
  SV_germline_clean$SYMBOL,
  SV_germline_clean$RS
)

common_patients_sv_germline <- intersect(
  names(hub_genes_by_patient),
  names(sv_germline_genes_by_patient)
)

sv_germline_hub_overlap <- lapply(common_patients_sv_germline, function(patient_id) {
  intersect(
    hub_genes_by_patient[[patient_id]],
    sv_germline_genes_by_patient[[patient_id]]
  )
})

names(sv_germline_hub_overlap) <- common_patients_sv_germline

#SV somatic TOP30 hub virsotnēm
SV_somatic_clean <- SV_somatic |>
  dplyr::select(RS, SYMBOL) |>
  dplyr::filter(!is.na(RS), !is.na(SYMBOL)) |>
  dplyr::distinct()

sv_somatic_genes_by_patient <- split(
  SV_somatic_clean$SYMBOL,
  SV_somatic_clean$RS
)

common_patients_sv_somatic <- intersect(
  names(hub_genes_by_patient),
  names(sv_somatic_genes_by_patient)
)

sv_somatic_hub_overlap <- lapply(common_patients_sv_somatic, function(patient_id) {
  intersect(
    hub_genes_by_patient[[patient_id]],
    sv_somatic_genes_by_patient[[patient_id]]
  )
})

names(sv_somatic_hub_overlap) <- common_patients_sv_somatic

#hub virsotņu degree ar sv somatic
sv_somatic_hub_degree <- lapply(names(sv_somatic_hub_overlap), function(patient_id) {
  
  genes <- sv_somatic_hub_overlap[[patient_id]]
  
  if (length(genes) == 0) return(NULL)
  
  deg_vec <- BALL_network_top30[[patient_id]][genes]
  
  data.frame(
    patient = patient_id,
    gene = names(deg_vec),
    degree = as.numeric(deg_vec),
    row.names = NULL
  )
})

sv_somatic_hub_degree <- do.call(rbind, sv_somatic_hub_degree)

#Mutācijas tīklu līmenī
mutation_tables <- list(
  SNV_germline = SNV_germline,
  SNV_somatic = SNV_somatic,
  SV_germline = SV_germline,
  SV_somatic = SV_somatic
)
get_mutated_nodes_degree <- function(mutation_table, graphs_list, mutation_type) {
  
  mutation_clean <- mutation_table |>
    filter(!is.na(RS), !is.na(SYMBOL)) |>
    distinct() |>
    mutate(RS = as.character(RS))
  
  mutation_genes_by_patient <- split(
    mutation_clean$SYMBOL,
    mutation_clean$RS
  )
  
  common_patients <- intersect(
    names(graphs_list),
    names(mutation_genes_by_patient)
  )
  
  map_dfr(common_patients, \(patient_id) {
    graph <- graphs_list[[patient_id]]
    mutated_genes <- mutation_genes_by_patient[[patient_id]]
    
    mutated_nodes <- intersect(V(graph)$name, mutated_genes)
    
    if (length(mutated_nodes) == 0) return(NULL)
    
    tibble(
      patient_id = patient_id,
      mutation_type = mutation_type,
      gene = mutated_nodes,
      degree = degree(graph, v = mutated_nodes)
    ) |>
      arrange(desc(degree))
  })
}
mutated_nodes_degree_results <- imap(
   mutation_tables,
   \(mutation_table, mutation_type) {
    get_mutated_nodes_degree(
    mutation_table,
    BALL_graphs,
    mutation_type)})

all_mutation_results <- bind_rows(mutated_nodes_degree_results)
all_mutation_results <- all_mutation_results |>
  dplyr::arrange(desc(degree))

View(all_mutation_results)
all_mutation_results_anon <- all_mutation_results %>%
  left_join(id_key, by = c("patient_id" = "sample")) %>%
  mutate(patient_id = anon_id) %>%
  select(-anon_id)

all_mutation_results_anon_clean <- all_mutation_results_anon %>%
  distinct(gene, patient_id, mutation_type)

library(ggplot2)

all_mutation_results_anon_clean_lv <- all_mutation_results_anon_clean %>%
  mutate(mutation_type = recode(mutation_type,
                                "SNV_somatic"  = "Somatisks SNV",
                                "SV_somatic"   = "Somatisks SV",
                                "SNV_germline" = "Pārmantots SNV",
                                "SV_germline"  = "Pārmantots SV"
  ))
# Aprēķina kopējo mutāciju skaitu katram pacientam
mutation_counts <- all_mutation_results_anon_clean_lv %>%
  group_by(patient_id) %>%
  summarise(total = n())

ggplot(all_mutation_results_anon_clean_lv, 
       aes(x = patient_id, fill = mutation_type)) +
  geom_bar(width = 0.7) +
  geom_text(data = mutation_counts,
            aes(x = patient_id, y = total, label = total),
            inherit.aes = FALSE,
            vjust = -0.5, size = 3.5, fontface = "bold") +
  scale_fill_manual(values = c(
    "Somatisks SNV" = "#08306B",  
    "Somatisks SV"  = "#2171B5",  
    "Pārmantots SNV"  = "#6BAED6",  
    "Pārmantots SV"   = "#C6DBEF"   
  )) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    axis.text.y = element_text(size = 9),
    axis.title = element_text(size = 11, face = "bold"),
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    panel.grid.major.x = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12)
  ) +
  labs(x = "Pacients", 
       y = "Mutāciju skaits", 
       fill = "Mutācijas veids",
       )

###
seed_after <- .Random.seed

