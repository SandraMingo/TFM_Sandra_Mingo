#!/usr/bin/env Rscript
# ======================================================================
# REDES DE CO-OCURRENCIA MICROBIANA (PRESENCIA/AUSENCIA Y FASTSPAR)
# ----------------------------------------------------------------------
# Descripción:
#   1. Carga tabla QIIME2 y metadata limpia.
#   2. Construye matriz de abundancias relativas y selecciona taxones
#      más abundantes para análisis de co-ocurrencia (CooccurrenceAffinity).
#   3. Genera matrices de presencia/ausencia, afinidad (alpha_mle) y
#      redes de co-ocurrencia con clustering.
#   4. Prepara matriz de abundancia colapsada a nivel de género para
#      FastSpar/SparCC.
#   5. Procesa matrices de correlación y p-valores de FastSpar para
#      generar una red de co-ocurrencia basada en correlaciones.
#
# Entradas:
#   - 1_input/aggregated_taxonomy_table_qiime.tsv
#   - 0_tables/metadata_clean.tsv
#   - 2_fastspar/correlation.tsv
#   - 2_fastspar/pvalues.tsv
#
# Salidas:
#   - 3_cooccurrence_network/cooccur_OTU_species_aggregated_top.tsv
#   - 3_cooccurrence_network/cooccur_top_PA.tsv
#   - 3_cooccurrence_network/cooccur_affinity.tsv
#   - 3_cooccurrence_network/cooccur_alpha_mle.tsv
#   - 3_cooccurrence_network/cooccur_pvalues.tsv
#   - 3_cooccurrence_network/cooccur_alpha_sign_selected.tsv
#   - 3_cooccurrence_network/net.sp.txt
#   - 0_tables/abundance_matrix.tsv (para FastSpar)
#   - 0_figs/cooccur_alpha_mle.svg
#   - 0_figs/cooccur_network_sign.svg
#   - 0_figs/cooccur_clustering_network.svg
#   - 0_figs/cooccur_network_sign_clusterized.svg
#   - 0_figs/correlation_network.svg
#
# Requisitos:
#   - R >= 4.0
#   - Paquetes: CooccurrenceAffinity, igraph, wesanderson, data.table,
#               dplyr, RColorBrewer, tidyverse, ggraph, tidygraph, tidyr
#
# Autor: Sandra Mingo-Ramirez
# Fecha: 2025
# ======================================================================

rm(list = ls())
set.seed(123)

PROJECT_DIR <- "/home/sandra/Projects/UCI_Variables/"
setwd(PROJECT_DIR)

dir.create("0_figs/", recursive = TRUE, showWarnings = FALSE)
dir.create("0_tables/", recursive = TRUE, showWarnings = FALSE)
dir.create("3_cooccurrence_network/", recursive = TRUE, showWarnings = FALSE)

# ----------------------------------------------------------------------
# 1. LIBRERÍAS
# ----------------------------------------------------------------------
library(CooccurrenceAffinity)
library(igraph)
library(wesanderson)
library(data.table)
library(dplyr)
library(RColorBrewer)
library(tidyverse)
library(ggraph)
library(tidygraph)
library(tidyr)

# ----------------------------------------------------------------------
# 2. CARGA Y SINCRONIZACIÓN DE DATOS
# ----------------------------------------------------------------------
metadata <- read.table(
  "0_tables/metadata_clean.tsv",
  header      = TRUE,
  row.names   = 1,
  sep         = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

table_qiime <- fread(
  "1_input/aggregated_taxonomy_table_qiime.tsv",
  header         = TRUE,
  data.table     = FALSE,
  fill           = TRUE,
  stringsAsFactors = FALSE
)

# Limpieza de nombres
rownames(table_qiime) <- gsub("^[a-z]__", "", table_qiime[, 1])
table_qiime <- table_qiime[, -1, drop = FALSE]
colnames(table_qiime) <- gsub("_L001", "", colnames(table_qiime))

# Sincronizar muestras
samples <- intersect(colnames(table_qiime), rownames(metadata))
data_tab <- table_qiime[, samples, drop = FALSE]
metadata <- metadata[samples, , drop = FALSE]

# ----------------------------------------------------------------------
# 3. CO-OCCURRENCEAFFINITY: MATRICES Y RED DE PRESENCIA/AUSENCIA
# ----------------------------------------------------------------------
# 3.1 Transformar a abundancias relativas (%)
data_tab_prop <- prop.table(as.matrix(data_tab), 2) * 100

# 3.2 Seleccionar taxones más abundantes
data_mean <- rowMeans(data_tab_prop)
cca.data  <- cbind(data_tab_prop, data_mean)
cca.data  <- cca.data[order(cca.data[, "data_mean"], decreasing = TRUE), ]

# Número de taxones top
n_top <- 25
topn  <- cca.data[1:n_top, -ncol(cca.data), drop = FALSE]

# Eliminar duplicados de nombre si existieran
topn <- topn[!duplicated(rownames(topn)), , drop = FALSE]

# Crear tabla taxonómica para topn
tax_strings <- rownames(topn)   # k__...|p__...|...

tax_split <- strsplit(tax_strings, "\\|")

tax_df <- do.call(rbind, lapply(tax_split, function(x) {
  x <- c(x, rep(NA, 7 - length(x)))  # hasta 7 niveles
  x
}))

colnames(tax_df) <- c("Kingdom","Phylum","Class","Order",
                      "Family","Genus","Species")
rownames(tax_df) <- tax_strings

# Simplificar nombre al último nivel taxonómico (sin prefijos)
simple_names <- sub(
  "^[a-z]__*", "",
  sapply(strsplit(tax_strings, "\\|"), tail, 1)
)

tax_table_topn <- data.frame(
  Taxon_simple = simple_names,
  Kingdom = tax_df[, "Kingdom"],
  Phylum  = tax_df[, "Phylum"],
  Class   = tax_df[, "Class"],
  Order   = tax_df[, "Order"],
  Family  = tax_df[, "Family"],
  Genus   = tax_df[, "Genus"],
  Species = tax_df[, "Species"],
  row.names = simple_names,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

rownames(topn) <- simple_names

write.table(
  topn,
  "3_cooccurrence_network/cooccur_OTU_species_aggregated_top.tsv",
  sep       = "\t",
  quote     = FALSE,
  col.names = FALSE
)

# 3.3 Matriz de presencia/ausencia
pa_data <- (topn > 0) * 1

write.table(
  pa_data,
  "3_cooccurrence_network/cooccur_top_PA.tsv",
  sep       = "\t",
  quote     = FALSE,
  col.names = NA
)

# 3.4 Análisis de afinidad de co-ocurrencia
myout <- affinity(
  data         = pa_data,
  row.or.col   = "row",
  squarematrix = "all"
)

write.table(
  myout$all,
  "3_cooccurrence_network/cooccur_affinity.tsv",
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE
)

alpha_mat <- myout$alpha_mle
pval_mat  <- myout$pval

write.table(
  alpha_mat,
  "3_cooccurrence_network/cooccur_alpha_mle.tsv",
  sep       = "\t",
  quote     = FALSE
)

write.table(
  pval_mat,
  "3_cooccurrence_network/cooccur_pvalues.tsv",
  sep       = "\t",
  quote     = FALSE
)

# 3.5 Selección de interacciones significativas
corr.set <- myout$alpha_mle_sig
write.table(
  corr.set,
  "3_cooccurrence_network/cooccur_alpha_sign_selected.tsv",
  sep       = "\t",
  quote     = FALSE,
  col.names = NA
)

# Matriz de adyacencia
corr.set <- as.matrix(corr.set)
corr.set[is.na(corr.set)] <- 0
corr.mat <- corr.set[, 1:ncol(corr.set), drop = FALSE]

# 3.6 Plot por defecto de CooccurrenceAffinity
plotgg(
  data        = myout,
  variable    = "alpha_mle",
  legendlimit = "datarange",
  show.value  = TRUE,
  text.size   = 2.5
)
ggsave(
  "0_figs/cooccur_alpha_mle.svg",
  width  = 10,
  height = 10
)

# 3.7 Red de co-ocurrencia (affinity)
network <- graph_from_adjacency_matrix(
  corr.mat,
  weighted = TRUE,
  mode     = "undirected",
  diag     = FALSE
)

write_graph(
  network,
  "3_cooccurrence_network/net.sp.txt",
  format = "edgelist"
)

# Eliminar nodos aislados
network <- delete_vertices(network, which(degree(network) == 0))

taxa    <- V(network)$name
phylum_fac <- factor(tax_table_topn[taxa, "Phylum"])
n_phyla <- length(levels(phylum_fac))
pal_phy <- brewer.pal(min(n_phyla, 12), "Set1")
V(network)$color <- pal_phy[ as.numeric(phylum_fac) ]

edge_vals        <- as.numeric(E(network)$weight)
E(network)$color <- ifelse(edge_vals > 0, "blue", "red")

svg(
  filename = "0_figs/cooccur_network_sign.svg",
  width    = 8,
  height   = 8
)
plot(
  network,
  vertex.label.family = "Helvetica",
  vertex.frame.color  = "transparent",
  vertex.label.color  = "black",
  vertex.color        = V(network)$color,
  edge.color          = E(network)$color,
  vertex.label.cex    = 0.5,
  edge.arrow.size     = 0.1
)
legend(
  "topleft",
  legend = levels(phylum_fac),
  col    = pal_phy[seq_along(levels(phylum_fac))],
  pch    = 16,
  pt.cex = 1.5,
  bty    = "n",
  title  = "Phylum"
)
dev.off()

# 3.8 Clustering de la red
ceb <- cluster_edge_betweenness(network)

svg(
  filename = "0_figs/cooccur_clustering_network.svg",
  width    = 8,
  height   = 8
)
plot_dendrogram(ceb, mode = "hclust")
dev.off()

svg(
  filename = "0_figs/cooccur_network_sign_clusterized.svg",
  width    = 8,
  height   = 8
)
plot(
  ceb,
  network,
  vertex.label.family = "Helvetica",
  vertex.frame.color  = "transparent",
  vertex.label.color  = "black",
  vertex.color        = V(network)$color,
  edge.color          = E(network)$color,
  vertex.label.cex    = 0.5,
  edge.arrow.size     = 0.1
)
dev.off()

# ----------------------------------------------------------------------
# 4. FASTSPAR: MATRIZ DE ABUNDANCIA A NIVEL GÉNERO
# ----------------------------------------------------------------------
# 4.1 Volver a cargar tabla QIIME (sin normalizar) para colapsar a género
table_qiime2 <- fread(
  "1_input/aggregated_taxonomy_table_qiime.tsv",
  header     = TRUE,
  data.table = FALSE,
  fill       = TRUE
)
rownames(table_qiime2) <- table_qiime2[, 1]
table_qiime2 <- table_qiime2[, -1, drop = FALSE]
colnames(table_qiime2) <- gsub("_L001$", "", colnames(table_qiime2))

tax_strings <- rownames(table_qiime2)

tax_clean <- tax_strings %>%
  gsub("d__", "", .) %>%
  gsub("p__", "", .) %>%
  gsub("c__", "", .) %>%
  gsub("o__", "", .) %>%
  gsub("f__", "", .) %>%
  gsub("g__", "", .)

tax_split <- strsplit(tax_clean, "\\|")

max_rank <- 6
tax_mat  <- t(vapply(
  tax_split,
  function(x) {
    length(x) <- max_rank
    x
  },
  character(max_rank)
))

colnames(tax_mat) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus")
tax_df <- as.data.frame(tax_mat, stringsAsFactors = FALSE)
rownames(tax_df) <- tax_strings

keep_genus <- !is.na(tax_df$Genus) &
  tax_df$Genus != "" &
  tax_df$Genus != "Unassigned"

otu_genus_level <- table_qiime2[keep_genus, , drop = FALSE]
tax_genus_level <- tax_df[keep_genus, , drop = FALSE]

otu_genus <- rowsum(
  otu_genus_level,
  group = tax_genus_level$Genus
)

otu_biom <- cbind(OTU_ID = rownames(otu_genus), otu_genus)
colnames(otu_biom)[1] <- "#OTU ID"

write.table(
  otu_biom,
  file      = "0_tables/abundance_matrix.tsv",
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE
)

# Ejecutar FastSpar

# ----------------------------------------------------------------------
# 5. PROCESAR RESULTADOS FASTSPAR Y GENERAR RED DE CORRELACIÓN
# ----------------------------------------------------------------------
# (Suponiendo que FastSpar ya se ejecutó y generó correlation.tsv y pvalues.tsv)

cor_mat <- as.matrix(read.table(
  "2_fastspar/correlation.tsv",
  header      = TRUE,
  row.names   = 1,
  check.names = FALSE,
  sep         = "\t",
  comment.char = ""
))

p_mat <- as.matrix(read.table(
  "2_fastspar/pvalues.tsv",
  header      = TRUE,
  row.names   = 1,
  check.names = FALSE,
  sep         = "\t",
  comment.char = ""
))

stopifnot(identical(rownames(cor_mat), rownames(p_mat)))
stopifnot(identical(colnames(cor_mat), colnames(p_mat)))

# 5.1 Crear tabla de aristas
r_cutoff <- 0.5
p_cutoff <- 0.05

cor_df <- cor_mat %>%
  as.data.frame() %>%
  rownames_to_column("Taxon1") %>%
  pivot_longer(
    cols      = -Taxon1,
    names_to  = "Taxon2",
    values_to = "r"
  )

p_df <- p_mat %>%
  as.data.frame() %>%
  rownames_to_column("Taxon1") %>%
  pivot_longer(
    cols      = -Taxon1,
    names_to  = "Taxon2",
    values_to = "p"
  )

edges <- cor_df %>%
  left_join(p_df, by = c("Taxon1", "Taxon2")) %>%
  filter(Taxon1 < Taxon2) %>%
  filter(abs(r) >= r_cutoff, p <= p_cutoff) %>%
  mutate(
    Taxon1 = as.character(Taxon1),
    Taxon2 = as.character(Taxon2)
  )

# 5.2 Tabla de nodos con taxonomía colapsada
tax_genus <- tax_genus_level %>%
  group_by(Genus) %>%
  summarize(
    Phylum = first(Phylum),
    Class  = first(Class),
    Order  = first(Order),
    Family = first(Family),
    .groups = "drop"
  )

nodes <- tibble(name = unique(c(edges$Taxon1, edges$Taxon2))) %>%
  left_join(tax_genus, by = c("name" = "Genus"))

nodes$Phylum[is.na(nodes$Phylum)] <- "Unknown"

# 5.3 Construcción del grafo
g <- graph_from_data_frame(
  d        = edges,
  vertices = nodes,
  directed = FALSE
)

E(g)$sign   <- ifelse(E(g)$r > 0, "positive", "negative")
E(g)$weight <- abs(E(g)$r)
V(g)$degree <- degree(g)

g_tbl <- as_tbl_graph(g)

# 5.4 Grafo de correlaciones
p_net <- ggraph(g_tbl, layout = "fr") +
  geom_edge_link(
    aes(color = sign, width = weight),
    alpha = 0.4
  ) +
  scale_edge_color_manual(
    values = c(positive = "firebrick3", negative = "steelblue3"),
    name   = "Correlación"
  ) +
  scale_edge_width(range = c(0.2, 2), guide = "none") +
  geom_node_point(
    aes(color = Phylum, size = degree),
    alpha = 0.9
  ) +
  scale_size_continuous(range = c(2, 8), guide = "none") +
  geom_node_text(
    aes(label = name),
    size  = 2,
    repel = TRUE
  ) +
  theme_graph() +
  labs(title = "Red de correlación con FastSpar")

ggsave(
  "0_figs/correlation_network.svg",
  p_net,
  width  = 8,
  height = 8
)
