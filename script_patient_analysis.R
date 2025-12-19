#!/usr/bin/env Rscript
# ======================================================================
# ANÁLISIS DE DIVERSIDAD Y COMPOSICIÓN POR PACIENTE (QIIME2)
# ----------------------------------------------------------------------
# Descripción:
#   1. Limpia y sincroniza tabla QIIME2 con metadata, generando metadata_clean.
#   2. Calcula métricas de diversidad alfa (Shannon, Chao1, dbp) por paciente.
#   3. Genera análisis de diversidad beta (Bray-Curtis, Jaccard) con NMDS por paciente.
#   4. Crea composiciones relativas por paciente y tiempo (top 10 taxones).
#   5. Visualiza evolución temporal de diversidad alfa por paciente (individuos y combinado).
#   6. Genera histogramas de distribución de muestras por paciente.
#
# Entradas:
#   - 1_input/aggregated_taxonomy_table_qiime.tsv
#   - 1_input/metadata.csv
#
# Salidas:
#   - 0_tables/metadata_clean.tsv
#   - 0_tables/alpha_metrics_output.tsv
#   - 0_figs/Paciente_NHC/muestras_por_paciente.svg
#   - 0_figs/Paciente_NHC/distribucion_muestras.svg
#   - 0_figs/Paciente_NHC/shannon_nhc.svg, chao1_nhc.svg, dbp_nhc.svg, alfa_nhc.svg
#   - 0_figs/Paciente_NHC/bray-curtis-nmds.svg, jaccard-nmds.svg
#   - 0_figs/Paciente_NHC/barplots_relab/abund_rel_*.svg
#   - 0_figs/Paciente_NHC/shannon_evolution/shannon_evol_*.svg (y similares para chao1, dbp)
#   - 0_figs/Paciente_NHC/evolucion_metricas/paciente_*.svg
#
# Requisitos:
#   - R >= 4.0
#   - Paquetes: dplyr, ggplot2, RColorBrewer, tidyr, stringr, tidyverse,
#               broom, rstatix, vegan, microbiome, phyloseq, ggrepel
#
# Autor: Sandra Mingo-Ramirez
# Fecha: 2025
# ======================================================================

set.seed(123)

PROJECT_DIR <- "/home/sandra/Projects/UCI_Variables/"
setwd(PROJECT_DIR)

# ----------------------------------------------------------------------
# 1. LIBRERÍAS
# ----------------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(tidyr)
library(stringr)
library(tidyverse)
library(broom)
library(rstatix)
library(vegan)
library(microbiome)
library(phyloseq)
library(ggrepel)

dir.create("0_figs/Paciente_NHC", recursive = TRUE, showWarnings = FALSE)
dir.create("0_tables", recursive = TRUE, showWarnings = FALSE)
dir.create("0_figs/Paciente_NHC/barplots_relab", recursive = TRUE, showWarnings = FALSE)
dir.create("0_figs/Paciente_NHC/shannon_evolution", recursive = TRUE, showWarnings = FALSE)
dir.create("0_figs/Paciente_NHC/chao1_evolution", recursive = TRUE, showWarnings = FALSE)
dir.create("0_figs/Paciente_NHC/dbp_evolution", recursive = TRUE, showWarnings = FALSE)
dir.create("0_figs/Paciente_NHC/evolucion_metricas", recursive = TRUE, showWarnings = FALSE)

# ----------------------------------------------------------------------
# 2. CARGA Y PREPROCESAMIENTO DE DATOS
# ----------------------------------------------------------------------
table_qiime <- fread(
  "1_input/aggregated_taxonomy_table_qiime.tsv",
  header     = TRUE,
  data.table = FALSE,
  fill       = TRUE
)
rownames(table_qiime) <- table_qiime[, 1]
table_qiime <- table_qiime[, -1, drop = FALSE]

# Limpieza taxonómica
rownames(table_qiime) <- gsub("^d__", "", rownames(table_qiime))
rownames(table_qiime) <- gsub("p__",  "", rownames(table_qiime))
rownames(table_qiime) <- gsub("c__",  "", rownames(table_qiime))
rownames(table_qiime) <- gsub("o__",  "", rownames(table_qiime))
rownames(table_qiime) <- gsub("f__",  "", rownames(table_qiime))
rownames(table_qiime) <- gsub("g__",  "", rownames(table_qiime))
colnames(table_qiime) <- gsub("_L001", "", colnames(table_qiime))

# Carga y limpieza metadata
metadata <- read.csv(
  "1_input/metadata.csv",
  header      = TRUE,
  row.names   = 1,
  check.names = FALSE
)
rownames(metadata) <- gsub("_L001", "", rownames(metadata))

# Corrección específica NHC 314681
metadata <- metadata %>%
  mutate(
    Tiempo = ifelse(
      NHC == "314681" & Tiempo == 0,
      paste0(Tiempo, "_", ave(Tiempo, NHC, Tiempo, FUN = seq_along)),
      Tiempo
    )
  ) %>%
  filter(!(rownames(metadata) %in% c("179_S93", "V1-V2-S-86-2024-9_S23")))

metadata$SampleID <- rownames(metadata)
metadata <- metadata %>% relocate(SampleID, .before = 1)

# Guardar metadata limpia
write.table(
  metadata,
  file       = "0_tables/metadata_clean.tsv",
  sep        = "\t",
  quote      = FALSE,
  row.names  = FALSE
)

# Sincronizar muestras
samples <- intersect(colnames(table_qiime), rownames(metadata))
table_qiime <- table_qiime[, samples, drop = FALSE]
metadata    <- metadata[samples, , drop = FALSE]

# ----------------------------------------------------------------------
# 3. DISTRIBUCIÓN DE MUESTRAS POR PACIENTE
# ----------------------------------------------------------------------
num_muestras <- table(metadata$NHC)
df_muestras  <- as.data.frame(num_muestras)
colnames(df_muestras) <- c("NHC", "Num_muestras")

# Histograma muestras por paciente
ggplot(df_muestras, aes(x = NHC, y = Num_muestras)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  scale_y_continuous(breaks = 0:max(df_muestras$Num_muestras)) +
  theme(axis.text.x = element_blank()) +
  labs(
    title = "Número de muestras por paciente",
    x     = "Paciente (NHC)",
    y     = "Número de muestras"
  ) +
  theme_minimal()
ggsave("0_figs/Paciente_NHC/muestras_por_paciente.svg", width = 5, height = 4)

# Distribución frecuencia
df_freq <- df_muestras %>%
  count(Num_muestras)

ggplot(df_freq, aes(x = factor(Num_muestras), y = n)) +
  geom_bar(stat = "identity", fill = "skyblue", color = "black") +
  geom_text(aes(label = n), vjust = -0.5) +
  theme_minimal() +
  labs(
    title    = "Distribución del número de muestras por paciente",
    x        = "Cantidad de muestras por paciente",
    y        = "Número de pacientes"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("0_figs/Paciente_NHC/distribucion_muestras.svg", width = 5, height = 5)

# ----------------------------------------------------------------------
# 4. DIVERSIDAD ALFA POR PACIENTE
# ----------------------------------------------------------------------
data_tab <- otu_table(table_qiime, taxa_are_rows = TRUE)
meta_tab <- sample_data(metadata)
pseq     <- merge_phyloseq(data_tab, meta_tab)

alpha_metrics <- microbiome::alpha(pseq, index = c("Shannon", "Chao1", "dbp"))

ps1.meta <- metadata
ps1.meta$Shannon <- alpha_metrics$diversity_shannon
ps1.meta$chao1   <- alpha_metrics$chao1
ps1.meta$dbp     <- alpha_metrics$dominance_dbp
ps1.meta$NHC     <- as.factor(ps1.meta$NHC)

write.table(
  ps1.meta,
  "0_tables/alpha_metrics_output.tsv",
  sep       = "\t",
  row.names = FALSE
)

# Boxplots por paciente
p1 <- ggplot(ps1.meta, aes(x = NHC, y = Shannon, fill = NHC)) +
  geom_boxplot() +
  scale_fill_viridis_d(option = "D") +
  theme_minimal() +
  labs(x = "Paciente", y = "Shannon") +
  theme(axis.text.x = element_blank(), legend.position = "none")

p2 <- ggplot(ps1.meta, aes(x = NHC, y = chao1, fill = NHC)) +
  geom_boxplot() +
  scale_fill_viridis_d(option = "D") +
  theme_minimal() +
  labs(x = "Paciente", y = "Chao1") +
  theme(axis.text.x = element_blank(), legend.position = "none")

p3 <- ggplot(ps1.meta, aes(x = NHC, y = dbp, fill = NHC)) +
  geom_boxplot() +
  scale_fill_viridis_d(option = "D") +
  theme_minimal() +
  labs(x = "Paciente", y = "dbp") +
  theme(axis.text.x = element_blank(), legend.position = "none")

ggsave("0_figs/Paciente_NHC/shannon_nhc.svg", p1, width = 7, height = 5)
ggsave("0_figs/Paciente_NHC/chao1_nhc.svg",   p2, width = 7, height = 5)
ggsave("0_figs/Paciente_NHC/dbp_nhc.svg",     p3, width = 7, height = 5)
(p1 | p2 | p3)
ggsave("0_figs/Paciente_NHC/alfa_nhc.svg", width = 10, height = 5)

# ----------------------------------------------------------------------
# 5. DIVERSIDAD BETA POR PACIENTE (Bray-Curtis y Jaccard)
# ----------------------------------------------------------------------
# Bray-Curtis
dist_bc <- phyloseq::distance(pseq, method = "bray")
nmds_bc <- vegan::metaMDS(dist_bc, k = 2, trymax = 100)

nmds_bc_df <- as.data.frame(nmds_bc$points)
nmds_bc_df$Sample <- rownames(nmds_bc_df)
nmds_bc_df$NHC    <- metadata[nmds_bc_df$Sample, "NHC"]
nmds_bc_df$NHC    <- as.factor(nmds_bc_df$NHC)
nmds_bc_df        <- nmds_bc_df %>% drop_na(MDS1, MDS2)

ellipse_groups <- nmds_bc_df %>%
  group_by(NHC) %>%
  filter(n() > 3, sd(MDS1) > 0, sd(MDS2) > 0)

ggplot(nmds_bc_df, aes(x = MDS1, y = MDS2, color = NHC)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(
    data         = ellipse_groups,
    aes(group = NHC, color = NHC),
    linetype     = 1,
    size         = 0.8,
    alpha        = 0.4,
    show.legend  = FALSE
  ) +
  scale_color_viridis_d(option = "D") +
  theme_minimal(base_size = 16) +
  labs(x = "NMDS1", y = "NMDS2", color = "Paciente") +
  theme(legend.position = "none")
ggsave("0_figs/Paciente_NHC/bray-curtis-nmds.svg", width = 10, height = 10)

# Jaccard
otu_pa <- data_tab
otu_pa@.Data <- ifelse(otu_pa@.Data > 0, 1, 0)
pseq_pa <- phyloseq(otu_pa, meta_tab)

dist_jac <- phyloseq::distance(pseq_pa, method = "jaccard", binary = TRUE)
nmds_jac <- vegan::metaMDS(dist_jac, k = 2, trymax = 100)

nmds_jac_df <- as.data.frame(nmds_jac$points)
nmds_jac_df$Sample <- rownames(nmds_jac_df)
nmds_jac_df$NHC    <- metadata[nmds_jac_df$Sample, "NHC"]
nmds_jac_df$NHC    <- as.factor(nmds_jac_df$NHC)
nmds_jac_df        <- nmds_jac_df %>% drop_na(MDS1, MDS2)

ellipse_groups <- nmds_jac_df %>%
  group_by(NHC) %>%
  filter(n() > 3, sd(MDS1) > 0, sd(MDS2) > 0)

p_jac <- ggplot(nmds_jac_df, aes(x = MDS1, y = MDS2, color = NHC)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(
    data         = ellipse_groups,
    aes(group = NHC, color = NHC),
    linetype     = 1,
    size         = 0.8,
    alpha        = 0.4,
    show.legend  = FALSE
  ) +
  scale_color_viridis_d(option = "D") +
  theme_minimal(base_size = 16) +
  labs(x = "NMDS1", y = "NMDS2", color = "Paciente") +
  theme(legend.position = "none")

ggsave("0_figs/Paciente_NHC/jaccard-nmds.svg", p_jac, width = 10, height = 10)

# Estadísticos beta
beta_bc <- betadisper(dist_bc, metadata$NHC)
permutest(beta_bc)
adonis2(dist_bc ~ NHC, data = metadata)

# ----------------------------------------------------------------------
# 6. COMPOSICIÓN RELATIVA POR PACIENTE Y TIEMPO
# ----------------------------------------------------------------------
prep_sample_composition <- function(abund_table, meta, top_n = 10, other_label = "Otros") {
  
  common_samples <- intersect(colnames(abund_table), rownames(meta))
  abund_table    <- as.data.frame(abund_table[, common_samples, drop = FALSE])
  meta_use       <- meta[common_samples, , drop = FALSE]
  
  col_sums <- colSums(abund_table)
  col_sums[col_sums == 0] <- 1
  abund_rel <- sweep(abund_table, 2, col_sums, "/")
  
  taxa_ids  <- rownames(abund_rel)
  last_term <- str_extract(taxa_ids, "[^|]+$")
  last_term[is.na(last_term) | last_term == ""] <- "Unassigned"
  
  abund_rel$TaxonLast <- last_term
  
  rel_agg <- abund_rel %>%
    group_by(TaxonLast) %>%
    summarise(across(where(is.numeric), sum), .groups = "drop")
  
  long_df <- rel_agg %>%
    pivot_longer(
      cols      = -TaxonLast,
      names_to  = "Sample",
      values_to = "AbundRel"
    ) %>%
    mutate(
      AbundRel = AbundRel * 100,
      NHC      = meta_use[Sample, "NHC"],
      Tiempo   = meta_use[Sample, "Tiempo"]
    )
  
  long_df <- long_df %>%
    group_by(NHC) %>%
    mutate(
      mean_tax_pac = tapply(AbundRel, TaxonLast, mean)[TaxonLast],
      Taxon = if_else(
        TaxonLast %in%
          (names(sort(tapply(AbundRel, TaxonLast, mean), decreasing = TRUE))[1:top_n]),
        TaxonLast,
        other_label
      )
    ) %>%
    group_by(Sample, NHC, Tiempo, Taxon) %>%
    summarise(AbundRel = sum(AbundRel), .groups = "drop")
  
  lev <- sort(unique(long_df$Taxon))
  if (other_label %in% lev) {
    lev <- c(setdiff(lev, other_label), other_label)
  }
  long_df$Taxon <- factor(long_df$Taxon, levels = lev)
  
  long_df
}

df_qiime_samples <- prep_sample_composition(table_qiime, metadata, top_n = 10)

pacientes <- unique(df_qiime_samples$NHC)

for (pac in pacientes) {
  df_pac <- df_qiime_samples %>%
    filter(NHC == pac) %>%
    mutate(
      Tiempo         = as.character(Tiempo),
      Tiempo_numeric = as.numeric(gsub("_.*", "", Tiempo))
    ) %>%
    arrange(Tiempo_numeric) %>%
    mutate(
      Tiempo = factor(Tiempo, levels = unique(Tiempo[order(Tiempo_numeric)]))
    ) %>%
    select(-Tiempo_numeric)
  
  p <- ggplot(df_pac, aes(x = Tiempo, y = AbundRel, fill = Taxon)) +
    geom_col(width = 0.9) +
    labs(
      title = paste("Composición en el paciente", pac),
      x     = "Tiempo",
      y     = "Abundancia relativa (%)"
    ) +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(
    paste0("0_figs/Paciente_NHC/barplots_relab/abund_rel_", pac, ".svg"),
    p,
    width  = 7,
    height = 5
  )
}

# ----------------------------------------------------------------------
# 7. EVOLUCIÓN TEMPORAL DIVERSIDAD ALFA (PACIENTES +=3 MUESTRAS)
# ----------------------------------------------------------------------
pacientes_con_3_muestras <- ps1.meta %>%
  group_by(NHC) %>%
  filter(n() >= 3) %>%
  ungroup() %>%
  pull(NHC) %>%
  unique()

df_alpha <- ps1.meta %>%
  filter(NHC %in% pacientes_con_3_muestras) %>%
  mutate(
    Tiempo      = as.character(Tiempo),
    Tiempo_num  = as.numeric(gsub("_.*", "", Tiempo)),
    Tiempo      = factor(Tiempo, levels = unique(Tiempo[order(Tiempo_num)]))
  )

# Shannon evolución completa
p <- ggplot(df_alpha, aes(x = Tiempo, y = Shannon, group = NHC)) +
  geom_point() +
  geom_line() +
  facet_wrap(~ NHC, scales = "free_x") +
  labs(
    title = "Evolución del índice de Shannon por paciente",
    x     = "Tiempo",
    y     = "Índice de Shannon"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("0_figs/Paciente_NHC/shannon_evolution_complete.svg", p, width = 10, height = 10)

# Shannon por paciente individual
for (pac in pacientes_con_3_muestras) {
  df_pac <- ps1.meta %>%
    filter(NHC == pac) %>%
    mutate(
      Tiempo         = as.character(Tiempo),
      Tiempo_num     = as.numeric(gsub("_.*", "", Tiempo)),
      Tiempo         = factor(Tiempo, levels = unique(Tiempo[order(Tiempo_num)]))
    ) %>%
    arrange(Tiempo)
  
  p <- ggplot(df_pac, aes(x = Tiempo, y = Shannon, group = 1)) +
    geom_point() +
    geom_line() +
    labs(
      title = paste("Evolución del índice de Shannon - Paciente", pac),
      x     = "Tiempo",
      y     = "Índice de Shannon"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(
    paste0("0_figs/Paciente_NHC/shannon_evolution/shannon_evol_", pac, ".svg"),
    p,
    width  = 7,
    height = 5
  )
}

# Chao1 evolución completa
p <- ggplot(df_alpha, aes(x = Tiempo, y = chao1, group = NHC)) +
  geom_point() +
  geom_line() +
  facet_wrap(~ NHC, scales = "free_x") +
  labs(
    title = "Evolución del índice de Chao1 por paciente",
    x     = "Tiempo",
    y     = "Índice de Chao1"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("0_figs/Paciente_NHC/chao1_evolution_complete.svg", p, width = 10, height = 10)

# Chao1 por paciente individual
for (pac in pacientes_con_3_muestras) {
  df_pac <- ps1.meta %>%
    filter(NHC == pac) %>%
    mutate(
      Tiempo         = as.character(Tiempo),
      Tiempo_num     = as.numeric(gsub("_.*", "", Tiempo)),
      Tiempo         = factor(Tiempo, levels = unique(Tiempo[order(Tiempo_num)]))
    ) %>%
    arrange(Tiempo)
  
  p <- ggplot(df_pac, aes(x = Tiempo, y = chao1, group = 1)) +
    geom_point() +
    geom_line() +
    labs(
      title = paste("Evolución del índice de Chao1 - Paciente", pac),
      x     = "Tiempo",
      y     = "Índice de Chao1"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(
    paste0("0_figs/Paciente_NHC/chao1_evolution/chao1_evol_", pac, ".svg"),
    p,
    width  = 7,
    height = 5
  )
}

# DBP evolución completa
p <- ggplot(df_alpha, aes(x = Tiempo, y = dbp, group = NHC)) +
  geom_point() +
  geom_line() +
  facet_wrap(~ NHC, scales = "free_x") +
  labs(
    title = "Evolución del índice de dominancia por paciente",
    x     = "Tiempo",
    y     = "dbp"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("0_figs/Paciente_NHC/dbp_evolution_complete.svg", p, width = 10, height = 10)

# DBP por paciente individual
for (pac in pacientes_con_3_muestras) {
  df_pac <- ps1.meta %>%
    filter(NHC == pac) %>%
    mutate(
      Tiempo         = as.character(Tiempo),
      Tiempo_num     = as.numeric(gsub("_.*", "", Tiempo)),
      Tiempo         = factor(Tiempo, levels = unique(Tiempo[order(Tiempo_num)]))
    ) %>%
    arrange(Tiempo)
  
  p <- ggplot(df_pac, aes(x = Tiempo, y = dbp, group = 1)) +
    geom_point() +
    geom_line() +
    labs(
      title = paste("Evolución del índice de dominancia - Paciente", pac),
      x     = "Tiempo",
      y     = "dbp"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(
    paste0("0_figs/Paciente_NHC/dbp_evolution/dbp_evol_", pac, ".svg"),
    p,
    width  = 7,
    height = 5
  )
}

# ----------------------------------------------------------------------
# 8. EVOLUCIÓN COMBINADA DE LAS 3 MÉTRICAS POR PACIENTE
# ----------------------------------------------------------------------
for (pac in pacientes_con_3_muestras) {
  df_pac <- ps1.meta %>%
    filter(NHC == pac) %>%
    mutate(
      Tiempo         = as.character(Tiempo),
      Tiempo_num     = as.numeric(gsub("_.*", "", Tiempo)),
      Tiempo         = factor(Tiempo, levels = unique(Tiempo[order(Tiempo_num)]))
    ) %>%
    arrange(Tiempo) %>%
    pivot_longer(
      cols      = c(Shannon, chao1, dbp),
      names_to  = "Metric",
      values_to = "Value"
    ) %>%
    mutate(Metric = factor(Metric, levels = c("Shannon", "chao1", "dbp")))
  
  p <- ggplot(df_pac, aes(x = Tiempo, y = Value, group = 1, color = Metric)) +
    geom_point() +
    geom_line() +
    facet_wrap(~ Metric, scales = "free_y") +
    labs(
      title = paste("Evolución de métricas de diversidad - Paciente", pac),
      x     = "Tiempo",
      y     = "Valor"
    ) +
    theme_minimal() +
    theme(
      axis.text.x        = element_text(angle = 45, hjust = 1),
      legend.position    = "none"
    )
  
  ggsave(
    paste0("0_figs/Paciente_NHC/evolucion_metricas/paciente_", pac, ".svg"),
    p,
    width  = 10,
    height = 5
  )
}
