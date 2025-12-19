#!/usr/bin/env Rscript
# ======================================================================
# BETA DIVERSITY (NMDS) VS VARIABLES CLÍNICAS
# ----------------------------------------------------------------------
# Descripción:
#   1. Carga tabla QIIME2 y metadata.
#   2. Calcula distancia Bray-Curtis y NMDS global.
#   3. Para cada variable categórica:
#      - Genera un NMDS coloreado por niveles de la variable.
#      - Calcula PERMANOVA (adonis2), betadisper y ANOSIM (opcional).
#   4. Analiza SAPS II (continua) con NMDS, PERMANOVA y betadisper.
#   5. Exporta una tabla resumen con todos los estadísticos, incluyendo
#      p-valores ajustados por FDR para cada test.
#
# Entradas:
#   - 1_input/aggregated_taxonomy_table_qiime.tsv
#   - 1_input/metadata.csv
#
# Salidas:
#   - 0_figs/NMDS_por_variable/NMDS_*.svg
#   - 0_tables/adonis2_estadisticos_resumen.csv
#
# Requisitos:
#   - R >= 4.0
#   - Paquetes: data.table, dplyr, vegan, tidyverse
#
# Autor: Sandra Mingo-Ramirez
# Fecha: 2025
# ======================================================================

PROJECT_DIR <- "/home/sandra/Projects/UCI_Variables/"
setwd(PROJECT_DIR)

# ----------------------------------------------------------------------
# 1. LIBRERÍAS
# ----------------------------------------------------------------------
library(data.table)
library(dplyr)
library(vegan)
library(tidyverse)

dir.create("0_figs/NMDS_por_variable", recursive = TRUE, showWarnings = FALSE)
dir.create("0_tables", showWarnings = FALSE)

# ----------------------------------------------------------------------
# 2. TABLA QIIME2 Y METADATA
# ----------------------------------------------------------------------
table_qiime <- fread(
  "1_input/aggregated_taxonomy_table_qiime.tsv",
  header     = TRUE,
  data.table = FALSE,
  fill       = TRUE
)

rownames(table_qiime) <- table_qiime[, 1]
table_qiime <- table_qiime[, -1]

rownames(table_qiime) <- rownames(table_qiime) |>
  str_replace("^d__", "") |>
  str_replace("p__", "") |>
  str_replace("c__", "") |>
  str_replace("o__", "") |>
  str_replace("f__", "") |>
  str_replace("g__", "")

colnames(table_qiime) <- gsub("_L001", "", colnames(table_qiime))

metadata <- read.csv(
  "1_input/metadata.csv",
  header      = TRUE,
  row.names   = 1,
  check.names = FALSE
)
rownames(metadata) <- gsub("_L001", "", rownames(metadata))

samples <- intersect(colnames(table_qiime), rownames(metadata))
table_qiime <- table_qiime[, samples]
metadata    <- metadata[samples, ]

# ----------------------------------------------------------------------
# 3. DISTANCIA BRAY-CURTIS Y NMDS GLOBAL
# ----------------------------------------------------------------------
dist_mat <- vegdist(t(table_qiime), method = "bray")
nmds     <- metaMDS(dist_mat, trymax = 100)

run_nmds_analysis <- function(
    variable,
    var_name,
    metadata,
    dist_mat,
    nmds,
    file_name,
    use_anosim = TRUE
) {
  factor_var <- as.factor(metadata[[variable]])
  cols       <- as.numeric(factor_var)
  
  # NMDS
  svg(file_name, width = 8, height = 8)
  plot(
    nmds$points,
    col  = cols,
    pch  = 19,
    main = paste("NMDS -", var_name)
  )
  legend(
    "topright",
    legend = levels(factor_var),
    col    = seq_along(levels(factor_var)),
    pch    = 19
  )
  dev.off()
  
  # PERMANOVA
  adon <- adonis2(dist_mat ~ factor_var, data = metadata)
  
  # Betadisper
  beta      <- betadisper(dist_mat, factor_var, type = "centroid")
  beta_perm <- permutest(beta)
  
  # ANOSIM (opcional)
  if (use_anosim) {
    anos   <- anosim(dist_mat, factor_var)
    anos_R <- anos$statistic
    anos_p <- anos$signif
  } else {
    anos_R <- NA
    anos_p <- NA
  }
  
  tibble(
    Variable     = var_name,
    adonis2_F    = adon$F[1],
    adonis2_R2   = adon$R2[1],
    adonis2_p    = adon$`Pr(>F)`[1],
    Betadisper_F = beta_perm$tab$F[1],
    Betadisper_p = beta_perm$tab$`Pr(>F)`[1],
    ANOSIM_R     = anos_R,
    ANOSIM_p     = anos_p
  )
}

# ----------------------------------------------------------------------
# 4. VARIABLES CATEGÓRICAS
# ----------------------------------------------------------------------
results <- bind_rows(
  run_nmds_analysis("Batch",                 "Batch",                  metadata, dist_mat, nmds, "0_figs/NMDS_por_variable/NMDS_batch.svg"),
  run_nmds_analysis("Sexo",                  "Sexo",                   metadata, dist_mat, nmds, "0_figs/NMDS_por_variable/NMDS_sexo.svg"),
  run_nmds_analysis("Foco",                  "Foco",                   metadata, dist_mat, nmds, "0_figs/NMDS_por_variable/NMDS_foco.svg"),
  run_nmds_analysis("Inmunosuprimido",       "Inmunosuprimido",        metadata, dist_mat, nmds, "0_figs/NMDS_por_variable/NMDS_inmunosuprimido.svg"),
  run_nmds_analysis("NHC",                   "Paciente (NHC)",         metadata, dist_mat, nmds, "0_figs/NMDS_por_variable/NMDS_NHC.svg"),
  run_nmds_analysis("Mortalidad.UCI",        "Mortalidad UCI",         metadata, dist_mat, nmds, "0_figs/NMDS_por_variable/NMDS_mortalidad_UCI.svg"),
  run_nmds_analysis("Mortalidad.hospitalaria","Mortalidad hospitalaria",metadata, dist_mat, nmds, "0_figs/NMDS_por_variable/NMDS_mortalidad_hospitalaria.svg"),
  run_nmds_analysis("Traqueostomia",         "Traqueostomía",          metadata, dist_mat, nmds, "0_figs/NMDS_por_variable/NMDS_traqueostomia.svg"),
  run_nmds_analysis("VM.durante.el.ingreso", "Ventilación mecánica",   metadata, dist_mat, nmds, "0_figs/NMDS_por_variable/NMDS_VM.svg"),
  run_nmds_analysis("GNAF.durante.elingreso","GNAF",                   metadata, dist_mat, nmds, "0_figs/NMDS_por_variable/NMDS_GNAF.svg"),
  run_nmds_analysis("VMNI.durante.el.ingreso","VMNI",                  metadata, dist_mat, nmds, "0_figs/NMDS_por_variable/NMDS_VMNI.svg"),
  run_nmds_analysis("Riesgo.RZ",             "Riesgo RZ",              metadata, dist_mat, nmds, "0_figs/NMDS_por_variable/NMDS_RiesgoRZ.svg"),
  run_nmds_analysis("Bicho",                 "Microorganismo",         metadata, dist_mat, nmds, "0_figs/NMDS_por_variable/NMDS_Bicho.svg"),
  run_nmds_analysis("X.Ingresa.séptico.",    "Ingreso séptico",        metadata, dist_mat, nmds, "0_figs/NMDS_por_variable/NMDS_Ingreso_septico.svg")
)

# ----------------------------------------------------------------------
# 5. VARIABLE CONTINUA SAPS II
# ----------------------------------------------------------------------
metadata$SAPS.II.score <- as.numeric(metadata$SAPS.II.score)
metadata_clean <- metadata[!is.na(metadata$SAPS.II.score), ]

dist_matrix <- as.matrix(dist_mat)
dist_clean  <- as.dist(dist_matrix[rownames(metadata_clean), rownames(metadata_clean)])

nmds_saps <- metaMDS(dist_clean, trymax = 100)

metadata_clean$SAPS.II.group <- cut(
  metadata_clean$SAPS.II.score,
  breaks = 4,
  labels = c("Bajo", "Medio-Bajo", "Medio-Alto", "Alto")
)

svg("0_figs/NMDS_por_variable/NMDS_SAPSII.svg", width = 10, height = 10)
plot(
  nmds_saps$points,
  col  = as.numeric(metadata_clean$SAPS.II.group),
  pch  = 19,
  main = "NMDS - SAPS II"
)
legend(
  "topright",
  legend = levels(metadata_clean$SAPS.II.group),
  col    = seq_along(levels(metadata_clean$SAPS.II.group)),
  pch    = 19
)
dev.off()

adon_saps       <- adonis2(dist_clean ~ SAPS.II.score, data = metadata_clean)
beta_saps       <- betadisper(dist_clean, metadata_clean$SAPS.II.group)
beta_saps_perm  <- permutest(beta_saps)

results <- bind_rows(
  results,
  tibble(
    Variable     = "SAPS II",
    adonis2_F    = adon_saps$F[1],
    adonis2_R2   = adon_saps$R2[1],
    adonis2_p    = adon_saps$`Pr(>F)`[1],
    Betadisper_F = beta_saps_perm$tab$F[1],
    Betadisper_p = beta_saps_perm$tab$`Pr(>F)`[1],
    ANOSIM_R     = NA,
    ANOSIM_p     = NA
  )
)

# ----------------------------------------------------------------------
# 6. AJUSTE POR FDR Y EXPORTACIÓN
# ----------------------------------------------------------------------
results_fdr <- results %>%
  mutate(
    adonis2_p_FDR    = p.adjust(adonis2_p,    method = "fdr"),
    Betadisper_p_FDR = p.adjust(Betadisper_p, method = "fdr"),
    ANOSIM_p_FDR     = p.adjust(ANOSIM_p,     method = "fdr")
  )

write.csv(
  results_fdr,
  "0_tables/adonis2_estadisticos_resumen.csv",
  row.names = FALSE
)
