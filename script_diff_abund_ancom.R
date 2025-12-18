#!/usr/bin/env Rscript
# ======================================================================
# ANCOM-BC2 EN DATOS DE MICROBIOMA LONGITUDINAL (QIIME2)
# ----------------------------------------------------------------------
# Descripción:
#   1. Carga tabla de abundancias QIIME2 y metadata limpia.
#   2. Construye un objeto phyloseq con recuentos y metadatos.
#   3. Ejecuta ANCOM-BC2 para variables clínicas (Sexo, Inmunosuprimido,
#      VMNI durante el ingreso).
#   4. Exporta tablas de resultados por variable y genera gráficos de
#      LFC para taxones significativos.
#
# Entradas:
#   - 1_input/aggregated_taxonomy_table_qiime.tsv
#   - 0_tables/metadata_clean.tsv
#
# Salidas:
#   - ANCOMBC2_Sexo.csv
#   - ANCOMBC2_Inmunosuprimido.csv
#   - ANCOMBC2_VMNI.csv
#   - ancombc-sexo.png
#   - ancombc-inmuno.png
#   - ancombc-vmni.png
#
# Requisitos:
#   - R >= 4.0
#   - Paquetes: ANCOMBC, phyloseq, data.table, dplyr, tidyverse, forcats
#
# Autor: Sandra Mingo-Ramirez
# Fecha: 2025
# ======================================================================

PROJECT_DIR <- "/home/sandra/Projects/UCI_Variables/"
setwd(PROJECT_DIR)

# ----------------------------------------------------------------------
# 1. LIBRERÍAS
# ----------------------------------------------------------------------
library(ANCOMBC)
library(phyloseq)
library(data.table)
library(dplyr)
library(tidyverse)
library(forcats)

set.seed(123)

# ----------------------------------------------------------------------
# 2. CARGA Y PREPARACIÓN DE DATOS
# ----------------------------------------------------------------------
metadata <- read.table(
  "0_tables/metadata_clean.tsv",
  header      = TRUE,
  row.names   = 1,
  check.names = FALSE,
  sep         = "\t"
)

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

samples     <- intersect(colnames(table_qiime), rownames(metadata))
abund_table <- table_qiime[, samples]
metadata    <- metadata[samples, ]

colnames(abund_table) <- rownames(metadata)

metadata$Sexo                      <- factor(metadata$Sexo)
metadata$Inmunosuprimido           <- factor(metadata$Inmunosuprimido)
metadata$VMNI.durante.el.ingreso   <- factor(metadata$VMNI.durante.el.ingreso)

otu_mat <- as.matrix(abund_table)
ps_otu  <- otu_table(otu_mat, taxa_are_rows = TRUE)
ps_meta <- sample_data(metadata)
ps      <- phyloseq(ps_otu, ps_meta)

# ----------------------------------------------------------------------
# 3. FUNCIÓN PARA EJECUTAR ANCOM-BC2 POR VARIABLE
# ----------------------------------------------------------------------
run_ancom_single <- function(ps_obj, variable) {
  
  var_levels <- levels(sample_data(ps_obj)[[variable]])
  ref_level  <- var_levels[1]
  
  res <- ancombc2(
    data         = ps_obj,
    fix_formula  = variable,
    rand_formula = NULL,
    p_adj_method = "fdr",
    prv_cut      = 0.10,
    lib_cut      = 1000,
    group        = variable,
    n_cl         = 1,
    verbose      = TRUE
  )
  
  contrast_levels <- var_levels[var_levels != ref_level]
  
  diff_cols <- paste0("diff_", variable, contrast_levels)
  p_cols    <- paste0("p_",   variable, contrast_levels)
  q_cols    <- paste0("q_",   variable, contrast_levels)
  lfc_cols  <- paste0("lfc_", variable, contrast_levels)
  
  out <- res$res %>%
    dplyr::select(
      taxon,
      dplyr::all_of(diff_cols),
      dplyr::all_of(p_cols),
      dplyr::all_of(q_cols),
      dplyr::all_of(lfc_cols)
    )
  
  out
}

# ----------------------------------------------------------------------
# 4. EJECUCIÓN POR VARIABLES
# ----------------------------------------------------------------------
res_sexo <- run_ancom_single(ps, "Sexo")
write.csv(res_sexo, "0_tables/ANCOMBC2_Sexo.csv", row.names = FALSE)

res_inmuno <- run_ancom_single(ps, "Inmunosuprimido")
write.csv(res_inmuno, "0_tables/ANCOMBC2_Inmunosuprimido.csv", row.names = FALSE)

res_vmni <- run_ancom_single(ps, "VMNI.durante.el.ingreso")
write.csv(res_vmni, "0_tables/ANCOMBC2_VMNI.csv", row.names = FALSE)

# ----------------------------------------------------------------------
# 5. FUNCIÓN DE GRÁFICO PARA RESULTADOS ANCOM-BC2
# ----------------------------------------------------------------------
plot_ancombc2 <- function(df, variable, q_threshold = 0.05, title = NULL) {
  
  var_levels <- unique(gsub(
    paste0(".*_", variable), "",
    grep(paste0("lfc_", variable), names(df), value = TRUE)
  ))
  contrast_levels <- var_levels[var_levels != var_levels[1]]
  
  lfc_cols <- paste0("lfc_", variable, contrast_levels)
  q_cols   <- paste0("q_",   variable, contrast_levels)
  
  df_long <- df %>%
    dplyr::select(taxon, dplyr::all_of(lfc_cols), dplyr::all_of(q_cols)) %>%
    tidyr::pivot_longer(
      cols         = -taxon,
      names_to     = c(".value", "contrast"),
      names_pattern = "(lfc|q)_(.+)"
    ) %>%
    dplyr::mutate(
      taxon    = sub(".*\\|", "", taxon),
      contrast = factor(contrast)
    )
  
  df_sig <- df_long %>% dplyr::filter(q < q_threshold)
  if (nrow(df_sig) == 0) {
    message("No hay taxones significativos para ", variable)
    return(NULL)
  }
  
  legend_labels <- c(
    "TRUE"  = paste("Mayor en", var_levels[2]),
    "FALSE" = paste("Mayor en", var_levels[1])
  )
  
  ggplot(df_sig, aes(
    x     = lfc,
    y     = forcats::fct_reorder(taxon, lfc),
    color = lfc > 0
  )) +
    geom_point(size = 3) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    scale_color_manual(
      values = c("TRUE" = "#1B9E77", "FALSE" = "#D95F02"),
      labels = legend_labels
    ) +
    labs(
      x     = "Log Fold Change (LFC)",
      y     = "Taxón",
      color = "Dirección",
      title = title
    ) +
    theme_minimal(base_size = 13)
}

# ----------------------------------------------------------------------
# 6. GRÁFICOS DE RESULTADOS
# ----------------------------------------------------------------------
plot_ancombc2(res_sexo, "Sexo", title = "ANCOM-BC2: Sexo")
ggsave("0_figs/ancombc-sexo.png", width = 12, height = 8)

plot_ancombc2(res_inmuno, "Inmunosuprimido", title = "ANCOM-BC2: Inmunosuprimido")
ggsave("0_figs/ancombc-inmuno.png", width = 12, height = 8)

plot_ancombc2(res_vmni, "VMNI.durante.el.ingreso", title = "ANCOM-BC2: VMNI")
ggsave("0_figs/ancombc-vmni.png", width = 12, height = 8)
