#!/usr/bin/env Rscript
# ======================================================================
# ALPHA, BETA Y COMPOSICIÓN POR VARIABLE CLÍNICA (QIIME2)
# ----------------------------------------------------------------------
# Descripción:
#   1. Carga tabla QIIME2 y metadata (cruda y limpia).
#   2. Calcula diversidad alfa (Shannon, Chao1, dbp) por nivel de variable.
#   3. Calcula beta diversidad (Bray/Jaccard), NMDS, PERMANOVA y betadisper.
#   4. Genera composiciones de abundancia relativa por nivel de variable.
#
# Entradas:
#   - 1_input/aggregated_taxonomy_table_qiime.tsv
#   - 1_input/metadata.csv
#   - 0_tables/metadata_clean.tsv (si se ha generado antes)
#
# Salidas:
#   - 0_figs/Paciente_<VAR>/alfa_<var>.png
#   - 0_figs/Paciente_<VAR>/bray-curtis-nmds.png
#   - 0_figs/Paciente_<VAR>/jaccard-nmds.png
#   - 0_figs/Paciente_<VAR>/abund_rel_por_<var>.png
#
# Requisitos:
#   - R >= 4.0
#   - Paquetes: dplyr, ggplot2, RColorBrewer, tidyr, stringr, tidyverse,
#               broom, rstatix, vegan, microbiome, phyloseq, ggrepel,
#               data.table, patchwork, viridis
#
# Autor: Sandra Mingo-Ramirez
# Fecha: 2025
# ======================================================================

suppressPackageStartupMessages({
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
  library(data.table)
  library(patchwork)
  library(viridis)
})

set.seed(123)

PROJECT_DIR <- "/home/sandra/Projects/UCI_Variables/"
setwd(PROJECT_DIR)

dir.create("0_figs", showWarnings = FALSE)
dir.create("0_tables", showWarnings = FALSE)

# ----------------------------------------------------------------------
# 1. FUNCIONES AUXILIARES COMUNES
# ----------------------------------------------------------------------

# 1.1 Cargar tabla QIIME2 ya normalizada
load_qiime_table <- function(path) {
  tab <- fread(path, header = TRUE, data.table = FALSE, fill = TRUE)
  rownames(tab) <- tab[, 1]
  tab <- tab[, -1, drop = FALSE]
  
  rownames(tab) <- gsub("^d__", "", rownames(tab))
  rownames(tab) <- gsub("p__",  "", rownames(tab))
  rownames(tab) <- gsub("c__",  "", rownames(tab))
  rownames(tab) <- gsub("o__",  "", rownames(tab))
  rownames(tab) <- gsub("f__",  "", rownames(tab))
  rownames(tab) <- gsub("g__",  "", rownames(tab))
  colnames(tab) <- gsub("_L001", "", colnames(tab))
  
  tab
}

# 1.2 Cargar metadata original y limpiar 
make_metadata_clean <- function(meta_path_in, meta_path_out) {
  metadata <- read.csv(meta_path_in, header = TRUE, row.names = 1, check.names = FALSE)
  rownames(metadata) <- gsub("_L001", "", rownames(metadata))
  
  # Ejemplo de corrección específica (como en tu script de NHC)
  metadata <- metadata %>%
    mutate(
      Tiempo = ifelse(
        NHC == "314681" & Tiempo == 0,
        paste0(Tiempo, "_", ave(Tiempo, NHC, Tiempo, FUN = seq_along)),
        Tiempo
      )
    ) %>%
    filter(
      !(rownames(metadata) %in% c("179_S93", "V1-V2-S-86-2024-9_S23"))
    )
  
  metadata$SampleID <- rownames(metadata)
  metadata <- metadata %>% relocate(SampleID, .before = 1)
  
  write.table(metadata, file = meta_path_out, sep = "\t",
              quote = FALSE, row.names = FALSE)
  
  metadata
}

# 1.3 Cargar metadata limpia (si ya existe)
load_metadata_clean <- function(meta_path_clean) {
  read.table(meta_path_clean, header = TRUE, row.names = 1,
             check.names = FALSE, sep = "\t")
}

# 1.4 Sincronizar tabla y metadata
sync_table_metadata <- function(abund, meta) {
  samples <- intersect(colnames(abund), rownames(meta))
  if (length(samples) == 0) {
    stop("No hay muestras comunes entre tabla y metadata.")
  }
  abund <- abund[, samples, drop = FALSE]
  meta  <- meta[samples, , drop = FALSE]
  list(abund = abund, meta = meta)
}

# 1.5 Preparar composiciones (genérico)
prep_sample_composition <- function(abund_table, meta, var_name,
                                    top_n = 10, other_label = "Otros") {
  
  common_samples <- intersect(colnames(abund_table), rownames(meta))
  if (length(common_samples) == 0) {
    stop("No hay muestras coincidentes entre abund_table y meta.")
  }
  
  abund_table <- as.data.frame(abund_table[, common_samples, drop = FALSE])
  meta_use    <- meta[common_samples, , drop = FALSE]
  
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
      AbundRel = 100 * AbundRel,
      Variable = meta_use[Sample, var_name],
      Tiempo   = meta_use[Sample, "Tiempo"]
    )
  
  tax_means <- long_df %>%
    group_by(TaxonLast) %>%
    summarise(mean_abund = mean(AbundRel), .groups = "drop") %>%
    arrange(desc(mean_abund))
  
  top_taxa <- head(tax_means$TaxonLast, top_n)
  
  long_df <- long_df %>%
    mutate(
      Taxon = if_else(TaxonLast %in% top_taxa, TaxonLast, other_label)
    ) %>%
    group_by(Sample, Variable, Tiempo, Taxon) %>%
    summarise(AbundRel = sum(AbundRel), .groups = "drop")
  
  lev <- unique(long_df$Taxon)
  lev <- c(setdiff(lev, other_label), other_label)
  long_df$Taxon <- factor(long_df$Taxon, levels = lev)
  
  long_df
}

# 1.6 Análisis alfa + beta + composición para una variable categórica
run_alpha_beta_for_factor <- function(
    abund,
    meta,
    var_name,
    outdir
) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  
  if (!var_name %in% colnames(meta)) {
    warning("Variable no encontrada en metadata: ", var_name)
    return(NULL)
  }
  
  # Filtrar NA
  meta_use <- meta[!is.na(meta[[var_name]]), , drop = FALSE]
  sync     <- sync_table_metadata(abund, meta_use)
  abund_f  <- sync$abund
  meta_f   <- sync$meta
  
  # Phyloseq
  data_tab <- otu_table(abund_f, taxa_are_rows = TRUE)
  meta_tab <- sample_data(meta_f)
  pseq     <- merge_phyloseq(data_tab, meta_tab)
  
  # ------------------------------------------------------------------
  # ALPHA DIVERSITY
  # ------------------------------------------------------------------
  alpha_metrics <- microbiome::alpha(pseq, index = c("Shannon", "Chao1", "dbp"))
  
  ps_meta <- meta_f
  ps_meta$Shannon <- alpha_metrics$diversity_shannon
  ps_meta$chao1   <- alpha_metrics$chao1
  ps_meta$dbp     <- alpha_metrics$dominance_dbp
  ps_meta[[var_name]] <- factor(ps_meta[[var_name]])
  
  # Boxplots
  p1 <- ggplot(ps_meta, aes(x = .data[[var_name]], y = Shannon,
                            fill = .data[[var_name]])) +
    geom_boxplot() +
    theme_minimal() +
    labs(x = var_name, y = "Shannon") +
    theme(legend.position = "none")
  
  p2 <- ggplot(ps_meta, aes(x = .data[[var_name]], y = chao1,
                            fill = .data[[var_name]])) +
    geom_boxplot() +
    theme_minimal() +
    labs(x = var_name, y = "Chao1") +
    theme(legend.position = "none")
  
  p3 <- ggplot(ps_meta, aes(x = .data[[var_name]], y = dbp,
                            fill = .data[[var_name]])) +
    geom_boxplot() +
    theme_minimal() +
    labs(x = var_name, y = "dbp") +
    theme(legend.position = "none")
  
  (p1 | p2 | p3)
  ggsave(
    file.path(outdir, paste0("alfa_", var_name, ".png")),
    width = 10, height = 5, dpi = 300
  )
  
  # ------------------------------------------------------------------
  # BETA DIVERSITY (Bray & Jaccard)
  # ------------------------------------------------------------------
  dist_bc <- phyloseq::distance(pseq, method = "bray")
  nmds_bc <- vegan::metaMDS(dist_bc, k = 2, trymax = 100)
  
  permanova_res_bc <- adonis2(dist_bc ~ .data[[var_name]],
                              data = meta_f, permutations = 999)
  bd_bc <- betadisper(dist_bc, meta_f[[var_name]])
  permutest(bd_bc)
  
  nmds_bc_df <- as.data.frame(nmds_bc$points)
  nmds_bc_df$Sample <- rownames(nmds_bc_df)
  nmds_bc_df$Variable <- meta_f[nmds_bc_df$Sample, var_name]
  nmds_bc_df$Variable <- factor(nmds_bc_df$Variable)
  nmds_bc_df <- nmds_bc_df %>% drop_na(MDS1, MDS2)
  
  ellipse_groups_bc <- nmds_bc_df %>%
    group_by(Variable) %>%
    filter(
      n() > 3,
      sd(MDS1) > 0,
      sd(MDS2) > 0
    )
  
  p_bc <- ggplot(nmds_bc_df, aes(x = MDS1, y = MDS2, color = Variable)) +
    geom_point(size = 3, alpha = 0.8) +
    stat_ellipse(
      data = ellipse_groups_bc,
      aes(group = Variable, color = Variable),
      linetype = 1, size = 0.8, alpha = 0.4,
      show.legend = FALSE
    ) +
    theme_minimal(base_size = 16) +
    labs(
      x = "NMDS1", y = "NMDS2", color = var_name,
      title = paste0("PERMANOVA = ",
                     signif(permanova_res_bc$`Pr(>F)`[1], 3))
    )
  
  ggsave(
    file.path(outdir, "bray-curtis-nmds.png"),
    p_bc,
    width = 10, height = 10, dpi = 300
  )
  
  # Jaccard
  otu_pa <- data_tab
  otu_pa@.Data <- ifelse(otu_pa@.Data > 0, 1, 0)
  pseq_pa <- phyloseq(otu_pa, meta_tab)
  
  dist_jac <- phyloseq::distance(pseq_pa, method = "jaccard", binary = TRUE)
  nmds_jac <- vegan::metaMDS(dist_jac, k = 2, trymax = 100)
  
  permanova_res_jac <- adonis2(dist_jac ~ .data[[var_name]],
                               data = meta_f, permutations = 999)
  
  nmds_jac_df <- as.data.frame(nmds_jac$points)
  nmds_jac_df$Sample <- rownames(nmds_jac_df)
  nmds_jac_df$Variable <- meta_f[nmds_jac_df$Sample, var_name]
  nmds_jac_df$Variable <- factor(nmds_jac_df$Variable)
  nmds_jac_df <- nmds_jac_df %>% drop_na(MDS1, MDS2)
  
  ellipse_groups_jac <- nmds_jac_df %>%
    group_by(Variable) %>%
    filter(
      n() > 3,
      sd(MDS1) > 0,
      sd(MDS2) > 0
    )
  
  p_jac <- ggplot(nmds_jac_df, aes(x = MDS1, y = MDS2, color = Variable)) +
    geom_point(size = 3, alpha = 0.8) +
    stat_ellipse(
      data = ellipse_groups_jac,
      aes(group = Variable, color = Variable),
      linetype = 1, size = 0.8, alpha = 0.4,
      show.legend = FALSE
    ) +
    theme_minimal(base_size = 16) +
    labs(
      x = "NMDS1", y = "NMDS2", color = var_name,
      title = paste0("PERMANOVA = ",
                     signif(permanova_res_jac$`Pr(>F)`[1], 3))
    )
  
  ggsave(
    file.path(outdir, "jaccard-nmds.png"),
    p_jac,
    width = 10, height = 10, dpi = 300
  )
  
  # ------------------------------------------------------------------
  # COMPOSICIÓN
  # ------------------------------------------------------------------
  df_composition <- prep_sample_composition(
    abund_table = abund_f,
    meta        = meta_f,
    var_name    = var_name,
    top_n       = 10
  )
  
  df_summary <- df_composition %>%
    group_by(Variable, Taxon) %>%
    summarise(AbundRel = mean(AbundRel), .groups = "drop")
  
  p_comp <- ggplot(df_summary, aes(x = Variable, y = AbundRel, fill = Taxon)) +
    geom_col(width = 0.9) +
    labs(
      title = paste("Composición media por", var_name),
      x = var_name,
      y = "Abundancia relativa media (%)"
    ) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank())
  
  ggsave(
    file.path(outdir, paste0("abund_rel_por_", var_name, ".png")),
    p_comp,
    width = 7, height = 5, dpi = 300
  )
  
  invisible(list(
    alpha = ps_meta,
    permanova_bray = permanova_res_bc,
    permanova_jac  = permanova_res_jac
  ))
}

# ----------------------------------------------------------------------
# 2. CARGA DE DATOS BASE
# ----------------------------------------------------------------------
message("[1/2] Cargando tabla QIIME2 y metadata limpia...")

table_qiime <- load_qiime_table("1_input/aggregated_taxonomy_table_qiime.tsv")

# Si ya está metadata_clean.tsv, usa:
metadata_clean_path <- "0_tables/metadata_clean.tsv"
if (file.exists(metadata_clean_path)) {
  metadata <- load_metadata_clean(metadata_clean_path)
} else {
  metadata <- make_metadata_clean(
    meta_path_in  = "1_input/metadata.csv",
    meta_path_out = metadata_clean_path
  )
}

sync <- sync_table_metadata(table_qiime, metadata)
table_qiime <- sync$abund
metadata    <- sync$meta

message(sprintf("  - Muestras comunes: %d", ncol(table_qiime)))

# ----------------------------------------------------------------------
# 3. EJECUCIÓN POR VARIABLES FACTORIALES
# ----------------------------------------------------------------------
message("[2/2] Ejecutando análisis por variable...")

vars_factoriales <- c(
  "Sexo",
  "Inmunosuprimido",
  "VMNI.durante.el.ingreso",
  "Inmunosuprimido"
)

results <- list()

for (v in vars_factoriales) {
  message("  - Variable: ", v)
  outdir_v <- file.path("0_figs", paste0("Paciente_", v))
  res_v <- run_alpha_beta_for_factor(
    abund   = table_qiime,
    meta    = metadata,
    var_name = v,
    outdir  = outdir_v
  )
  results[[v]] <- res_v
}

cat("\n=== Análisis por variables (alpha, beta, composición) completado ===\n")
