#!/usr/bin/env Rscript
# ======================================================================
# ANÁLISIS LONGITUDINAL DEL MICROBIOMA (QIIME2)
# ----------------------------------------------------------------------
# Descripción:
#   1. Carga tabla QIIME2 y metadata limpia.
#   2. Calcula dinámicas temporales de los taxones más abundantes.
#   3. Evalúa cambios temporales con Trendyspliner (por taxón y alfa).
#   4. Compara trayectorias por grupos clínicos (permuspliner).
#   5. Ajusta modelos GAM para efectos clínicos (Inmunosuprimido,
#      Traqueostomía, VMNI, GNAF) a nivel global y por taxón.
#   6. Entrena un mapa autoorganizado (SOM) para patrones de composición
#      y describe perfiles temporales por clúster.
#
# Entradas:
#   - 1_input/aggregated_taxonomy_table_qiime.tsv
#   - 0_tables/metadata_clean.tsv
#   - 0_tables/alpha_metrics_output.tsv
#
# Salidas principales:
#   - 0_figs/Paciente_longitudinal/*.png / *.svg
#   - 0_tables/trendy_result.tsv
#   - 0_tables/permu_sexo.tsv
#   - 0_tables/som_model.RData
#   - 0_tables/resultados_gam_*_por_taxon.tsv
#   - 0_tables/resultados_gam_*_global.tsv
#   - 0_tables/resultados_gam_*_temporal.tsv
#
# Requisitos:
#   - R >= 4.0
#   - Paquetes: paletteer, tidyverse, tidyr, splinectomeR, phyloseq,
#               microbiome, patchwork, mgcv, kohonen
#
# Autor: Sandra Mingo-Ramirez
# Fecha: 2025
# ======================================================================

rm(list = ls())
set.seed(123)

PROJECT_DIR <- "/home/sandra/Projects/UCI_Variables/"
setwd(PROJECT_DIR)

dir.create("0_figs/Paciente_longitudinal/", recursive = TRUE, showWarnings = FALSE)
dir.create("0_tables/", recursive = TRUE, showWarnings = FALSE)

# ----------------------------------------------------------------------
# 1. LIBRERÍAS
# ----------------------------------------------------------------------
library(paletteer)
library(tidyverse)
library(tidyr)
library(splinectomeR)
library(phyloseq)
library(microbiome)
library(patchwork)
library(mgcv)
library(kohonen)

# ----------------------------------------------------------------------
# 2. PALETA DE COLORES
# ----------------------------------------------------------------------
comp_colors <- c("#FFCC00", "#FF7000FF", "#920000FF", "#FF0000", "#006DDBFF",
                 "#6DB6FFFF", "#70A520", "#004949FF", "#B87EF2", "#490092FF")

# ----------------------------------------------------------------------
# 3. METADATA Y TABLA QIIME2
# ----------------------------------------------------------------------
metadata <- read.table(
  "0_tables/metadata_clean.tsv",
  header      = TRUE,
  row.names   = 1,
  check.names = FALSE,
  sep         = "\t"
)

metadata$Tiempo <- as.numeric(trimws(metadata$Tiempo))

metadata$Tiempo_grupo <- cut(
  metadata$Tiempo,
  breaks = c(-Inf, 0, 7, 21, Inf),
  labels = c("T0", "Temprano", "Intermedio", "Tardío")
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

samples <- intersect(colnames(table_qiime), rownames(metadata))
table_qiime <- table_qiime[, samples]
metadata    <- metadata[samples, ]

table_qiime <- table_qiime %>%
  rownames_to_column(var = "Taxon")

metadata <- metadata %>%
  rownames_to_column(var = "SampleID")

tax_long <- table_qiime %>%
  pivot_longer(
    cols      = -Taxon,
    names_to  = "SampleID",
    values_to = "abundance"
  )

tax_long <- tax_long %>%
  group_by(SampleID) %>%
  mutate(rel_abund = abundance / sum(abundance)) %>%
  ungroup()

df <- tax_long %>%
  left_join(metadata, by = "SampleID")

patients_keep <- df %>%
  group_by(NHC) %>%
  summarise(n = n_distinct(SampleID), .groups = "drop") %>%
  filter(n >= 3) %>%
  pull(NHC)

df <- df %>% filter(NHC %in% patients_keep)

top_taxa <- df %>%
  group_by(Taxon) %>%
  summarise(mean_ab = mean(rel_abund), .groups = "drop") %>%
  arrange(desc(mean_ab)) %>%
  slice(1:10) %>%
  pull(Taxon)

df_top <- df %>%
  filter(Taxon %in% top_taxa) %>%
  mutate(Taxon_simple = str_extract(Taxon, "[^|]+$")) %>%
  filter(rel_abund > 0)
  
# ----------------------------------------------------------------------
# 4. DINÁMICA DE LOS 10 TAXONES MÁS ABUNDANTES
# ----------------------------------------------------------------------
df_top_plot <- df_top %>%
  group_by(Tiempo, Taxon_simple) %>%
  summarise(
    Mean_rel_abund = mean(rel_abund, na.rm = TRUE),
    SD_rel_abund   = sd(rel_abund, na.rm = TRUE),
    .groups        = "drop"
  )

top10_abund_plot <- ggplot(
  df_top_plot,
  aes(x = Tiempo, y = Mean_rel_abund * 100,
      color = Taxon_simple, group = Taxon_simple)
) +
  geom_smooth(method = "loess", se = FALSE, size = 1) +
  geom_smooth(
    aes(fill = Taxon_simple),
    method      = "loess",
    se          = TRUE,
    alpha       = 0.075,
    size        = 0.55,
    show.legend = FALSE
  ) +
  labs(
    x     = "Tiempo (días)",
    y     = "Abundancia relativa media (%)",
    color = "Taxón",
    title = "Dinámica temporal de los 10 géneros bacterianos más abundantes"
  ) +
  scale_color_manual(values = comp_colors) +
  scale_fill_manual(values = comp_colors) +
  scale_x_continuous(
    breaks = seq(0, max(df_top_plot$Tiempo, na.rm = TRUE), by = 2)
  ) +
  coord_cartesian(ylim = c(0, 60)) +
  theme_minimal() +
  theme(legend.text = element_text(face = "italic"))

print(top10_abund_plot)
svg("0_figs/Paciente_longitudinal/top10_tax_abundance.svg", width = 10, height = 10)
print(top10_abund_plot)
dev.off()

# ----------------------------------------------------------------------
# 5. TRENDYSPLINER POR TAXÓN
# ----------------------------------------------------------------------
splinectome_df <- df_top %>%
  mutate(Taxon = Taxon_simple)

trendyspliner_top_species <- function(long_data) {
  results_df <- data.frame(
    Taxon = character(),
    pval  = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (s in unique(long_data$Taxon)) {
    s__df <- dplyr::filter(long_data, Taxon == s)
    cat(s)
    s.result <- trendyspliner(
      as.data.frame(s__df),
      xvar  = "Tiempo",
      yvar  = "rel_abund",
      cases = "SampleID",
      quiet = TRUE,
      perms = 999
    )
    cat(paste0(", p = ", s.result$pval, "\n"))
    results_df <- rbind(
      results_df,
      data.frame(
        Taxon = s,
        pval  = s.result$pval,
        stringsAsFactors = FALSE
      )
    )
  }
  
  results_df$pval_adj <- p.adjust(results_df$pval, method = "fdr")
  results_df$signif <- sapply(results_df$pval_adj, function(p) {
    if (p < 0.001) {
      "***"
    } else if (p < 0.01) {
      "**"
    } else if (p < 0.05) {
      "*"
    } else {
      ""
    }
  })
  
  results_df
}

trendy_result <- trendyspliner_top_species(splinectome_df)

write.table(
  trendy_result,
  "0_tables/trendy_result.tsv",
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE
)

# ----------------------------------------------------------------------
# 6. PLOT DE TAXONES SIGNIFICATIVOS (TRENDYSPLINER)
# ----------------------------------------------------------------------
sign_taxa <- trendy_result %>%
  filter(pval_adj < 0.05) %>%
  pull(Taxon)

df_sign <- df_top %>%
  filter(Taxon_simple %in% sign_taxa) %>%
  mutate(Tiempo_grupo = factor(
    Tiempo_grupo,
    levels = c("T0", "Temprano", "Intermedio", "Tardío")
  ))

df_sign_plot <- df_sign %>%
  group_by(Tiempo_grupo, Taxon_simple) %>%
  summarise(
    Mean_rel_abund = mean(rel_abund, na.rm = TRUE),
    SD_rel_abund   = sd(rel_abund, na.rm = TRUE),
    .groups        = "drop"
  )

sign_taxa_plot <- ggplot(
  df_sign_plot,
  aes(x = Tiempo_grupo, y = Mean_rel_abund * 100,
      color = Taxon_simple, group = Taxon_simple)
) +
  geom_smooth(method = "loess", se = FALSE, size = 1) +
  geom_smooth(
    aes(fill = Taxon_simple),
    method      = "loess",
    se          = TRUE,
    alpha       = 0.075,
    size        = 0.55,
    show.legend = FALSE
  ) +
  labs(
    x     = "Grupo de tiempo",
    y     = "Abundancia relativa media (%)",
    color = "Taxón",
    title = "Dinámica temporal de taxones con cambios significativos (Trendyspliner)"
  ) +
  facet_wrap(~ Taxon_simple, scales = "free_y", ncol = 4) +
  coord_cartesian(ylim = c(0, 60)) +
  scale_color_manual(values = comp_colors) +
  scale_fill_manual(values = comp_colors) +
  theme_minimal() +
  theme(
    legend.text = element_text(face = "italic"),
    strip.text  = element_text(face = "italic", size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(sign_taxa_plot)
png("0_figs/Paciente_longitudinal/sign_tax_abundance.png",
    width = 1000, height = 700, res = 120)
print(sign_taxa_plot)
dev.off()

# ----------------------------------------------------------------------
# 7. DINÁMICA DE DIVERSIDAD ALFA
# ----------------------------------------------------------------------
meta_alpha <- read.table(
  "0_tables/alpha_metrics_output.tsv",
  header = TRUE,
  sep    = "\t"
)

patients_keep <- meta_alpha %>%
  group_by(NHC) %>%
  summarise(n_samples = n_distinct(SampleID), .groups = "drop") %>%
  filter(n_samples >= 3) %>%
  pull(NHC)

meta_alpha <- meta_alpha %>%
  filter(NHC %in% patients_keep)

meta_alpha$Tiempo <- as.numeric(trimws(meta_alpha$Tiempo))

adiv_plot <- meta_alpha %>%
  group_by(Tiempo) %>%
  summarise(
    Mean_Shannon = mean(Shannon, na.rm = TRUE),
    Mean_Chao1   = mean(chao1,   na.rm = TRUE),
    Mean_dbp     = mean(dbp,     na.rm = TRUE),
    .groups      = "drop"
  ) %>%
  pivot_longer(
    cols      = starts_with("Mean_"),
    names_to  = "Index",
    values_to = "Mean"
  ) %>%
  mutate(
    Index = str_replace(Index, "Mean_", ""),
    Index = factor(Index, levels = c("Shannon", "Chao1", "dbp"))
  )

alpha_div_plot <- ggplot(adiv_plot, aes(x = Tiempo, y = Mean, group = Index)) +
  geom_smooth(
    method = "loess",
    se     = TRUE,
    color  = "#1f77b4",
    fill   = "#1f77b4",
    alpha  = 0.2,
    size   = 1
  ) +
  facet_wrap(~ Index, scales = "free_y") +
  labs(
    x     = "Tiempo (días)",
    y     = "Media de diversidad alfa",
    title = "Dinámica temporal de la diversidad alfa (LOESS)"
  ) +
  scale_x_continuous(
    breaks = seq(0, max(adiv_plot$Tiempo, na.rm = TRUE), by = 2)
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 11, face = "bold"),
    plot.title = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 12)
  )

print(alpha_div_plot)
svg("0_figs/Paciente_longitudinal/alpha_temporal.svg",
    width = 12, height = 7)
print(alpha_div_plot)
dev.off()

# Trendyspliner Shannon
alpha_trendresult_all <- trendyspliner(
  data   = as.data.frame(meta_alpha),
  xvar   = "Tiempo",
  yvar   = "Shannon",
  cases  = "NHC",
  perms  = 999,
  quiet  = TRUE
)
cat("Shannon p-value =", alpha_trendresult_all$pval, "\n")

shannon_trendy_plot <- trendyspliner.plot.perms(
  alpha_trendresult_all,
  xlabel = "Tiempo (días)",
  ylabel = "Diversidad de Shannon"
)

# Trendyspliner Chao1
alpha_trendresult_all <- trendyspliner(
  data   = as.data.frame(meta_alpha),
  xvar   = "Tiempo",
  yvar   = "chao1",
  cases  = "NHC",
  perms  = 999,
  quiet  = TRUE
)
cat("Chao1 p-value =", alpha_trendresult_all$pval, "\n")

chao_trendy_plot <- trendyspliner.plot.perms(
  alpha_trendresult_all,
  xlabel = "Tiempo (días)",
  ylabel = "Diversidad de Chao1"
)

# Trendyspliner dbp
alpha_trendresult_all <- trendyspliner(
  data   = as.data.frame(meta_alpha),
  xvar   = "Tiempo",
  yvar   = "dbp",
  cases  = "NHC",
  perms  = 999,
  quiet  = TRUE
)
cat("dbp p-value =", alpha_trendresult_all$pval, "\n")

dbp_trendy_plot <- trendyspliner.plot.perms(
  alpha_trendresult_all,
  xlabel = "Tiempo (días)",
  ylabel = "Dominancia dbp"
)

alpha_trendy_all <- (shannon_trendy_plot / chao_trendy_plot / dbp_trendy_plot) +
  plot_layout(ncol = 1)

svg("0_figs/Paciente_longitudinal/alpha_trendy_all.svg",
    width = 8, height = 10)
print(alpha_trendy_all)
dev.off()

# ----------------------------------------------------------------------
# 8. PERMUSPLINER POR VARIABLE (SEXO)
# ----------------------------------------------------------------------
permuspliner_by_variable <- function(long_data, variable, min_samples_per_group = 8) {
  
  long_data <- long_data %>%
    filter(
      !is.na(.data[[variable]]),
      !is.na(Tiempo),
      !is.na(rel_abund),
      !is.na(SampleID),
      rel_abund > 0,
      Tiempo > 0,
      is.finite(Tiempo)
    )
  
  groups <- unique(long_data[[variable]])
  g1 <- groups[1]
  g2 <- groups[2]
  message("Comparando: ", g1, " vs ", g2)
  
  taxon_counts <- long_data %>%
    count(Taxon, .data[[variable]], SampleID) %>%
    count(Taxon, .data[[variable]], name = "n_samples") %>%
    pivot_wider(
      names_from  = all_of(variable),
      values_from = n_samples,
      values_fill = 0
    ) %>%
    filter(if_all(all_of(c(g1, g2)), ~ . >= min_samples_per_group))
  
  taxones_validos <- taxon_counts$Taxon
  message("Taxones válidos (≥", min_samples_per_group,
          " muestras/grupo): ", length(taxones_validos))
  
  results_df <- data.frame(
    Taxon      = character(),
    Comparison = character(),
    pval       = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (s in taxones_validos) {
    s__df <- long_data %>%
      filter(
        Taxon == s,
        complete.cases(Tiempo, rel_abund, SampleID, .data[[variable]])
      )
    
    n_cases_g1 <- length(unique(s__df$SampleID[s__df[[variable]] == g1]))
    n_cases_g2 <- length(unique(s__df$SampleID[s__df[[variable]] == g2]))
    
    if (n_cases_g1 >= 3 & n_cases_g2 >= 3) {
      tryCatch({
        s.result <- permuspliner(
          as.data.frame(s__df),
          xvar    = "Tiempo",
          yvar    = "rel_abund",
          cases   = "SampleID",
          category = variable,
          groups  = c(g1, g2),
          quiet   = TRUE,
          perms   = 999
        )
        results_df <- rbind(
          results_df,
          data.frame(
            Taxon      = s,
            Comparison = paste0(variable, ":", g1, "_vs_", g2),
            pval       = s.result$pval,
            stringsAsFactors = FALSE
          )
        )
      }, error = function(e) {
        message("SKIPPED ", s, ": ", e$message)
      })
    }
  }
  
  results_df$pval_adj <- p.adjust(results_df$pval, "fdr")
  results_df$signif <- cut(
    results_df$pval_adj,
    breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
    labels = c("***", "**", "*", "")
  )
  results_df
}

splinectome_df$Sexo[splinectome_df$Sexo == ""] <- NA

set.seed(123)

permu_sexo <- permuspliner_by_variable(splinectome_df, "Sexo")
write.table(
  permu_sexo,
  "0_tables/permu_sexo.tsv",
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE
)

sig_taxa <- permu_sexo %>%
  filter(pval_adj < 0.05) %>%
  pull(Taxon)

for (taxon_i in sig_taxa) {
  df_i <- splinectome_df %>%
    filter(Taxon == taxon_i, !is.na(Sexo))
  
  res_i <- permuspliner(
    data       = as.data.frame(df_i),
    xvar       = "Tiempo",
    yvar       = "rel_abund",
    cases      = "SampleID",
    category   = "Sexo",
    groups     = c("Hombre", "Mujer"),
    perms      = 999,
    retain_perm = TRUE,
    quiet      = TRUE
  )
  
  p <- permuspliner.plot.permsplines(
    data = res_i,
    xvar = "Tiempo",
    yvar = "rel_abund"
  )
  
  print(p)
  fname <- paste0("0_figs/Paciente_longitudinal/permusplines_sexo_", taxon_i, ".svg")
  ggsave(fname, plot = p, width = 6, height = 5)
}

# Permuspliner para alfa-diversidad por sexo
res_shannon <- permuspliner(
  data       = meta_alpha,
  xvar       = "Tiempo",
  yvar       = "Shannon",
  cases      = "SampleID",
  category   = "Sexo",
  groups     = c("Hombre", "Mujer"),
  perms      = 999,
  retain_perm = TRUE,
  quiet      = TRUE
)

p_shannon <- permuspliner.plot.permsplines(
  data = res_shannon,
  xvar = "Tiempo",
  yvar = "Shannon"
)

res_chao <- permuspliner(
  data       = meta_alpha,
  xvar       = "Tiempo",
  yvar       = "chao1",
  cases      = "SampleID",
  category   = "Sexo",
  groups     = c("Hombre", "Mujer"),
  perms      = 999,
  retain_perm = TRUE,
  quiet      = TRUE
)

p_chao <- permuspliner.plot.permsplines(
  data = res_chao,
  xvar = "Tiempo",
  yvar = "chao1"
)

res_dbp <- permuspliner(
  data       = meta_alpha,
  xvar       = "Tiempo",
  yvar       = "dbp",
  cases      = "SampleID",
  category   = "Sexo",
  groups     = c("Hombre", "Mujer"),
  perms      = 999,
  retain_perm = TRUE,
  quiet      = TRUE
)

p_dbp <- permuspliner.plot.permsplines(
  data = res_dbp,
  xvar = "Tiempo",
  yvar = "dbp"
)

alpha_permu_sexo <- (p_shannon / p_chao / p_dbp) +
  plot_layout(ncol = 1)

svg("0_figs/Paciente_longitudinal/alpha_permu_sexo.svg",
    width = 8, height = 10)
print(alpha_permu_sexo)
dev.off()

# ----------------------------------------------------------------------
# 9. GAM: EFECTO INMUNOSUPRIMIDO (GLOBAL Y POR TAXÓN)
# ----------------------------------------------------------------------
splinectome_df$Inmunosuprimido[splinectome_df$Inmunosuprimido == ""] <- NA
splinectome_df <- splinectome_df %>%
  filter(!is.na(rel_abund), !is.na(Tiempo))

gam_model_inmuno <- gam(
  rel_abund ~ Inmunosuprimido +
    s(Tiempo, k = 4, bs = "cr") +
    s(NHC, bs = "re"),
  data   = splinectome_df,
  method = "REML"
)

summary(gam_model_inmuno)

coefs_inmuno <- summary(gam_model_inmuno)$p.table %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Termino") %>%
  filter(Termino != "s(Tiempo)") %>%
  mutate(
    Termino = case_when(
      Termino == "(Intercept)"         ~ "Intercepto",
      Termino == "InmunosuprimidoSí"   ~ "Inmunosuprimido (Sí vs No)",
      TRUE                             ~ Termino
    )
  ) %>%
  select(
    Termino,
    Estimate    = Estimate,
    Std.Error   = `Std. Error`,
    t.value     = `t value`,
    p.value     = `Pr(>|t|)`
  )

smooths_inmuno <- summary(gam_model_inmuno)$s.table %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Smooth") %>%
  mutate(p.value = `p-value`) %>%
  select(Smooth, edf, Ref.df = `Ref.df`, F = F, p.value)

cat("=== RESULTADOS FINALES INMUNOSUPRIMIDOS ===\n")
cat("n total:", nrow(splinectome_df), "\n")
cat("Inmunosuprimidos:",
    sum(splinectome_df$Inmunosuprimido == "Sí", na.rm = TRUE), "\n")
cat(sprintf(
  "Efecto inmunosupresión: β = %.3f (SE = %.3f), p = %.4f\n",
  coefs_inmuno$Estimate[2], coefs_inmuno$Std.Error[2], coefs_inmuno$p.value[2]
))
cat(sprintf(
  "Efecto tiempo: edf = %.2f, F = %.1f, p = %.4f\n",
  smooths_inmuno$edf[1], smooths_inmuno$F[1], smooths_inmuno$p.value[1]
))
cat("R² ajustado:", round(summary(gam_model_inmuno)$r.sq, 4), "\n")

write_tsv(coefs_inmuno, "0_tables/resultados_gam_inmunodeprimido_global.tsv")
write_tsv(smooths_inmuno, "0_tables/resultados_gam_inmunodeprimido_temporal.tsv")

df_top$Inmunosuprimido[df_top$Inmunosuprimido == ""] <- NA

resultados_taxones <- purrr::map_dfr(
  unique(df_top$Taxon_simple),
  function(taxon) {
    datos_taxon <- df_top %>%
      filter(Taxon_simple == taxon, !is.na(Inmunosuprimido)) %>%
      mutate(Inmunosuprimido = factor(Inmunosuprimido, levels = c("No", "Sí")))
    
    n_inmuno <- sum(datos_taxon$Inmunosuprimido == "Sí", na.rm = TRUE)
    n_total  <- nrow(datos_taxon)
    
    if (n_total > 20 & n_inmuno > 2) {
      modelo <- gam(
        rel_abund ~ Inmunosuprimido +
          s(Tiempo, k = 3, bs = "cr") +
          s(NHC, bs = "re"),
        data   = datos_taxon,
        method = "REML",
        family = gaussian()
      )
      
      coef_table <- summary(modelo)$p.table %>%
        as.data.frame() %>%
        tibble::rownames_to_column("Term") %>%
        filter(Term == "InmunosuprimidoSí")
      
      coef_inmuno <- tibble(
        Taxon    = taxon,
        n_total  = n_total,
        n_inmuno = n_inmuno,
        beta     = coef_table$Estimate,
        se       = coef_table$`Std. Error`,
        t_value  = coef_table$`t value`,
        p_value  = coef_table$`Pr(>|t|)`
      )
      
      assign(paste0("gam_", taxon), modelo, envir = .GlobalEnv)
      return(coef_inmuno)
    }
    NULL
  }
)

tabla_final <- resultados_taxones %>%
  arrange(p_value) %>%
  mutate(
    p_value_adj = p.adjust(p_value, method = "BH"),
    sig = case_when(
      p_value_adj < 0.05 ~ "***",
      p_value     < 0.05 ~ "**",
      p_value     < 0.1  ~ "*",
      TRUE              ~ ""
    ),
    beta        = round(beta, 4),
    se          = round(se, 4),
    p_value     = round(p_value, 4),
    p_value_adj = round(p_value_adj, 4)
  ) %>%
  select(Taxon, n_total, n_inmuno, beta, se, p_value, p_value_adj, sig)

write_tsv(tabla_final, "0_tables/resultados_gam_inmunosupresion_por_taxon.tsv")

df_top %>%
  filter(Taxon_simple %in% c("Enterococcus", "Bacteroides")) %>%
  ggplot(aes(Inmunosuprimido, rel_abund, fill = Inmunosuprimido)) +
  geom_boxplot(alpha = 0.7) +
  scale_y_log10(labels = scales::percent) +
  facet_wrap(~ Taxon_simple, scales = "free_y") +
  labs(title = "Efectos significativos de inmunosupresión") +
  theme_bw()

ggsave("0_figs/Paciente_longitudinal/boxplot_sign_taxa_inmuno.svg",
       width = 10, height = 8)

# ----------------------------------------------------------------------
# 10. GAM: TRAQUEOSTOMÍA, VMNI, GNAF (RESUMEN)
# ----------------------------------------------------------------------
# TRAQUEOSTOMÍA 
splinectome_df$Traqueostomia[splinectome_df$Traqueostomia == ""] <- NA
splinectome_df$Traqueostomia <- factor(splinectome_df$Traqueostomia,
                                       levels = c("No", "Sí"))

datos_traqueo <- splinectome_df %>%
  filter(!is.na(rel_abund), !is.na(Tiempo))

gam_model_traqueo <- gam(
  rel_abund ~ Traqueostomia +
    s(Tiempo, k = 4, bs = "cr") +
    s(NHC, bs = "re"),
  data   = datos_traqueo,
  method = "REML"
)

summary(gam_model_traqueo)

coefs_traqueo <- summary(gam_model_traqueo)$p.table %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Termino") %>%
  filter(Termino != "s(Tiempo)") %>%
  mutate(
    Termino = case_when(
      Termino == "(Intercept)"          ~ "Intercepto",
      Termino == "TraqueostomiaSí"      ~ "Traqueostomía (Sí vs No)",
      TRUE                              ~ Termino
    )
  ) %>%
  select(
    Termino,
    Estimate    = Estimate,
    Std.Error   = `Std. Error`,
    t.value     = `t value`,
    p.value     = `Pr(>|t|)`
  )

smooths_traqueo <- summary(gam_model_traqueo)$s.table %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Smooth") %>%
  mutate(p.value = `p-value`) %>%
  select(Smooth, edf, Ref.df = `Ref.df`, F = F, p.value)

write_tsv(coefs_traqueo, "0_tables/resultados_gam_traqueostomia_global.tsv")
write_tsv(smooths_traqueo, "0_tables/resultados_gam_traqueostomia_temporal.tsv")

df_top$Traqueostomia[df_top$Traqueostomia == ""] <- NA
df_top$Traqueostomia <- factor(df_top$Traqueostomia,
                               levels = c("No", "Sí"))

resultados_traqueo_taxones <- purrr::map_dfr(
  unique(df_top$Taxon_simple),
  function(taxon) {
    datos_taxon <- df_top %>%
      filter(Taxon_simple == taxon, !is.na(Traqueostomia)) %>%
      mutate(Traqueostomia = factor(Traqueostomia, levels = c("No", "Sí")))
    
    n_traq <- sum(datos_taxon$Traqueostomia == "Sí", na.rm = TRUE)
    n_total  <- nrow(datos_taxon)
    
    if (n_total > 20 & n_traq > 2) {
      modelo <- gam(
        rel_abund ~ Traqueostomia +
          s(Tiempo, k = 3, bs = "cr") +
          s(NHC, bs = "re"),
        data   = datos_taxon,
        method = "REML",
        family = gaussian()
      )
      
      coef_table <- summary(modelo)$p.table %>%
        as.data.frame() %>%
        tibble::rownames_to_column("Term") %>%
        filter(Term == "TraqueostomiaSí")
      
      coef_traq <- tibble(
        Taxon    = taxon,
        n_total  = n_total,
        n_traq   = n_traq,
        beta     = coef_table$Estimate,
        se       = coef_table$`Std. Error`,
        t_value  = coef_table$`t value`,
        p_value  = coef_table$`Pr(>|t|)`
      )
      
      assign(paste0("gam_traq_", taxon), modelo, envir = .GlobalEnv)
      return(coef_traq)
    }
    NULL
  }
)

tabla_traqueo <- resultados_traqueo_taxones %>%
  arrange(p_value) %>%
  mutate(
    p_value_adj = p.adjust(p_value, method = "BH"),
    sig = case_when(
      p_value_adj < 0.05 ~ "***",
      p_value     < 0.05 ~ "**",
      p_value     < 0.1  ~ "*",
      TRUE              ~ ""
    ),
    beta        = round(beta, 4),
    se          = round(se, 4),
    p_value     = round(p_value, 4),
    p_value_adj = round(p_value_adj, 4)
  ) %>%
  select(Taxon, n_total, n_traq, beta, se, p_value, p_value_adj, sig)

write_tsv(tabla_traqueo,
          "0_tables/resultados_gam_traqueostomia_por_taxon.tsv")

# VENTILACIÓN MECÁNICA NO INVASIVA (VMNI)
splinectome_df$VMNI.durante.el.ingreso[splinectome_df$VMNI.durante.el.ingreso == ""] <- NA
splinectome_df$VMNI.durante.el.ingreso <- factor(splinectome_df$VMNI.durante.el.ingreso,
                                       levels = c("No", "Sí"))

datos_vmni <- splinectome_df %>%
  filter(!is.na(rel_abund), !is.na(Tiempo))

gam_model_vmni <- gam(
  rel_abund ~ VMNI.durante.el.ingreso +
    s(Tiempo, k = 4, bs = "cr") +
    s(NHC, bs = "re"),
  data   = datos_vmni,
  method = "REML"
)

summary(gam_model_vmni)

coefs_vmni <- summary(gam_model_vmni)$p.table %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Termino") %>%
  filter(Termino != "s(Tiempo)") %>%
  mutate(
    Termino = case_when(
      Termino == "(Intercept)"          ~ "Intercepto",
      Termino == "VMNI.durante.el.ingresoSí"      ~ "VMNI.durante.el.ingreso (Sí vs No)",
      TRUE                              ~ Termino
    )
  ) %>%
  select(
    Termino,
    Estimate    = Estimate,
    Std.Error   = `Std. Error`,
    t.value     = `t value`,
    p.value     = `Pr(>|t|)`
  )

smooths_vmni <- summary(gam_model_vmni)$s.table %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Smooth") %>%
  mutate(p.value = `p-value`) %>%
  select(Smooth, edf, Ref.df = `Ref.df`, F = F, p.value)

write_tsv(coefs_vmni, "0_tables/resultados_gam_vmni_global.tsv")
write_tsv(smooths_vmni, "0_tables/resultados_gam_vmni_temporal.tsv")

df_top$VMNI.durante.el.ingreso[df_top$VMNI.durante.el.ingreso == ""] <- NA
df_top$VMNI.durante.el.ingreso <- factor(df_top$VMNI.durante.el.ingreso,
                               levels = c("No", "Sí"))

resultados_vmni_taxones <- purrr::map_dfr(
  unique(df_top$Taxon_simple),
  function(taxon) {
    datos_taxon <- df_top %>%
      filter(Taxon_simple == taxon, !is.na(VMNI.durante.el.ingreso)) %>%
      mutate(VMNI.durante.el.ingreso = factor(VMNI.durante.el.ingreso, levels = c("No", "Sí")))
    
    n_vmni <- sum(datos_taxon$VMNI.durante.el.ingreso == "Sí", na.rm = TRUE)
    n_total  <- nrow(datos_taxon)
    
    if (n_total > 20 & n_vmni > 2) {
      modelo <- gam(
        rel_abund ~ VMNI.durante.el.ingreso +
          s(Tiempo, k = 3, bs = "cr") +
          s(NHC, bs = "re"),
        data   = datos_taxon,
        method = "REML",
        family = gaussian()
      )
      
      coef_table <- summary(modelo)$p.table %>%
        as.data.frame() %>%
        tibble::rownames_to_column("Term") %>%
        filter(Term == "VMNI.durante.el.ingresoSí")
      
      coef_vmni <- tibble(
        Taxon    = taxon,
        n_total  = n_total,
        n_vmni   = n_vmni,
        beta     = coef_table$Estimate,
        se       = coef_table$`Std. Error`,
        t_value  = coef_table$`t value`,
        p_value  = coef_table$`Pr(>|t|)`
      )
      
      assign(paste0("gam_vmni_", taxon), modelo, envir = .GlobalEnv)
      return(coef_vmni)
    }
    NULL
  }
)

tabla_vmni <- resultados_vmni_taxones %>%
  arrange(p_value) %>%
  mutate(
    p_value_adj = p.adjust(p_value, method = "BH"),
    sig = case_when(
      p_value_adj < 0.05 ~ "***",
      p_value     < 0.05 ~ "**",
      p_value     < 0.1  ~ "*",
      TRUE              ~ ""
    ),
    beta        = round(beta, 4),
    se          = round(se, 4),
    p_value     = round(p_value, 4),
    p_value_adj = round(p_value_adj, 4)
  ) %>%
  select(Taxon, n_total, n_vmni, beta, se, p_value, p_value_adj, sig)

write_tsv(tabla_vmni,
          "0_tables/resultados_gam_vmni_por_taxon.tsv")

df_top %>%
  filter(Taxon_simple %in% c("Clostridium_innocuum_group", "Prevotella_9")) %>%
  ggplot(aes(VMNI.durante.el.ingreso, rel_abund, fill = VMNI.durante.el.ingreso)) +
  geom_boxplot(alpha = 0.7) +
  scale_y_log10(labels = scales::percent) +
  facet_wrap(~ Taxon_simple, scales = "free_y") +
  labs(title = "Efectos significativos de ventilación mecánica no invasiva") +
  theme_bw()

ggsave("0_figs/Paciente_longitudinal/boxplot_sign_taxa_vmni.svg",
       width = 10, height = 8)

# GAFAS NASALES DE ALTO FLUJO (GNAF)
splinectome_df$GNAF.durante.elingreso[splinectome_df$GNAF.durante.elingreso == ""] <- NA
splinectome_df$GNAF.durante.elingreso <- factor(splinectome_df$GNAF.durante.elingreso,
                                                levels = c("No", "Sí"))

datos_gnaf <- splinectome_df %>%
  filter(!is.na(rel_abund), !is.na(Tiempo))

gam_model_gnaf <- gam(
  rel_abund ~ GNAF.durante.elingreso +
    s(Tiempo, k = 4, bs = "cr") +
    s(NHC, bs = "re"),
  data   = datos_gnaf,
  method = "REML"
)

summary(gam_model_gnaf)

coefs_gnaf <- summary(gam_model_gnaf)$p.table %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Termino") %>%
  filter(Termino != "s(Tiempo)") %>%
  mutate(
    Termino = case_when(
      Termino == "(Intercept)"          ~ "Intercepto",
      Termino == "GNAF.durante.elingresoSí"      ~ "GNAF.durante.elingreso (Sí vs No)",
      TRUE                              ~ Termino
    )
  ) %>%
  select(
    Termino,
    Estimate    = Estimate,
    Std.Error   = `Std. Error`,
    t.value     = `t value`,
    p.value     = `Pr(>|t|)`
  )

smooths_gnaf <- summary(gam_model_gnaf)$s.table %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Smooth") %>%
  mutate(p.value = `p-value`) %>%
  select(Smooth, edf, Ref.df = `Ref.df`, F = F, p.value)

write_tsv(coefs_gnaf, "0_tables/resultados_gam_gnaf_global.tsv")
write_tsv(smooths_gnaf, "0_tables/resultados_gam_gnaf_temporal.tsv")

df_top$GNAF.durante.elingreso[df_top$GNAF.durante.elingreso == ""] <- NA
df_top$GNAF.durante.elingreso <- factor(df_top$GNAF.durante.elingreso,
                                         levels = c("No", "Sí"))

resultados_gnaf_taxones <- purrr::map_dfr(
  unique(df_top$Taxon_simple),
  function(taxon) {
    datos_taxon <- df_top %>%
      filter(Taxon_simple == taxon, !is.na(GNAF.durante.elingreso)) %>%
      mutate(GNAF.durante.elingreso = factor(GNAF.durante.elingreso, levels = c("No", "Sí")))
    
    n_gnaf <- sum(datos_taxon$GNAF.durante.elingreso == "Sí", na.rm = TRUE)
    n_total  <- nrow(datos_taxon)
    
    if (n_total > 20 & n_gnaf > 2) {
      modelo <- gam(
        rel_abund ~ GNAF.durante.elingreso +
          s(Tiempo, k = 3, bs = "cr") +
          s(NHC, bs = "re"),
        data   = datos_taxon,
        method = "REML",
        family = gaussian()
      )
      
      coef_table <- summary(modelo)$p.table %>%
        as.data.frame() %>%
        tibble::rownames_to_column("Term") %>%
        filter(Term == "GNAF.durante.elingresoSí")
      
      coef_gnaf <- tibble(
        Taxon    = taxon,
        n_total  = n_total,
        n_gnaf   = n_gnaf,
        beta     = coef_table$Estimate,
        se       = coef_table$`Std. Error`,
        t_value  = coef_table$`t value`,
        p_value  = coef_table$`Pr(>|t|)`
      )
      
      assign(paste0("gam_gnaf_", taxon), modelo, envir = .GlobalEnv)
      return(coef_gnaf)
    }
    NULL
  }
)

tabla_gnaf <- resultados_gnaf_taxones %>%
  arrange(p_value) %>%
  mutate(
    p_value_adj = p.adjust(p_value, method = "BH"),
    sig = case_when(
      p_value_adj < 0.05 ~ "***",
      p_value     < 0.05 ~ "**",
      p_value     < 0.1  ~ "*",
      TRUE              ~ ""
    ),
    beta        = round(beta, 4),
    se          = round(se, 4),
    p_value     = round(p_value, 4),
    p_value_adj = round(p_value_adj, 4)
  ) %>%
  select(Taxon, n_total, n_gnaf, beta, se, p_value, p_value_adj, sig)

write_tsv(tabla_gnaf,
          "0_tables/resultados_gam_gnaf_por_taxon.tsv")

# ----------------------------------------------------------------------
# 11. SOM: SELF-ORGANIZING MAP
# ----------------------------------------------------------------------
df_som <- df %>%
  filter(rel_abund >= 0) %>%
  group_by(Taxon) %>%
  filter(sum(rel_abund > 0) >= 5) %>%
  ungroup()

som_mat <- df_som %>%
  select(SampleID, Taxon, rel_abund) %>%
  pivot_wider(
    names_from  = Taxon,
    values_from = rel_abund,
    values_fill = 0
  ) %>%
  column_to_rownames("SampleID")

taxa_ranked <- df_som %>%
  group_by(Taxon) %>%
  summarise(
    total_abund = sum(rel_abund),
    prevalence  = sum(rel_abund > 0),
    .groups     = "drop"
  ) %>%
  mutate(score = total_abund * log1p(prevalence)) %>%
  arrange(desc(score))

n_taxa_som <- min(50, ncol(som_mat))
top_taxa_som <- taxa_ranked$Taxon[1:n_taxa_som]

som_mat_filtered <- som_mat[, top_taxa_som]

som_mat_clean <- som_mat_filtered[complete.cases(som_mat_filtered), ]
som_scaled    <- scale(som_mat_clean)

som_grid <- somgrid(xdim = 4, ydim = 4, topo = "hex")

som_model <- som(
  X         = som_scaled,
  grid      = som_grid,
  rlen      = 500,
  alpha     = c(0.05, 0.01),
  radius    = c(1, 0.5),
  keep.data = TRUE
)

save(som_model, som_scaled, som_mat_clean,
     file = "0_tables/som_model.RData")

codes_scaled <- scale(som_model$codes[[1]])
som_cluster  <- cutree(hclust(dist(codes_scaled)), k = 3)

# Visualizaciones básicas SOM
svg("0_figs/Paciente_longitudinal/som_visualizations.svg",
    width = 10, height = 10)

par(mfrow = c(3, 3))

plot(som_model, type = "counts", main = "Muestras por nodo")
plot(som_model, type = "dist.neighbours", main = "Distancia entre nodos")
plot(
  som_model,
  type  = "mapping",
  bgcol = rainbow(3)[som_cluster],
  main  = "Clusters de nodos"
)
add.cluster.boundaries(som_model, som_cluster)

top6      <- colnames(som_scaled)[1:6]
top6_name <- str_extract(top6, "[^|]+$")

for (i in 1:6) {
  plot(
    som_model,
    type     = "property",
    property = som_model$codes[[1]][, i],
    main     = top6_name[i],
    palette.name = heat.colors
  )
}

dev.off()
par(mfrow = c(1, 1))

write.table(
  som_mat_clean,
  "0_tables/som_matrix.tsv",
  sep       = "\t",
  quote     = FALSE
)

# ----------------------------------------------------------------------
# 12. TABLA SAMPLE–NODE–CLUSTER Y PERFILES TEMPORALES
# ----------------------------------------------------------------------
df_cluster <- tibble(
  SampleID = rownames(som_scaled),
  node     = som_model$unit.classif
) %>%
  mutate(cluster = som_cluster[node])

if ("cluster" %in% names(metadata)) {
  metadata <- metadata %>% rename(cluster_meta = cluster)
}

df_cluster <- df_cluster %>%
  left_join(metadata, by = "SampleID")

top6        <- colnames(som_scaled)[1:6]
top6_pretty <- str_extract(top6, "[^|]+$")

som_df <- som_scaled[, top6, drop = FALSE] %>%
  as.data.frame() %>%
  tibble::rownames_to_column("SampleID")

df_long_som <- df_cluster %>%
  left_join(som_df, by = "SampleID") %>%
  mutate(
    Tiempo_grupo = factor(
      Tiempo_grupo,
      levels = c("T0", "Temprano", "Intermedio", "Tardío")
    )
  ) %>%
  pivot_longer(
    cols      = all_of(top6),
    names_to  = "Taxon_full",
    values_to = "valor_z"
  ) %>%
  mutate(Taxon = str_extract(Taxon_full, "[^|]+$"))

perfiles_multi <- df_long_som %>%
  filter(!is.na(Tiempo_grupo)) %>%
  group_by(Taxon, cluster, Tiempo_grupo) %>%
  summarise(
    mean_val = mean(valor_z, na.rm = TRUE),
    se_val   = sd(valor_z, na.rm = TRUE) / sqrt(n()),
    n        = n(),
    .groups  = "drop"
  )

som_temporal_plot <- ggplot(
  perfiles_multi,
  aes(x = Tiempo_grupo, y = mean_val, group = 1)
) +
  geom_line() +
  geom_point(size = 1.8) +
  geom_errorbar(
    aes(
      ymin = mean_val - se_val,
      ymax = mean_val + se_val
    ),
    width = 0.15
  ) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  facet_grid(Taxon ~ cluster, scales = "free_y") +
  labs(
    x     = "Tiempo",
    y     = "Respuesta SOM (z-score)",
    title = "Perfiles temporales por clúster y taxón (top 6)"
  ) +
  theme_bw()

ggsave(
  "0_figs/Paciente_longitudinal/som_temporal.svg",
  som_temporal_plot,
  width  = 10,
  height = 12
)
