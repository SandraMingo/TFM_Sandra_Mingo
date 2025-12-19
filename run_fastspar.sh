#!/bin/bash
# ======================================================================
# FASTSPAR: CORRELACIONES COMPOSICIONALES PARA MICROBIOMA (SparCC)
# ----------------------------------------------------------------------
# Descripción:
#   1. Calcula correlaciones y covarianzas en la tabla de abundancias
#      original usando FastSpar (implementación optimizada de SparCC).
#   2. Genera 100 tablas bootstrap de la matriz de OTUs colapsada a género.
#   3. Calcula matrices de correlación para cada bootstrap (500 iteraciones).
#   4. Estima p-valores pseudo mediante comparación con distribuciones
#      bootstrap (100 permutaciones).
#
# SparCC corrige el sesgo composicional inherente a datos de secuenciación
# de microbioma, proporcionando estimaciones de correlación más fiables
# que Pearson/Spearman para tablas de abundancia relativa.
#
# Entradas:
#   - 0_tables/abundance_matrix.tsv (generada por script R de redes)
#
# Salidas:
#   - 2_fastspar/correlation.tsv
#   - 2_fastspar/covariance.tsv
#   - 2_fastspar/bootstrap_counts/boot_*.tsv (100 tablas)
#   - 2_fastspar/bootstrap_corr/cor_boot_*.tsv (100 correlaciones)
#   - 2_fastspar/pvalues.tsv
#
# Requisitos:
#   - Conda: cograph
#   - FastSpar instalado
#   - 0_tables/abundance_matrix.tsv (géneros, formato BIOM)
#
# Autor: Sandra Mingo-Ramirez
# Fecha: 2025
# ======================================================================

set -euo pipefail  # Exit on error, undefined vars, pipe failures

# Configuración
PROJECT_DIR="/home/sandra/Projects/UCI_Variables/"
cd "$PROJECT_DIR"

# Verificar entorno conda
if ! conda info --envs | grep -q "^cograph\s"; then
    echo "ERROR: Activar entorno 'cograph' primero: conda activate cograph"
    exit 1
fi

echo "=== FASTSPAR: Análisis de correlaciones composicionales ==="
echo "Directorio: $PROJECT_DIR"
echo "Fecha: $(date)"

# Verificar tabla de entrada
if [[ ! -f "0_tables/abundance_matrix.tsv" ]]; then
    echo "ERROR: 0_tables/abundance_matrix.tsv no existe."
    echo "Ejecutar primero script R de redes de co-ocurrencia."
    exit 1
fi

# Crear directorios de salida
mkdir -p 2_fastspar{,/bootstrap_counts,/bootstrap_corr}

echo "[1/4] Estimando correlaciones/covarianzas en tabla original..."
fastspar \
    --otu_table 0_tables/abundance_matrix.tsv \
    --correlation 2_fastspar/correlation.tsv \
    --covariance 2_fastspar/covariance.tsv \
    --iterations 500

echo "[2/4] Generando 100 tablas bootstrap..."
fastspar_bootstrap \
    --otu_table 0_tables/abundance_matrix.tsv \
    --number 100 \
    --prefix 2_fastspar/bootstrap_counts/boot_ \
    --threads 4

echo "[3/4] Calculando correlaciones para bootstraps (100 archivos)..."
i=0
for f in 2_fastspar/bootstrap_counts/boot_*.tsv; do
    if [[ -f "$f" ]]; then
        echo "Bootstrap $i/100: $f"
        fastspar \
            --otu_table "$f" \
            --correlation "2_fastspar/bootstrap_corr/cor_boot_${i}.tsv" \
            --covariance "2_fastspar/bootstrap_corr/cov_boot_${i}.tsv" \
            --iterations 500 \
            --threads 4
        i=$((i + 1))
    fi
done

echo "[4/4] Calculando p-valores pseudo (100 permutaciones)..."
fastspar_pvalues \
    --otu_table 0_tables/abundance_matrix.tsv \
    --correlation 2_fastspar/correlation.tsv \
    --prefix 2_fastspar/bootstrap_corr/cor_boot_ \
    --permutations 100 \
    --outfile 2_fastspar/pvalues.tsv

echo "=== ANÁLISIS FASTSPAR COMPLETADO ==="
echo "Archivos generados:"
echo "  - 2_fastspar/correlation.tsv      (correlaciones originales)"
echo "  - 2_fastspar/covariance.tsv       (covarianzas originales)"
echo "  - 2_fastspar/pvalues.tsv          (p-valores pseudo)"
echo "  - 2_fastspar/bootstrap_counts/    (100 tablas bootstrap)"
echo "  - 2_fastspar/bootstrap_corr/      (100 matrices correlación)"
echo ""
echo "Próximo paso: ejecutar script R para visualizar redes."
echo "Listo para análisis downstream con igraph/ggraph."

