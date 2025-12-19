# Caracterización de la microbiota intestinal de pacientes críticos y su relación con el estado de salud

Este repositorio contiene el proyecto desarrollado como parte del Trabajo de Fin de Máster en Bioinformática y Biología Computacional en la Universidad Autónoma de Madrid. El objetivo general del proyecto es entender la relación de la microbiota intestinal con el desarrollo y progresión de las neumonías intrahospitalarias en la UCI. Los objetivos específicos del TFM son:

- Analizar la evolución de la composición y diversidad de la microbiota intestinal de pacientes críticos seguidos longitudinalmente. La caracterización se centrará en el análisis de datos de secuenciación de 16S rRNA (bacterias y arqueas).
- Estimar las principales interacciones microbianas de la microbiota a través del análisis de redes de cooccurrencia.
- Determinar la asociación de la microbiota (diversidad y composición) con variables clínicas recogidas de los pacientes para entender el impacto de la microbiota en la evolución clínica (metadata). 

## Estructura de trabajo

```
├── 0_figs/                  # Figuras generadas a lo largo del pipeline
├── 0_tables/                # Tablas generadas a lo largo del pipeline
├── 1_input/                 # Metadatos y asignación taxonómica realizada con QIIME2
├── 2_fastspar/              # Estimación de la correlación de datos composicionales con FastSpar
│   ├── bootstrap_corr/      # Correlación de bootstrap para FastSpar
│   └── bootstrap_counts/    # Conteo de bootstrap para FastSpar
├── 3_cooccurrence_network/  # Resultados de la red de coocurrencia de CooccurrenceAffinity R
│   ├── affinity.tsv        
│   ├── alpha_mle.tsv     
│   ├── alpha_sign_selected.tsv
│   ├── OTU_species_aggregated_top.tsv
│   ├── pvalues.tsv
│   ├── top_PA.tsv
│   └── net.sp.tsv
├── environment-cograph.yml  # Entorno conda para FastSpar
├── run_fastspar.sh
├── script_adonis.R
├── script_diff_abund_ancom.R
├── script_diversity_clinical_vars.R
├── script_longitudinal_analysis.R
├── script_network_analysis.R
├── script_patient_analysis.R
└── README.md
```

El orden en el que se han creado y ejecutado los scripts, y el orden en el que se deberían leer, es:

- `script_adonis.R`
- `script_patient_analysis.R`
- `script_diversity_clinical_vars.R`
- `script_longitudinal_analysis.R`
- `script_diff_abund_ancom.R`
- `script_network_analysis.R` (secciones 1-4)
- `run_fastspar.sh`
- `script_network_analysis.R` (sección 5)

## Requisitos y entornos

El trabajo se ha realizado en Linux Mint 22.2 Cinnamon mediante el uso de entornos de **conda** y en R 4.3.3. En el repositorio se incluye el fichero `yml` con los requisitos concretos para FastSpar. 

```
conda env create -f environment-cograph.yml
conda activate cograph
```

Los paquetes de R necesarios para ejecutar los scripts se encuentran especificados al inicio de cada uno.
