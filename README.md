# Multi-Omics Analysis of Acute Myeloid Leukemia Using MOFA

## 📌 Description

This project applies **Multi-Omics Factor Analysis (MOFA)** to integrate and analyze multiple biological data types associated with **Acute Myeloid Leukemia (AML)**.

The analysis was performed on **200 AML patient samples**, using three omics datasets:

* **DNA Methylation**
* **RNA-seq**
* **miRNA**

The datasets were obtained from **TCGA through LinkedOmics** and analyzed using **R** and the **MOFA framework**.

MOFA was used to identify **latent factors** that explain major sources of variation across the different omics datasets. The resulting factors were investigated in relation to clinical characteristics such as **disease status, gender, and race**.

## 🎯 Objectives

* Integrate different types of omics data into a unified analysis.
* Identify hidden patterns and major sources of variation.
* Identify important molecular features contributing to specific factors.
* Investigate relationships between molecular factors and clinical characteristics.
* Explore potential sample subgroups and biological patterns.
* Visualize and interpret multi-omics results.

## 🔬 Methodology

The workflow included:

1. Collecting AML multi-omics data from TCGA/LinkedOmics.
2. Preparing DNA methylation, RNA-seq, and miRNA datasets.
3. Applying **Multi-Omics Factor Analysis (MOFA)**.
4. Analyzing the variance explained by different factors.
5. Investigating correlations between factors and clinical variables.
6. Identifying features with high factor weights.
7. Visualizing the results using heatmaps and scatter plots.

MOFA provided a low-dimensional representation of the integrated datasets and helped identify shared and dataset-specific sources of variation.

## 📊 Visualizations

The project includes several visualizations:

* Factor distributions according to clinical characteristics.
* Factor–clinical variable relationships.
* Feature-weight visualizations.
* DNA methylation heatmaps.
* miRNA heatmaps.
* RNA-seq heatmaps.
* Methylation scatter plots.
* miRNA scatter plots.
* RNA-seq scatter plots.

## 🛠️ Technologies

* **R**
* **MOFA (Multi-Omics Factor Analysis)**
* **Bioinformatics**
* **Multi-Omics Data Integration**
* **DNA Methylation Analysis**
* **RNA-seq Analysis**
* **miRNA Analysis**
* **Data Visualization**
* **TCGA / LinkedOmics**

## 🔎 Key Findings

The analysis showed that different latent factors captured variation from different omics modalities. For example, **Factor 1 showed significant variation in DNA methylation**, while **Factor 2 showed moderate RNA-seq variation**, and **Factor 3 showed similar contributions from RNA-seq and miRNA**.

Feature-weight analysis was also used to identify molecular features strongly associated with specific factors across the different omics datasets.

## 🚀 Future Work

Planned extensions include:

* Integrating **mutation data** into the existing multi-omics analysis.
* Creating biological networks for each omics dataset to further investigate relationships between molecular features.
