# Bee-Nest Soil Microbiome

This repository contains the R script and dataset used to analyse soil chemistry, bacterial and fungal communities associated with bee nests, surrounding soils, and bees using physicochemical measurements and metabarcoding data.

---

## Associated publication

**Title:** *To be added upon publication*

**Authors:** *To be added*

**Journal:** *Under review / in preparation*

**Year:** *2026*

**DOI:** *To be added*

**Repository:** *To be added (e.g. Dryad or Zenodo)*

---

## Necessary files to run the analyses

* `Script.R` → Main script for the complete analytical workflow
* `Dataset_Lmal2023.xlsx` → Dataset containing metadata, soil chemistry, bacterial (16S), and fungal (ITS) sequencing data

**Important:** `Script.R` and `Dataset_Lmal2023.xlsx` must be placed in the **same directory**.

When the script is executed, the following folders will be automatically created:

* `results/` → Statistical outputs, model summaries, tables, and supplementary datasets
* `figures/` → Publication-ready figures and graphical outputs

---

## Workflow overview

1. Download the repository files
2. Place `Script.R` and `Dataset_Lmal2023.xlsx` in the same folder
3. Open the script in RStudio
4. Run the script sequentially following the numbered sections

The analysis pipeline is organized into modular sections.

---

## Script structure

### **0 – Prepare the environment**

Installs and loads all required R packages, creates output directories, imports metadata and sequencing datasets, formats variables, defines plotting palettes, and prepares the working environment.

---

### **1 – Soil chemistry analyses**

Performs descriptive and inferential analyses of soil physicochemical variables.

Analyses include:

* Summary statistics (mean ± SD)
* Linear mixed-effects models
* Type-II ANOVA
* Tukey post-hoc comparisons
* Principal Component Analysis (PCA)
* PERMANOVA
* Multivariate dispersion analyses
* Coefficient of variation comparisons
* Publication-quality boxplots and PCA visualizations

---

### **2 – Bacterial community analyses (16S)**

Characterizes bacterial communities associated with bee nests, bees, and control soils.

Analyses include:

* Relative abundance barplots at multiple taxonomic levels
* Heatmaps of dominant bacterial genera
* Alpha-diversity analyses (Richness, Shannon diversity, Pielou evenness)
* Linear mixed-effects models
* Differential abundance analysis using Maaslin2
* Venn diagrams of shared OTUs
* NMDS ordination
* Beta-diversity analyses
* PERMANOVA and dispersion tests

---

### **3 – Fungal community analyses (ITS)**

Characterizes fungal communities associated with bee nests, bees, and control soils.

Analyses include:

* Relative abundance barplots at multiple taxonomic levels
* Heatmaps of dominant fungal genera
* Alpha-diversity analyses
* Linear mixed-effects models
* Differential abundance analysis using Maaslin2
* Venn diagrams of shared OTUs
* NMDS ordination
* Beta-diversity analyses
* PERMANOVA and dispersion tests

---

### **4 – Supplementary outputs**

Generates additional tables summarizing taxonomic occurrence across sample types.

Outputs include:

* Presence/absence tables for bacterial genera
* Presence/absence tables for fungal genera
* Shared taxonomic information among bees, nests, and surrounding soils

---

## Main outputs

### Results

* Summary statistics
* Linear mixed model outputs
* ANOVA tables
* Post-hoc pairwise comparisons
* PERMANOVA results
* Alpha-diversity metrics
* Beta-diversity statistics
* Maaslin2 differential abundance results
* Shared taxonomy tables

### Figures

* Soil chemistry boxplots
* PCA ordinations
* Relative abundance barplots
* Heatmaps
* Alpha-diversity boxplots
* NMDS ordinations
* Dispersion plots
* Venn diagrams
* Maaslin2 effect size plots

---

## Requirements

The script automatically installs missing CRAN packages.

The following Bioconductor package must also be available:

* `DESeq2`

If necessary, install it manually using:

```r
BiocManager::install("DESeq2")
```

The script was developed using **R version > 4.2**, although recent versions of R are recommended.

---

## Contact

For questions regarding the analyses, data, or code:

**Andrea Ferrari**

*a.ferrari.research@gmail.com*

---

## Citation

If you use this code or dataset, please cite the associated publication once it becomes available.
