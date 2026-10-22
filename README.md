# Bee-Nest Soil Microbiome

This repository contains the R script and dataset used to analyse physicochemical and microbial properties of nest-entrance soil, control soil, and bees in Lasioglossum malachurum (2023).

---

## Associated publication

**Title:** *Overlooked soil engineers: digger bees induce localised shifts in soil chemical properties and microbiome*

**Authors:** *Andrea Ferrari, Raphael Tiziani, Alessandro Borsotti, Federico Ronchetti, Ilaria Fracasso, Oussama Bouaicha, Tanja Mimmo, Luigimaria Borruso, Carlo Polidori

**Journal:** *Soil Biology and Biochemistry*

**Year:** *2026*

**DOI:** *https://doi.org/10.1016/j.soilbio.2026.110341*

**Repository:** *Raw dataset is available as supplementary file*

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

* **0 – Prepare the environment** → Installs and loads required packages, creates output directories, imports the dataset and metadata, formats variables, and prepares the R environment.

* **1 – Soil chemistry analyses** → Performs descriptive and inferential analyses of soil physicochemical variables.

* **2 – Bacterial community analyses (16S)** → Characterizes bacterial communities associated with bee nests, bees, and control soils.

* **3 – Fungal community analyses (ITS)** → Characterizes fungal communities associated with bee nests, bees, and control soils.

* **4 – Supplementary outputs** → Generates additional tables summarizing taxonomic occurrence across sample types.

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

**Andrea Ferrari** *a.ferrari.research@gmail.com*

---

## Citation

If you use this code or dataset, please cite the associated publication.
