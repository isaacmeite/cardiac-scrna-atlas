# A Newly Curated Single Cell Sequencing Atlas Predicts Efferocytic-dependent Crosstalk in Hearts

[![DOI](https://zenodo.org/badge/1402466168.svg)](https://doi.org/10.5281/zenodo.23112514)

**Running title:** CardioImmunology Mertk Atlas

**Isaac Meite¹#, Lorenzo Pesce², Matthew Feinstein⁵, Edward B. Thorp¹³†, Connor Lantz¹⁴†**

\# First Author | † Co-Corresponding & Co-Senior Authors

**Affiliations**

¹ Department of Pathology, Feinberg School of Medicine, Northwestern University, Chicago, IL, USA  
² Institute for Artificial Intelligence in Medicine, Center for Deep Phenotyping and Precision Therapeutics, Northwestern University, Chicago, IL, USA  
³ Lurie Children's Hospital, Heart Center, Chicago, IL, USA  
⁴ Comprehensive Transplant Center, Department of Surgery, Northwestern University Feinberg School of Medicine, Chicago, IL, USA  
⁵ Indiana University, Division of Cardiovascular Medicine, Department of Medicine, Indiana University School of Medicine, Indianapolis, IN, USA

Published in *American Journal of Physiology - Heart and Circulatory Physiology* (2026)

---

## Overview

This repository contains all analysis code for:

- Construction of a curated 22-cell-type murine cardiac single-cell RNA sequencing atlas integrating 36,121 cells from 8 publicly available datasets
- Training and validation of a cardiac-specific k-nearest neighbor immune classifier (weighted F1 = 0.87 vs 0.32 for ImmGen)
- Pseudotime trajectory analysis of myeloid and stromal compartments using Monocle 3
- Differential abundance and cell-cell communication analysis of MERTK knockout versus wild-type adult hearts at 7 days post-myocardial infarction using CellChat
- Within-cell-type differential expression and pathway module scoring (KO vs WT)

---

## Repository Structure

```
cardiac-scrna-atlas/
├── atlas/
│   └── V2_Atlas_Creation.Rmd              
├── classifier/
│   ├── 01_data_preparation.Rmd
│   ├── 02_harmony_integration.Rmd
│   ├── 03_data_splitting.Rmd
│   ├── 04_smote_balancing.Rmd
│   ├── 05_hyperparameter_tuning.Rmd
│   ├── 06_final_validation.Rmd
│   ├── 07_user_query_classifier.Rmd
│   └── classify_cardiac_immune.R          
└── mertk_analysis/
    ├── MERTK2.Rmd                         
    ├── V2_Adult_Classification.Rmd        
    └── V2_Immgen_v_Classifier.Rmd        
```

---
**Key files:**
- `V2_Atlas_Creation.Rmd` - Atlas integration, annotation, and pseudotime analysis
- `classify_cardiac_immune.R` - Standalone classifier function for deployment
- `MERTK2.Rmd` - Adult-only dataset processing and Harmony integration
- `V2_Adult_Classification.Rmd` - Classifier application and publication figures
- `V2_Immgen_v_Classifier.Rmd` - ImmGen vs cardiac classifier comparison
---


## Data Availability

Published datasets used in atlas construction are available under the following accessions:

- **E-MTAB-13264** - ArrayExpress
- **E-MTAB-9816** - ArrayExpress
- **E-MTAB-7895** - ArrayExpress
- **GSE119355** - NCBI GEO

Raw and processed single-cell RNA sequencing data for the Mertk-deficient and control adult mice generated in this study are deposited in the Gene Expression Omnibus (GEO) - [accession number will be provided at publication].

Supplemental figures and tables are available at:  
https://doi.org/10.6084/m9.figshare.33868957

Large processed objects (integrated Seurat atlas, trained classifier model, CellChat WT/KO objects) will be hosted on Zenodo/Globus. Links will be added to this repository upon publication.

Any additional data are available from the corresponding authors upon reasonable request.

---

## Requirements

R version 4.4.0 or later. Key packages:

```r
Seurat >= 5.3
harmony
SingleCellExperiment
SingleR
celldex
monocle3
CellChat >= 2.0
scuttle
scater
glmGamPoi
ggplot2
patchwork
dplyr
stringr
ggrepel
pheatmap
clustree
```

---

## Usage

Scripts should be run in the following order:

1. `classifier/01_data_preparation.Rmd` through `classifier/06_final_validation.Rmd`
2. `atlas/V2_Atlas_Creation.Rmd`
3. `mertk_analysis/MERTK2.Rmd`
4. `mertk_analysis/V2_Adult_Classification.Rmd`
5. `mertk_analysis/V2_Immgen_v_Classifier.Rmd`

Input paths in each script point to the Northwestern University Quest HPC cluster (`/projects/b1246/Isaac/`). Update these to match your local paths before running.

---

## Citation

If you use this code or atlas in your work, please cite:

> Meite I, Pesce L, Feinstein M, Thorp EB, Lantz C. A Newly Curated Single  
> Cell Sequencing Atlas Predicts Efferocytic-dependent Crosstalk in Hearts.  
> *Am J Physiol Heart Circ Physiol.* 2026. doi: [will be updated at the time of publication]

---

## Contact

**Isaac Meite** — isaacnediembomeite@gmail.com  
**Edward B. Thorp** — ebthorp@northwestern.edu  
**Connor Lantz** — connor.lantz@northwestern.edu
