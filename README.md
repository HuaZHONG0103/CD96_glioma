This repository contains the complete, modular R bioinformatics pipeline designed to **reproduce the computational analyses and publication-grade figures** from the published study:

> **Zhang Q, Zhong H, Fan Y, Liu Q, Song J, Yao S, Cao F.** (2020). *Immune and Clinical Features of CD96 Expression in Glioma by in silico Analysis*. *Frontiers in Bioengineering and Biotechnology*, 8, 592. [DOI: 10.3389/fbioe.2020.00592](https://doi.org/10.3389/fbioe.2020.00592)

---

## 📊 Project Overview

Glioma is the most common primary malignant brain tumor. This project investigates the role of the immune checkpoint molecule **CD96** as an immune-context biomarker linked to glioma malignancy, immune infiltration, and patient prognosis. Utilizing public cohorts from **TCGA** and the **Chinese Glioma Genome Atlas (CGGA)**, the pipeline covers:
- **Clinical Phenotype Association**: Examining $CD96$ expression across WHO grades, IDH mutation status, and molecular subtypes.
- **Biomarker & ROC Performance**: Evaluating diagnostic capacity using ROC and AUC metrics (`pROC`).
- **Functional Enrichment**: Performing Gene Ontology (GO) and GSVA pathway analysis on $CD96$-correlated genes.
- **Immune Microenvironment**: Analyzing immune checkpoint co-expression networks and immune cell infiltration (TIMER).
- **Prognostic Modeling**: Constructing Kaplan-Meier survival curves and multivariate Cox proportional hazards models.
