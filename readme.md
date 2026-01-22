# 2025 Mitochondria Hypoxia Project

## Overview

This repository contains the code and analysis pipeline for the publication:

**"Hypoxia leads to reduced mito-nuclear gene expression and increased mtDNA transcriptional pausing in human cells"**
by Noam Shtolz, Sarah Dadon & Dan Mishmar
*Communications Biology* (2025)
DOI: [10.1038/s42003-025-09457-y](https://doi.org/10.1038/s42003-025-09457-y)

The project investigates how hypoxia affects the coordinated expression of mitochondrial (mtDNA) and nuclear DNA (nDNA)-encoded oxidative phosphorylation (OXPHOS) genes across diverse human cell lines.

---

## Analysis Notebooks

The analysis is structured in sequential notebooks that should be run in order:

### 01_integ_rna_seq_analysis.ipynb
**Purpose:** Data integration and DESeq2 database creation

Processes all RNA-seq datasets and creates a unified DESeq2 results database. Handles 275 RNA-seq samples from 25 datasets representing 18 human cell lines.

**Input:**
- Publicly available RNA-seq data from ENA/GEO (listed in Supplementary Data S1)
- Newly generated RNA-seq data: HCT-116, MCF7, U87 cells (1% O2, 24h)

**Output:**
- Consolidated parquet file containing all DESeq2 results used by subsequent notebooks

---

### 02_mitonuclear_expr_and_hypx_index.ipynb
**Purpose:** Mito-nuclear expression and hypoxia index analysis

Analyzes coordinated mito-nuclear gene expression response to hypoxia.

**Generates:**
- Figure 1a: OXPHOS gene expression across datasets
- Figure 1b: Mito-ribosome gene expression patterns
- Figure 1c: Correlation between hypoxia index and OXPHOS expression
- Figure 1d: Correlation between hypoxia index and mito-ribosome expression

---

### 03_nuclear_response.ipynb
**Purpose:** Nuclear-encoded mitochondrial gene response analysis

Analyzes nuclear-encoded mitochondrial genes and metabolic pathways.

**Generates:**
- Figure 2a: TCA cycle gene expression
- Figure 2b: Fatty acid oxidation pathway response
- Figure 2c: Glycolysis gene expression
- Figure 2d: LDHA vs LDHB expression comparison
- Figure 2e: Correlation between metabolic pathway response and mito-nuclear coordination score

---

### 04_mitochondrial_response.ipynb
**Purpose:** mtDNA-encoded gene response analysis

Analyzes mtDNA-encoded genes and mitochondrial gene regulation pathways.

**Generates:**
- Figure 3a: Mitochondrial transcription genes expression (GO:0006390)
- Figure 3b: Mitochondrial RNA processing genes (GO:0000963)
- Figure 3c: Mitochondrial replication genes (GO:0032042)

---

### 05_sarah_hyp.ipynb
**Purpose:** PRO-seq nascent transcription analysis

Analyzes nascent RNA transcription using PRO-seq data processed through PEPPRO pipeline.

**Generates:**
- Supplementary Figure S11a: Hypoxia marker transcription validation
- Supplementary Figure S11b: nDNA-encoded OXPHOS gene transcription
- Supplementary Figure S11c: nDNA-encoded mito-ribosome gene transcription

**Input:**
- PRO-seq data from HeLa, U87, and D407 cells (1% O2, 24h)
- PEPPRO processed outputs (Supplementary Data S12)

---

### 06_pausing_index_z_method.ipynb
**Purpose:** mtDNA transcriptional pausing site analysis

Analyzes mtDNA transcriptional pausing patterns using PRO-seq data with a Z-score based method.

**Generates:**
- Figure 5a-c: mtDNA PRO-seq read coverage patterns in HeLa, U87, and D407 cells
- Figure 6a: Distribution of pausing sites in heavy and light strands
- Figure 6b-d: Pausing site intensity analysis

---

### 07_peppro_analysis.ipynb
**Purpose:** Additional PRO-seq analysis

Supporting analysis for nascent transcription patterns and quality control of PRO-seq data.

---

## Python Modules

### consts.py
**Purpose:** Project constants and gene definitions

**Defines:**
- `MT_GENES`: mtDNA-encoded genes (13 protein-coding + 2 rRNA)
- `GENELIST_DF`: Master gene list with complex assignments (OXPHOS: 88 nDNA + 13 mtDNA; Mito-ribosome: 78 nDNA + 2 mtDNA)
- `MTDNA_GENE_ORDER`: Gene order along mtDNA
- `COMPLEX_ORDER`: Display order for OXPHOS complexes
- `PEPPRO_RANGES`: Quality control thresholds for PRO-seq analysis

---

### hypoxia_index.py
**Purpose:** Hypoxia index calculation using position-weighted metagene signature

Based on Lombardi et al. (2022) Cell Reports methodology.

**Functions:**
- `generate_hypx_weights(topN)`: Creates position-weighted dictionary from sorted hypoxia genes
- `hypx_index(df, weights)`: Calculates hypoxia index from DESeq2 output

---

### utils.py
**Purpose:** General utility functions for analysis

**Functions:**
- `genelist()` / `genelist_on_folder()`: Load gene lists from files
- Boxplot and strip plot generation with statistical annotations
- Statistical analysis helpers (Mann-Whitney U tests)
- Figure formatting and layout functions

---

### correlation_analysis.py
**Purpose:** Statistical correlation and mito-nuclear score calculations

**Functions:**
- `extract_resids()`: Extract regression residuals
- `correlation_analysis()`: Spearman correlation analysis
- Mito-nuclear score quantification (captures proportion of coordinated DEGs)
- Pathway enrichment statistics

---

### annotation.py
**Purpose:** Gene annotation and GenBank sequence processing

**Functions:**
- `annotate()`: Creates DataFrame of gene features and positions from GenBank records
- `strand_coverage()`: Per-strand coverage analysis for mtDNA
- `TPM_calc()`: Transcripts per million calculation
- `record_check()`: GenBank file management and retrieval from NCBI

---

### remove_batch_effects.py
**Purpose:** Batch effect correction using ComBat algorithm

**Functions:**
- `retrieve_batches()`: Codify batch information from sample metadata
- `remove_batch_effects()`: Wrapper for pycombat batch correction
- `qc_plots()`: Pre/post correction visualization

---

### run_deseq.py
**Purpose:** R/DESeq2 interface wrapper using diffexpr package

**Functions:**
- `run_deseq()`: Calls DESeq2 for differential expression analysis
- Handles contrasts and generates normalized counts
- Returns log2 fold-change and adjusted p-value tables

---

### pausing_index_calculator.py
**Purpose:** Calculate mtDNA transcriptional pausing index from PRO-seq data

Uses Z-score method to identify pausing sites as described in the manuscript.

**Functions:**
- `mtdna_region()`: Circular DNA window definition handling mtDNA boundaries
- `z_score()`: Standardization function
- `iter_over_windows()`: Rolling Z-score calculation across mtDNA positions

---

### proseq_mito_peak_caller.py
**Purpose:** Locate transcription initiation and termination sites from PRO-seq data

**Functions:**
- `pro_peaks()`: Identifies transcription initiation sites (TIS) and termination sites (TERM)
- Handles logging and circular mtDNA geometry

---

### flanking_regions.py
**Purpose:** Filter annotation output by flanking regions for circular mtDNA

**Functions:**
- `genomic_ranges()`: Extract range coordinates from annotations
- `wrap()`: Handle mtDNA circular boundaries for coordinate calculations

---

### dna_features_translator_class.py
**Purpose:** Custom visualization classes for GenBank mtDNA features

**Classes:**
- `MitoTranslator`: Custom translator for dna_features_viewer with color coding by OXPHOS complex (I-V)
- Feature filtering and labeling for publication-quality mtDNA maps

---

### scale_axis.py
**Purpose:** Custom asymmetric scale for matplotlib axes

**Classes:**
- `AsymScale`: Custom axis transformation allowing different scaling for positive and negative values

---

## Data Acquisition

### RNA-seq Data
Publicly available RNA-seq datasets were retrieved from ENA (www.ebi.ac.uk/ena) and GEO (ncbi.nlm.nih.gov/geo/). The complete list of datasets with accession numbers is available in Supplementary Data S1 on FigShare.

### PRO-seq Data
PRO-seq data generated in this study has been deposited in ArrayExpress under accession E-MTAB-14981.

### Supplementary Data Files
All supplementary data files are available at:
[https://figshare.com/projects/Shtolz_2025_Mitochondria_Hypoxia_Project/237023](https://figshare.com/projects/Shtolz_2025_Mitochondria_Hypoxia_Project/237023)

**Key files:**
- **S1**: RNA-seq dataset metadata
- **S2**: Hypoxia marker genes (48 genes with HN-scores from Lombardi et al. 2022)
- **S3**: OXPHOS gene lists by complex
- **S4**: Mito-ribosome gene lists
- **S9**: DESeq2 results (log2 fold-change and p-values)
- **S10**: DESeq2 normalized expression values
- **S11**: qPCR raw CT values for mtDNA copy number
- **S12**: PEPPRO PRO-seq analysis outputs

---

## Installation

### Prerequisites
```
Python 3.8+
Jupyter Notebook
R 4.0+ (for DESeq2)
```

### Install Python Dependencies
```bash
pip install -r requirements.txt
```

---

## Citation

```bibtex
@article{shtolz2025hypoxia,
  title={Hypoxia leads to reduced mito-nuclear gene expression and increased
         mtDNA transcriptional pausing in human cells},
  author={Shtolz, Noam and Dadon, Sarah and Mishmar, Dan},
  journal={Communications Biology},
  year={2025},
  doi={10.1038/s42003-025-09457-y}
}
```

---

## Contact

**Corresponding Author:**
Dan Mishmar, PhD
Department of Life Sciences
Ben-Gurion University of the Negev
Email: dmishmar@bgu.ac.il

---

## License

This article is licensed under a Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International License.
