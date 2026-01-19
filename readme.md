# Hypoxia Effects on Mito-Nuclear Coordination

## Overview

This repository contains the code and analysis pipeline for the publication:

**"Hypoxia leads to reduced mito-nuclear gene expression and increased mtDNA transcriptional pausing in human cells"**
by Noam Shtolz, Sarah Dadon & Dan Mishmar
*Communications Biology* (2025)
DOI: [10.1038/s42003-025-09457-y](https://doi.org/10.1038/s42003-025-09457-y)

The project investigates how hypoxia affects the coordinated expression of mitochondrial (mtDNA) and nuclear DNA (nDNA)-encoded oxidative phosphorylation (OXPHOS) genes across diverse human cell lines.

## Key Findings

1. **Coordinated MNGE reduction**: Hypoxia induces coordinated down-regulation of mito-nuclear OXPHOS gene expression in 78.5% (11/14) of responsive cell line datasets
2. **RNA-level regulation**: The reduction occurs at the RNA level, not through changes in mtDNA copy number
3. **HIF1α involvement**: HIF1α is essential for coordinated MNGE response, but mtDNA genes respond to hypoxia even in HIF1α knockout cells
4. **Metabolic shift**: Broader mitochondrial metabolic response including TCA cycle and fatty acid oxidation downregulation, with glycolysis upregulation
5. **mtDNA transcriptional pausing**: Hypoxia increases mtDNA transcriptional pausing sites and intensity, suggesting a novel regulatory mechanism

## Analysis Workflow

The analysis is structured in **5 sequential notebooks** that should be run in order:

### 01. Data Integration and DESeq2 Database Creation
**Notebook:** [01_integ_rna_seq_analysis.ipynb](01_integ_rna_seq_analysis.ipynb)

**Purpose:** Processes all RNA-seq datasets and creates a unified DESeq2 results database

**Input:**
- 275 RNA-seq samples from 25 datasets representing 18 human cell lines
- Publicly available data from ENA/GEO (see Supplementary Data S1)
- Newly generated data: HCT-116, MCF7, U87 cells (1% O₂, 24h)

**Processing:**
- DESeq2 differential expression analysis for each dataset
- Identifies "responder" datasets (≥25% DEGs in OXPHOS genes)
- Validates hypoxia response using marker genes

**Output:**
- Consolidated parquet file containing all DESeq2 results
- This database is used by all subsequent analysis notebooks

**Note:** This notebook contains exploratory code; the main purpose is database creation.

---

### 02. Mito-Nuclear Expression and Hypoxia Index Analysis
**Notebook:** [02_mitonuclear_expr_and_hypx_index.ipynb](02_mitonuclear_expr_and_hypx_index.ipynb)

**Purpose:** Analyzes coordinated mito-nuclear gene expression response to hypoxia

**Generates:**
- **Figure 1a**: OXPHOS gene expression across datasets (coordinated down-regulation in 11/14 datasets)
- **Figure 1b**: Mito-ribosome gene expression patterns
- **Figure 1c**: Correlation between hypoxia index and OXPHOS expression (ρ=-0.63, p=0.0008)
- **Figure 1d**: Correlation between hypoxia index and mito-ribosome expression (ρ=-0.69, p=0.0001)

**Key analyses:**
- Hypoxia index calculation using Lombardi et al. (2022) position-weighted metagene signature
- Mito-nuclear score quantification
- Classification of responsive vs non-responsive cell lines

**Depends on:** DESeq2 database from notebook 01

---

### 03. Nuclear-Encoded Gene Response
**Notebook:** [03_nuclear_response.ipynb](03_nuclear_response.ipynb)

**Purpose:** Analyzes nuclear-encoded mitochondrial genes and metabolic pathways

**Generates:**
- **Figure 2a**: TCA cycle gene expression (coordinated down-regulation)
- **Figure 2b**: Fatty acid oxidation pathway response (down-regulation)
- **Figure 2c**: Glycolysis gene expression (up-regulation - metabolic shift)
- **Figure 2d**: LDHA vs LDHB expression comparison (LDHA upregulated, LDHB not)
- **Figure 2e**: Correlation between metabolic pathway response and mito-nuclear coordination score

**Key analyses:**
- GO pathway enrichment (TCA cycle, fatty acid oxidation, glycolysis)
- Metabolic shift characterization
- Relationship between metabolic response and MNGE coordination

**Depends on:** DESeq2 database from notebook 01

---

### 04. Mitochondrial DNA-Encoded Gene Response
**Notebook:** [04_mitochondrial_response.ipynb](04_mitochondrial_response.ipynb)

**Purpose:** Analyzes mtDNA-encoded genes and mitochondrial gene regulation pathways

**Generates:**
- **Figure 3a**: Mitochondrial transcription genes expression (GO:0006390) - down-regulated
- **Figure 3b**: Mitochondrial RNA processing genes (GO:0000963) - down-regulated
- **Figure 3c**: Mitochondrial replication genes (GO:0032042) - down-regulated but less consistent

**Key analyses:**
- mtDNA-encoded OXPHOS gene expression patterns
- Analysis of genes involved in mtDNA transcription, processing, and replication
- Evidence that MNGE reduction occurs at RNA level, not through mtDNA copy number changes

**Depends on:** DESeq2 database from notebook 01

**Note:** Figure 3d-g (qPCR mtDNA copy number) code will be added later

---

### 05. PRO-seq Nascent Transcription Analysis
**Notebook:** [05_peppro_analysis.ipynb](05_peppro_analysis.ipynb)

**Purpose:** Analyzes nascent RNA transcription using PRO-seq data

**Generates:**
- **Supplementary Figure S11b**: nDNA-encoded OXPHOS gene transcription (no significant change)
- **Supplementary Figure S11c**: nDNA-encoded mito-ribosome gene transcription (no significant change)

**Key analyses:**
- PEPPRO pipeline outputs for nuclear gene transcription
- Demonstrates that nDNA OXPHOS downregulation is post-transcriptional
- Hypoxia marker transcription validation (Supplementary Fig S11a)

**Input data:**
- PRO-seq data from HeLa, U87, and D407 cells (1% O₂, 24h)
- Requires PEPPRO processed outputs (Supplementary Data S12)

**Note:** Full PRO-seq analysis including Figure 5 (mtDNA transcriptional patterns) and Figure 6 (pausing site analysis) will be added later

---

## Python Scripts

### Core Analysis Modules

#### **utils.py**
General utility functions used across all notebooks

**Key functions:**
- `genelist()` / `genelist_on_folder()`: Load gene lists from files
- Custom boxplot and strip plot generation
- Statistical analysis helpers (Mann-Whitney U tests)
- Figure formatting and layout functions

#### **consts.py**
Project constants and gene definitions

**Defines:**
- `MT_GENES`: mtDNA-encoded genes (13 protein-coding + 2 rRNA)
- `GENELIST_DF`: Master gene list with complex assignments
  - OXPHOS genes: 88 nDNA + 13 mtDNA
  - Mito-ribosome genes: 78 nDNA + 2 mtDNA
  - Complex assignments: I, II, III, IV, V
- `MTDNA_GENE_ORDER`: Gene order along mtDNA
- `COMPLEX_ORDER`: Display order for complexes
- `PEPPRO_RANGES`: Quality control thresholds

#### **hypoxia_index.py**
Hypoxia index calculation using position-weighted metagene signature

**Functions:**
- `generate_hypx_weights(topN)`: Creates position-weighted dictionary from sorted hypoxia genes
- `hypx_index(df, weights)`: Calculates hypoxia index from DESeq2 output

**Method:** Based on Lombardi et al. (2022) - uses top N hypoxia marker genes with position-based weighting. Paper uses top 48 genes.

#### **correlation_analysis.py**
Statistical correlation and mito-nuclear score calculations

**Functions:**
- Spearman correlation analysis
- Mito-nuclear score calculation (captures proportion of coordinated DEGs)
- Pathway enrichment statistics

#### **annotation.py**
Gene annotation and name standardization utilities

**Functions:**
- Gene symbol standardization
- GO term associations
- Complex and pathway membership assignments

#### **remove_batch_effects.py**
Batch effect correction and normalization

**Functions:**
- Dataset normalization procedures
- Batch effect assessment
- Combat correction (if applicable)

#### **run_deseq.py**
R/DESeq2 interface wrapper

**Purpose:**
- Calls R-based DESeq2 analysis
- Generates log2 fold-change and adjusted p-value tables

---

## Data Files

### Download from FigShare
All supplementary data files are available at:
[https://figshare.com/projects/Shtolz_2025_Mitochondria_Hypoxia_Project/237023](https://figshare.com/projects/Shtolz_2025_Mitochondria_Hypoxia_Project/237023)

**Key Data Files:**
- **S1**: RNA-seq dataset metadata (275 samples, 25 datasets, 18 cell lines, hypoxia conditions)
- **S2**: Hypoxia marker genes from Lombardi et al. (2022) - 48 genes with HN-scores
- **S3**: OXPHOS gene lists (by complex)
- **S4**: Mito-ribosome gene lists
- **S9**: DESeq2 results - log2 fold-change and p-values for all genes across all datasets
- **S10**: DESeq2 normalized expression values
- **S11**: qPCR raw CT values for mtDNA copy number measurements
- **S12**: PEPPRO PRO-seq analysis outputs (gene-level transcription quantification)

### Gene Lists Directory
Located in `data/genes/` (TSV format with Name, Complex, Description, Protein_type columns):

**OXPHOS genes:**
- Complex I (NADH dehydrogenase): nDNA and mtDNA subunits
- Complex II (Succinate dehydrogenase): nDNA only
- Complex III (Cytochrome bc1): nDNA and mtDNA subunits
- Complex IV (Cytochrome c oxidase): nDNA and mtDNA subunits
- Complex V (ATP synthase): nDNA and mtDNA subunits

**Other gene lists:**
- Mito-ribosome genes (MRPS*, MRPL*)
- Mitochondrial transcription (GO:0006390)
- Mitochondrial RNA processing (GO:0000963)
- Mitochondrial replication (GO:0032042)
- TCA cycle (GO:0006099)
- Fatty acid oxidation (GO:0019395)
- Glycolysis (GO:0006110)
- Mitophagy (GO:0000423)

### Hypoxia Markers
Located in `data/hyp_genes/`:
- `human18_hn2.tsv`: 48 hypoxia marker genes with HN-scores from Lombardi et al. (2022)

---

## Installation and Setup

### Prerequisites
```bash
Python 3.8+
Jupyter Notebook
R 4.0+ (for DESeq2)
```

### Install Python Dependencies
```bash
pip install -r requirements.txt
```

**Required packages:**
```
pandas >= 1.3.0
numpy >= 1.21.0
scipy >= 1.7.0
seaborn >= 0.11.0
matplotlib >= 3.4.0
statannot >= 0.2.3
statannotations >= 0.5.0
```

### Download Data
1. Download supplementary data from [FigShare](https://figshare.com/projects/Shtolz_2025_Mitochondria_Hypoxia_Project/237023)
2. Place data files in the `data/` directory

---

## Running the Analysis

### Quick Start
Run notebooks in sequential order:

```bash
jupyter notebook 01_integ_rna_seq_analysis.ipynb  # Creates DESeq2 database
jupyter notebook 02_mitonuclear_expr_and_hypx_index.ipynb  # Figures 1
jupyter notebook 03_nuclear_response.ipynb  # Figure 2
jupyter notebook 04_mitochondrial_response.ipynb  # Figure 3a-c
jupyter notebook 05_peppro_analysis.ipynb  # Supplementary Figure S11
```

### Example Usage

**Calculate hypoxia index:**
```python
from hypoxia_index import generate_hypx_weights, hypx_index
import pandas as pd

# Load DESeq2 results
deseq_df = pd.read_parquet('deseq_results.parquet')

# Generate weights for top 48 genes (as used in paper)
weights = generate_hypx_weights(topN=48)

# Calculate hypoxia index
index = hypx_index(deseq_df, weights,
                   gene_col='gene',
                   log2fc_col='log2FoldChange')
print(f"Hypoxia index: {index:.3f}")
```

**Load gene lists:**
```python
from consts import GENELIST_DF, MT_GENES

# Get all Complex I genes
complex_i = GENELIST_DF[GENELIST_DF['Complex'] == 'I']

# Separate by genome
nuc_genes = complex_i[complex_i['Genome'] == 'NUC']
mt_genes = complex_i[complex_i['Genome'] == 'MT']

print(f"Complex I: {len(nuc_genes)} nuclear + {len(mt_genes)} mitochondrial genes")
```

**Generate publication plots:**
```python
import utils
import pandas as pd

# Load expression data
expr_df = pd.read_parquet('normalized_expression.parquet')

# Create boxplot with statistical annotations
utils.make_boxplot_with_stats(
    data=expr_df,
    x='condition',
    y='expression',
    hue='genome',
    pairs=[('Normoxia', 'Hypoxia')]
)
```

---

## Experimental Methods Summary

### RNA-seq Data Analysis
- **Datasets**: 18 cell lines, 25 datasets, 275 samples
- **Hypoxia conditions**: 0.2-2% O₂ for 16-24 hours
- **Alignment**: STAR aligner v2.7.11 to GRCh38
- **Counting**: HTSeq v0.6.1
- **Differential expression**: DESeq2 v1.46
- **Significance threshold**: FDR-adjusted p-value < 0.1
- **Responder criteria**: ≥25% DEGs in OXPHOS genes (22/88 nDNA or 3/13 mtDNA)

### PRO-seq Analysis
- **Cell lines**: HeLa, U87, D407
- **Condition**: 1% O₂ for 24 hours
- **Method**: Biotin-11-CTP run-on labeling
- **Sequencing**: Paired-end 150bp (Illumina NextSeq)
- **Pipeline**: PEPPRO v0.10.1 for nDNA genes
- **Custom analysis**: mtDNA transcription patterns and pausing sites

### mtDNA Copy Number (qPCR)
- **Cell lines**: HeLa, U87, MCF7, HCT-116
- **Reference genes**: JUN and ANKRD37 (or ACTB)
- **mtDNA targets**: 4 regions (RNR2, CO3, CO1, ND5)
- **Replicates**: 4 biological replicates, technical triplicates
- **Calculation**: Relative copy number = 2 × 2^ΔCT

---

## Key Results Summary

### Coordinated MNGE Response
- **11/14 datasets (78.5%)** show coordinated down-regulation of OXPHOS genes
- **Strong correlation** between hypoxia index and MNGE (ρ=-0.63 for OXPHOS, ρ=-0.69 for ribosome)
- Response is **genome-wide**, affecting both mtDNA and nDNA-encoded genes

### Metabolic Shift
- **TCA cycle**: Down-regulated
- **Fatty acid oxidation**: Down-regulated (correlated with mito-nuclear score, ρ=-0.58)
- **Glycolysis**: Up-regulated
- **LDHA/LDHB ratio**: Increased (shift to glycolytic metabolism)

### Mechanism
- **RNA-level regulation**: No significant mtDNA copy number changes
- **HIF1α-dependent**: Essential for coordinated nDNA OXPHOS response
- **HIF1α-independent**: mtDNA genes respond even in HIF1α knockout cells
- **Post-transcriptional**: nDNA OXPHOS transcription unchanged (PRO-seq), suggesting post-transcriptional regulation

### Candidate Regulators
- **MTERF3** and **LONP1**: HIF1α targets differentially expressed in responders
- **c-Jun (JUN)**: HIF1α-independent candidate for mtDNA regulation
- **KLF2/ETV1**: Mitochondrial biogenesis repressors correlated with MNGE response

---

## Repository Structure

```
.
├── README.md                                    # This file
├── requirements.txt                             # Python dependencies
├── RELEASE_v1.0.0.md                           # Release notes
│
├── Analysis Notebooks (run in order)
│   ├── 01_integ_rna_seq_analysis.ipynb         # DESeq2 database creation
│   ├── 02_mitonuclear_expr_and_hypx_index.ipynb # Figures 1
│   ├── 03_nuclear_response.ipynb               # Figure 2
│   ├── 04_mitochondrial_response.ipynb         # Figure 3a-c
│   └── 05_peppro_analysis.ipynb                # Supplementary Fig S11
│
├── Core Python Modules
│   ├── utils.py                                # General utilities
│   ├── consts.py                               # Constants and gene lists
│   ├── hypoxia_index.py                        # Hypoxia index calculation
│   ├── correlation_analysis.py                 # Statistical analysis
│   ├── annotation.py                           # Gene annotations
│   ├── remove_batch_effects.py                 # Batch correction
│   └── run_deseq.py                           # DESeq2 interface
│
├── Exploratory Notebooks (not in paper)
│   └── analyze_vhl_deficiency.ipynb           # VHL analysis
│
└── data/                                       # Data directory
    ├── genes/                                  # Gene lists (TSV)
    │   ├── OXPHOS genes by complex
    │   ├── Mito-ribosome genes
    │   └── Pathway gene lists
    ├── hyp_genes/                             # Hypoxia markers
    │   └── human18_hn2.tsv
    └── figures/                               # Output figures
```

---

## Citation

If you use this code or data, please cite:

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

## Data Availability

- **Raw and processed data**: ArrayExpress accession E-MTAB-14981
- **Supplementary data**: [FigShare Project](https://figshare.com/projects/Shtolz_2025_Mitochondria_Hypoxia_Project/237023)
- **Code repository**: [GitHub](https://github.com/NoamSt/2025_Mitochondria_Hypoxia_Project)
- **Code release**: Zenodo DOI: 10.5281/zenodo.17845989

---

## Contact

**Corresponding Author:**
Dan Mishmar, PhD
Department of Life Sciences
Ben-Gurion University of the Negev
Email: dmishmar@bgu.ac.il
Tel: +972-8-6461355

**First Author:**
Noam Shtolz
Department of Life Sciences
Ben-Gurion University of the Negev

---

## Funding

This study was funded by:
- Israel Science Foundation (grants 372-17, 404-21)
- US Army Life Science Division (LS 80581-BB)
- Myles Thaler Genetics and Genomics Research Cathedra (DM)
- Negev Foundation scholarship for excellent students (NS)

---

## License

This article is licensed under a Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International License, which permits any non-commercial use, sharing, distribution and reproduction in any medium or format, as long as you give appropriate credit to the original author(s) and the source.

---

**Last Updated:** January 2026
**Version:** 1.0.0
