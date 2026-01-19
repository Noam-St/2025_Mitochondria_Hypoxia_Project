# 2025 Mitochondria Hypoxia Project - GitHub project for the paper [[link_here](https://www.nature.com/articles/s42003-025-09457-y)]
## Overview
In this project, Shtolz et al. investigated the effect of hypoxia on mito-nuclear expression coordination in a large compendium of publicly available RNA-seq data. It includes analysis of RNA sequencing data, hypoxia-related gene expression, and mitonuclear coordination. The project is structured to facilitate the exploration of mitochondrial and nuclear gene interactions under hypoxic conditions.

## Project Structure

- **Notebooks**: Jupyter notebooks for data analysis and visualization.
  - `mitonuclear_expr_and_hypx_index.ipynb`: Main notebook for analyzing mitonuclear expression and hypoxia index.
  - `integ_rna_seq_analysis.ipynb`: Integrative RNA sequencing analysis.
  - `peppro_analysis.ipynb`: Analysis of promoter regions.

- **Scripts**: Python scripts for data preprocessing and analysis.
  - `annotation.py`, `correlation_analysis.py`: Scripts for annotation and correlation analysis.
  - `hypoxia_index.py`: Functions for calculating hypoxia indices.
  - `remove_batch_effects.py`: Script for batch effect removal.

## Key Features

- **Hypoxia Index Calculation**: Implements Buffa, Winters, and Lombardi hypoxia indices.
- **Mitonuclear Coordination Analysis**: Identifies coordinated and uncoordinated datasets for OXPHOS and ribosome genes.
- **Custom Visualization**: Generates boxplots and strip plots for gene expression data.

## Usage

1. **Setup**:
   - Ensure all dependencies are installed (e.g., `pandas`, `numpy`, `seaborn`, `matplotlib`).
   - Place raw data files downloaded from [figshare](https://figshare.com/account/home#/projects/237023).

2. **Run Notebooks**:
   - Open and execute the Jupyter notebooks in the specified order for analysis.

3. **Scripts**:
   - Use the Python scripts for preprocessing and additional analysis.

## Dependencies

- Python 3.8+
- Required libraries: `pandas`, `numpy`, `seaborn`, `matplotlib`, `scipy`, `statannot`

## Contact
For questions or collaboration, please contact the project team.
