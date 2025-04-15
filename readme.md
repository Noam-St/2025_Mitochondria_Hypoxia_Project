# 2025 Mitochondria Hypoxia Project

## Overview
This project investigates the mitonuclear expression and hypoxia index in various datasets. It includes analysis of RNA sequencing data, hypoxia-related gene expression, and mitonuclear coordination. The project is structured to facilitate the exploration of mitochondrial and nuclear gene interactions under hypoxic conditions.

## Project Structure

- **Data**: Contains raw and processed data files, including RNA sequencing results, gene lists, and metadata.
  - `all_hct_and_ko.csv`, `all_mcf7_and_ko.csv`: Aggregated datasets for specific cell lines.
  - `HCT-116_hypx_2.tsv`, `HPMEC_hypx.tsv`: Hypoxia-related datasets.
  - `oxphos_genes.txt`, `hypoxia_associated_genes.xlsx`: Gene lists for analysis.

- **Notebooks**: Jupyter notebooks for data analysis and visualization.
  - `mitonuclear_expr_and_hypx_index.ipynb`: Main notebook for analyzing mitonuclear expression and hypoxia index.
  - `integ_rna_seq_analysis.ipynb`: Integrative RNA sequencing analysis.
  - `peppro_analysis.ipynb`: Analysis of promoter regions.

- **Scripts**: Python scripts for data preprocessing and analysis.
  - `annotation.py`, `correlation_analysis.py`: Scripts for annotation and correlation analysis.
  - `hypoxia_index.py`: Functions for calculating hypoxia indices.
  - `remove_batch_effects.py`: Script for batch effect removal.

- **Figures**: Contains generated plots and visualizations.
  - `oxphos_by_cell_line_assigned_all_coordinated_uncoordinated_log2foldchange.svg`: Visualization of OXPHOS coordination.
  - `ribosome_by_cell_line_unassigned_all_coordinated_uncoordinated_log2foldchange.svg`: Ribosome coordination plots.

## Key Features

- **Hypoxia Index Calculation**: Implements Buffa, Winters, and Lombardi hypoxia indices.
- **Mitonuclear Coordination Analysis**: Identifies coordinated and uncoordinated datasets for OXPHOS and ribosome genes.
- **Custom Visualization**: Generates boxplots and strip plots for gene expression data.

## Usage

1. **Setup**:
   - Ensure all dependencies are installed (e.g., `pandas`, `numpy`, `seaborn`, `matplotlib`).
   - Place raw data files in the `data` directory.

2. **Run Notebooks**:
   - Open and execute the Jupyter notebooks in the specified order for analysis.

3. **Scripts**:
   - Use the Python scripts for preprocessing and additional analysis.

## Dependencies

- Python 3.8+
- Required libraries: `pandas`, `numpy`, `seaborn`, `matplotlib`, `scipy`, `statannot`

## Contact
For questions or collaboration, please contact the project team.