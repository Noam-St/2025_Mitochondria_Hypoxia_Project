"""
This code removes batch effects from a given matrix and sampleData using ComBat. (Designed for expression matrices)
"""

# Import packages
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
from combat.pycombat import pycombat

def retrieve_batches(expression_matrix, count_data, batch_cols, sample_col):
    """
    Receive a sampleData dataframe and a batch column name and return a codified list of batches.
    """
    sampleData = count_data.groupby(sample_col).agg({i : 'first' for i in batch_cols}).reset_index(drop=False)
    sampleData['batch'] = sampleData[batch_cols].apply(lambda x: '_'.join(x.astype(str)), axis=1)
    sampleData['batch_digit'] = sampleData['batch'].astype('category').cat.codes
    batch_list = sampleData['batch_digit'].tolist()
    return batch_list, sampleData

def remove_batch_effects(expression_matrix, sampleData, batch_list):
    """
    Receive an expression matrix and a sampleData dataframe and return a corrected expression matrix.
    """

    # Remove batch effects
    corrected_expression_matrix = pycombat(expression_matrix, batch_list)
    return corrected_expression_matrix

def qc_plots(expression_matrix, expression_matrix_corrected, sampleData, id_col, sample_col, value_col, batch_cols):
    """
    Receive an expression matrix, a corrected expression matrix, a sampleData dataframe, and the names of the columns that contain the batch information and return a pre- and post- batch correction boxplot of the distribution of counts per gene to visualize the batch effect correction
    """
    expression_matrix_corrected_long = expression_matrix_corrected.reset_index(drop=False).melt(id_vars = id_col, var_name = sample_col, value_name = value_col)
    expression_matrix_long = expression_matrix.reset_index(drop=False).melt(id_vars = id_col, var_name = sample_col, value_name = value_col)
    # Add dataset and cell_type columns
    expression_matrix_corrected_long = expression_matrix_corrected_long.merge(sampleData[batch_cols + [sample_col, 'batch']], on = sample_col, how = 'left')
    expression_matrix_long = expression_matrix_long.merge(sampleData[batch_cols + [sample_col, 'batch']], on = sample_col, how = 'left')

    # Create a pre-combat boxplot of the distribution of counts per gene
    _, axes = plt.subplots(2, 1, figsize=(7, 3), sharex = True, sharey = True)
    sns.boxplot(ax = axes[0], data = expression_matrix_long, x = sample_col, y = value_col, hue = 'batch', linewidth = 0.5, dodge = False, showfliers = False)
    sns.boxplot(ax = axes[1], data = expression_matrix_corrected_long, x = sample_col, y = value_col, hue = 'batch', linewidth = 0.5, dodge = False, showfliers = False)
    # Remove legend
    axes[0].get_legend().remove()
    axes[1].get_legend().remove()
    axes[0].set_title('Before Batch Correction', fontsize = 10)
    axes[1].set_title('After Batch Correction', fontsize = 10)
    axes[0].set_xlabel('')
    axes[1].set_xlabel('')
    axes[0].set_ylabel('')
    axes[1].set_ylabel('')
    # Plot a line at mean y value
    axes[0].axhline(y = expression_matrix_long[value_col].median(), color = 'black', linestyle = '--', linewidth = 0.5)
    axes[1].axhline(y = expression_matrix_corrected_long[value_col].median(), color = 'black', linestyle = '--', linewidth = 0.5)
    # Despine
    sns.despine()
    # Remove x ticklabels and tick marks
    axes[0].set_xticklabels([])
    axes[1].set_xticklabels([])
    axes[0].tick_params(axis = 'x', length = 0)
    axes[1].tick_params(axis = 'x', length = 0)

def main(expression_matrix, count_data, id_col = 'gene', value_col = 'counts', sample_col = 'run_accession', batch_cols = ['dataset', 'cell_type'], qc = True):
    """
    Receive an expression matrix and a sampleData dataframe and return a corrected expression matrix.
    """
    # Check if the samples are the rows or columns of expression_matrix
    if expression_matrix.shape[0] == len(count_data[sample_col].unique()):
        expression_matrix = expression_matrix.T
    # Drop rows with all zeros
    expression_matrix = expression_matrix.loc[~(expression_matrix == 0).all(axis=1)]
    # Retrieve batches
    batch_list, sampleData = retrieve_batches(expression_matrix, count_data, batch_cols, sample_col)

    # Remove batch effects
    expression_matrix_corrected = remove_batch_effects(expression_matrix, sampleData, batch_list)
    # QC plots
    if qc:
        qc_plots(expression_matrix, expression_matrix_corrected, sampleData, id_col, sample_col, value_col, batch_cols)
    return expression_matrix_corrected, sampleData

def main_split_by(expression_matrix, count_data, id_col = 'gene', value_col = 'counts', sample_col = 'run_accession', batch_cols = ['dataset'], split_by = 'cell_type'):
    """
    Receive an expression matrix, split it by a column, do a batch correction for batch_cols in each split and return a concatenated corrected expression matrix.
    """
    expr_matrix_list = []
    sampleData_list = []
    for split in count_data[split_by].unique():
        cur_expr_matrix = expression_matrix.loc[:, count_data[count_data[split_by] == split]['sample'].unique()]
        cur_count_data = count_data[count_data[split_by] == split]
        # Dont do batch correction if there is only one batch
        if len(count_data[count_data[split_by] == split][batch_cols].drop_duplicates()) == 1:
            print('Only one batch in split ' + split + '. No batch correction done.')
            _, cur_sampleData = retrieve_batches(cur_expr_matrix, cur_count_data, batch_cols, sample_col)
        else:
            cur_expr_matrix, cur_sampleData = main(
                cur_expr_matrix, cur_count_data, id_col = id_col, value_col = value_col, sample_col = sample_col, batch_cols = batch_cols, qc = False)
        if type(cur_expr_matrix) != pd.DataFrame:
            print(f'Empty expression matrix for split {split}. Skipping...')
            continue
        expr_matrix_list.append(cur_expr_matrix)
        sampleData_list.append(cur_sampleData)
    # Concatenate expression matrices and sampleData
    expression_matrix_corrected = pd.concat(expr_matrix_list, axis = 1)
    sampleData = pd.concat(sampleData_list, axis = 0)
    # QC plots
    qc_plots(expression_matrix, expression_matrix_corrected, sampleData, id_col, sample_col, value_col, batch_cols)
    return expression_matrix_corrected, sampleData

    