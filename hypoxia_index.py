"""
This script calculates the hypoxia index for a given dataset. The index is calculated as the weighted sum of all known hypoxia related genes. 
The weights are defined as the position of the genes in the hypoxia related gene list, taken from: "Meta-Analysis of Hypoxic Transcriptomes from Public Databases"
"""

import pandas as pd
import os

PATH = os.getcwd()
HYPX_GENE_PATH = os.path.join(PATH, 'data', 'hyp_genes', 'human18_hn2.tsv')
HYPX_GENE_LIST = list(pd.read_csv(os.path.join(PATH, 'data', 'hyp_genes', 'human18_hn2.tsv'),sep = '\t').sort_values(by = 'HN-score', ascending = False)['Gene'])

def generate_hypx_weights(hypx_sorted_list = HYPX_GENE_LIST, topN = 100):
    """
    Recieve a sorted list of genes by their upregulation level in hypoxia.
    Return a dictionary of the topN genes and their weights.
    """

    hypx_weights = {}
    for i in range(topN):
        # The weight of the gene is the position of the gene in the list, divided by the total number of genes in the list.
        hypx_weights[hypx_sorted_list[i]] = (topN - i)/topN
    return hypx_weights

def hypx_index(df, hypx_weights, gene_col = 'gene', log2fc_col = 'log2FoldChange'):
    """
    Recieve a dataframe that represents the log2foldChange output from DESeq2  and a dictionary of genes and their weights.
    Return the hypoxia index of the dataframe.
    """
    hypx_index = 0
    c = 0
    for gene in hypx_weights.keys():
        if gene in df[gene_col].unique():
            c += 1
            hypx_index += df.loc[df[gene_col] == gene, log2fc_col].values * hypx_weights[gene]
    hypx_index = float(hypx_index/c) 
    return hypx_index