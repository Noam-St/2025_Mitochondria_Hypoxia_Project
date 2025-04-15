import pandas as pd 
import os
from diffexpr.py_deseq import py_DESeq2

def run_deseq(
    counts : pd.DataFrame,
    metadata : pd.DataFrame,
    out_name : str,
    counts_format = 'long',
    design_formula = '~Temperature',
    contrast = None,
    id_col = 'Name',
    sample_col = 'Sample',
    value_col = 'NumReads',
    ):
    """
    DESeq2 wrapper for DESeq2.py
    """
    if counts_format == 'long': counts = counts.pivot_table(index = id_col, columns = sample_col, values = value_col).reset_index()
    metadata.sample_name = metadata.sample_name.str.replace('Sample_', '')
    metadata = metadata.sort_values('sample_name')
    metadata.index = metadata.sample_name
    counts.loc[:, 'A':'I'] = counts.loc[:, 'A':'I'].astype(int)
    dds = py_DESeq2(count_matrix = counts, design_matrix = metadata, design_formula = design_formula, gene_column= id_col)

    dds.run_deseq()
    if contrast:
        dds.get_deseq_result(contrasts = contrast)
        res = dds.deseq_result
        res.to_csv(os.path.join('data', 'counts', os.path.split(out_name)[-1].replace('.csv','') + contrast[1] + '_vs_' + contrast[2] +  '.csv'))
    norm = dds.normalized_count() #DESeq2 normalized count
    norm.to_csv(os.path.join(out_name), index=True) #Save the normalized count:
    return norm