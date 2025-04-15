import pandas as pd
import numpy as np
import statsmodels.api as sm
import random
from patsy import dmatrices
from importlib import reload
import consts
reload(consts)

def extract_resids(df, var_col, predic_cols):
    """
    Recieve a dataframe where each row is an observation of var_col and predic_cols.
    Perform a regression of var_col on predic_cols and return the residuals to regress out the effect of predic_cols on var_col

    Parameters
    ----------
    df : pd.DataFrame
        A dataframe where each row is an observation of var_col and predic_cols
    var_col : str
        The name of the column containing the variable to regress on
    predic_cols : list
        A list of the names of the columns containing the predictors
    
    Returns
    -------
    pd.Series
        A series containing the residuals from the regression
    """
    mod_format = var_col + ' ~ ' + ' + '.join(predic_cols)
    y, X = dmatrices(mod_format, data=df, return_type='dataframe')
    model = sm.OLS(y, X)
    results = model.fit()
    return results.resid

def resid_per_grp(df, grp_cols, var_col, predic_cols):
    """
    Recieve a dataframe where each row is an observation of var_col and predic_cols. Regress var_col on predic_cols per group(s) and return the residuals to regress out the effect of predic_cols on var_col

    Parameters
    ----------
    df : pd.DataFrame
        A dataframe where each row is an observation of var_col and predic_cols
    grp_cols : str or list
        The name(s) of the column(s) to group by
    var_col : str
        The name of the column containing the variable to regress on
    predic_cols : list
        A list of the names of the columns containing the predictors
    
    Returns
    -------
    pd.DataFrame
        A dataframe containing the residuals from the regression per group(s)
    """
    if type(grp_cols) == str:
        grp_cols = [grp_cols]
    if type(predic_cols) == str:
        predic_cols = [predic_cols]
    
    df = df.copy()
    small_df = df[grp_cols + [var_col] + predic_cols]
    resid_df = small_df.groupby(grp_cols, sort = False).apply(lambda x: extract_resids(x, var_col, predic_cols)).reset_index().set_index(f'level_{len(grp_cols)}')
    resid_df.columns = grp_cols + ['resid']
    df = df.merge(resid_df[['resid']], left_index = True, right_index = True, how = 'left', )
    return df

def downsample(df, samp_col, grp_cols, n = 'min'):
    """
    Randomly downsample a dataframe to n observations per group(s). n being the number of observations in the smallest group by default (stratified sampling)

    Parameters
    ----------
    df : pd.DataFrame
        A dataframe to downsample
    samp_col : str
        The name of the column containing the sample names
    grp_cols : str or list
        The name(s) of the column(s) to group by
    n : int or str
        The number of observations to downsample to per group(s). 'min' by default, which is the number of observations in the smallest group in df
    
    Returns
    -------
    pd.DataFrame
        A dataframe containing the downsampled observations
    """
    df = df.copy()
    if type(grp_cols) == str:
        grp_cols = [grp_cols]
    samples = df.groupby(samp_col).agg({i : 'first' for i in grp_cols})
    min_n = samples.value_counts(grp_cols).min()
    if n == 'min':
        n = min_n
        replace = False
    else:
        replace = True if n > min_n else False
    samples = samples.groupby(grp_cols, group_keys = False).apply(lambda x: x.sample(n = n, replace = replace)).reset_index(samp_col)[samp_col].to_list()
    return df.loc[df[samp_col].isin(samples), :]

def correlation_analysis(df, grp_cols, var_col = 'resid', gene_col = 'Name', samp_col = None, n = 'min', corr_type = 'spearman', iterations = 1, agg_type = 'median', genome_col = 'Genome', nuc_name = 'NUC', mt_name = 'MT', sort = True):
    """
    Returns a dataframe with the correlation between the residuals from a regression of var_col on predic_cols and the residuals from a regression of var_col on predic_cols, per group(s)

    If samp_col is provided, the dataframe is downsampled to n observations per group(s) before the analysis. n being the number of observations in the smallest group by default

    Parameters
    ----------
    df : pd.DataFrame
        A dataframe where each row is an observation of var_col and predic_cols
    grp_cols : str or list
        The name(s) of the column(s) to group by
    var_col : str
        The name of the column containing the continuous variable
    gene_col : str
        The name of the column containing the gene names
    samp_col : str
        The name of the column containing the sample names
    n : int or str
        The number of observations to downsample to per group(s). 'min' by default, which is the number of observations in the smallest group in df
    corr_type : str
        The type of correlation to calculate. 'spearman' by default
    iterations : int
        The number of times to calculate the correlation. 1 by default
    agg_type : str
        The type of aggregation to perform on the correlations. 'median' by default
    genome_col : str
        The name of the column containing the genome names
    nuc_name : str
        The name of the nuclear genome
    mt_name : str
        The name of the mitochondrial genome
    sort : bool
        Whether to sort the dataframe by the correlation values. True by default

    Returns
    -------
    pd.DataFrame
        A dataframe containing the correlations between the residuals from a regression of var_col on predic_cols and the residuals from a regression of var_col on predic_cols, per group(s)
    """
    df = df.copy()
    dfs = []
    for _ in range(iterations):
        if samp_col:
            df = downsample(df, samp_col, grp_cols, n = n)
        corr_df = df.pivot_table(index = samp_col, columns = gene_col, values = var_col).dropna(axis = 1).corr(method = corr_type)
        corr_df = corr_df.stack()
        corr_df.index = corr_df.index.rename([gene_col + '_1', gene_col + '_2'])
        corr_df = corr_df.reset_index().rename(columns = {0 : 'corr'})
    dfs.append(corr_df)
    corr_df = pd.concat(dfs, axis = 0).groupby([gene_col + '_1', gene_col + '_2']).agg(agg_type).reset_index()
    corr_df = corr_df.pivot_table(index = gene_col + '_1', columns = gene_col + '_2', values = 'corr')
    if sort:
        all_genes = corr_df.columns.to_list()
        nuc_genes = sorted([i for i in df[df[genome_col] == nuc_name][gene_col].unique() if i in all_genes], key = lambda x: corr_df.loc[x,:].mean(), reverse = True)
        mt_genes = consts.MT_GENES
        corr_df = corr_df.loc[mt_genes + nuc_genes, mt_genes + nuc_genes]
    return corr_df

def randomly_pick_n_genes(df, n = None, appended_genes = consts.MT_GENES, cond_col = 'Complex', filter_value = None, gene_col = 'gene', genome_col = 'Genome', nuc_name = 'NUC'):
    """
    Recieve a df containing genes and expression in long format. Return a df containing n random genes from the nuclear genome and the appended genes

    Parameters
    ----------
    df : pd.DataFrame
        A dataframe containing genes and expression in long format
    n : int
        The number of genes to pick. If None, n will be equal to the number of oxphos-related nuclear genes
    appended_genes : list
        A list of genes to append to the randomly picked genes (always included)
    cond_col : str
        The name of the column containing the condition
    filter_value : str
        The value to filter the condition column by
    gene_col : str
        The name of the column containing the gene names
    genome_col : str
        The name of the column containing the genome names
    nuc_name : str
        The name of the nuclear genome

    Returns
    -------
    pd.DataFrame
        A dataframe containing n random genes from the nuclear genome and the appended genes
    """
    if filter_value is None:
        genes = df[df[cond_col].isna()][gene_col].unique()
    else:
        genes = df[df[cond_col] == filter_value][gene_col].unique()
    if n is None and filter_value is None:
        n = len(df[(df[cond_col].notna()) & (df[genome_col] == nuc_name)][gene_col].unique())
    elif n is None and filter_value is not None:
        n = len(df[(df[cond_col] == filter_value) & (df[genome_col] == nuc_name)][gene_col].unique())
    rnd_genes = list(np.random.choice(genes, n, replace=False)) + appended_genes
    rnd_genes_df = df[df[gene_col].isin(rnd_genes)]
    return rnd_genes_df

def randomize_and_correlate(df,
                            cond_col = 'Complex', filter_value = None, gene_col = 'gene', genome_col = 'Genome', nuc_name = 'NUC', appended_genes = consts.MT_GENES, n_random = None,
                            grp_cols = ['gene', 'treatment'], var_col = 'counts', predic_cols = ['cell_type'], use_resids = True, corr_type = 'spearman', iterations = 1, agg_type = 'median', sort = True, treatment_col = 'treatment', treatment_name = 'Hypoxia', control_name = 'Normoxia', samp_col = 'sample'):
    """
    Randomly pick n genes from the nuclear genome and the appended genes, and correlate the residuals from a regression of var_col on predic_cols, per group(s)
    """
    rnd_genes_df = randomly_pick_n_genes(df, n = n_random, appended_genes = appended_genes, cond_col = cond_col, filter_value = filter_value, gene_col = gene_col, genome_col = genome_col, nuc_name = nuc_name)
    if use_resids:
        rnd_genes_df = resid_per_grp(rnd_genes_df, grp_cols = grp_cols, var_col = var_col, predic_cols = predic_cols)
        var_col = 'resid'
    treat_df = rnd_genes_df[rnd_genes_df[treatment_col] == treatment_name]
    ctrl_df = rnd_genes_df[rnd_genes_df[treatment_col] == control_name]
    treat_corr = correlation_analysis(treat_df, grp_cols = grp_cols, var_col = var_col, gene_col = gene_col, corr_type = corr_type, iterations = iterations, agg_type = agg_type, sort = sort, samp_col = samp_col)
    ctrl_corr = correlation_analysis(ctrl_df, grp_cols = grp_cols, var_col = var_col, gene_col = gene_col, corr_type = corr_type, iterations = iterations, agg_type = agg_type, sort = sort, samp_col = samp_col)
    
    return treat_corr, ctrl_corr
    