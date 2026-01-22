"""
Given a dataframe where each row is a position, with data on forward and reverse RPM coverage, calculate the pausing index for each position based on:
https://www.biorxiv.org/content/10.1101/2023.10.17.562821v1.full

The steps for calculating the pausing index are as follows:
1. Per-strand, convert the RPM data to a Z-score relative to the mean and standard deviation of read-end coverage at a 201 base window centered on the position.
2. Base positions with a Z-score >= 3 and at least 1 read-end coverage on the same strand are considered peaks.
3. Peaks that appear across samples are considered pausing sites.
"""
import pandas as pd
import numpy as np

def mtdna_region(pos, window, total, left, inclusive = True) -> str:
    """
    Designed for circular DNA, return a list of positions (INDEX 1 BASED) to the left or right side of pos

    Parameters
    ----------
    pos : int
        The current position to return a window around
    window : int
        The size of the window
    total : int
        The total size of the DNA
    left : bool
        True if the window should be to the left of the
    inclusive :
        (Default value = True)

    Returns
    -------
    positions : list
        A list of positions (INDEX 1 BASED) to the left or right side of pos
    """

    if type(left) != bool: raise TypeError(f'Parameter left must be either True (left side window) or False (right side window! Given {left} instead')
    if total < window or total < pos:
        raise ValueError(f'The total size of the DNA must be smaller than both the window and the position!\nParameters given:\nTotal = {total}\nPosition = {pos}\nWindow = {window}')
    
    positions = []
    if left:
        if pos - window < 1:
            positions += list(range(total - abs(window - pos) + (1 if inclusive else 0), total + 1))
            positions += list(range(1, pos + (1 if inclusive else 0)))
        else:
            positions += list(range(pos - window + (1 if inclusive else 0), pos + (1 if inclusive else 0)))
    else: #right
        if pos + window > total:
            positions += list(range(1, (pos + window) - total + (0 if inclusive else 1)))
            positions += list(range(pos + (0 if inclusive else 1), total + 1))
        else:
            positions += list(range(pos + (0 if inclusive else 1), pos + window + 1))
    return positions

def z_score(n, mean, std) -> float:
    """
    Compute Z-score.

    Parameters
    ----------
    n : int 
        The number to compute the Z-score for
    mean : float
        The mean of the distribution
    std : float
        The standard deviation of the distribution
    
    Returns
    -------
    float
        The Z-score
    
    References
    ----------
    https://en.wikipedia.org/wiki/Standard_score
    """
    try: return (n - mean) / std
    except ZeroDivisionError: return 0

def iter_over_windows(coverage, window_size=201, cov_thresh = 0.1) -> list:
    """
    Iterate over a list of values and compute Z-scores for each window.

    Parameters
    ----------
    coverage : list
        The list of values to compute Z-scores for
    window_size : int
        The size of the window to compute the Z-scores for
    
    Returns
    -------
    list
        The list of Z-scores
    """
    z_scores = []
    for i in range(len(coverage)):
        left_positions = mtdna_region(pos = i, window = window_size//2, total = 16569, left = True, inclusive = True)
        right_positions = mtdna_region(pos = i, window = window_size//2, total = 16569, left = False, inclusive = True)
        window_positions = left_positions + right_positions
        window_values = list(np.take(coverage, [i-1 for i in window_positions]))
        mean_cov = np.mean(window_values)
        std_cov = np.std(window_values)
        z_scores.append(z_score(coverage[i], mean_cov ,std_cov))
    return z_scores

def pausing_index_calculator(df: pd.DataFrame, window_size: int=201, sample_col_name = 'sample', overlap_range = 100, y = 'read_end_RPM', cov_y = 'RPM',  z_score_thresh = 3, RPM_thresh = 1) -> pd.DataFrame:
    """
    Calculate the pausing index for each position in a given dataframe.

    Parameters
    ----------
    df : pd.DataFrame
        The dataframe to calculate the pausing index for
    threshold : int
        The threshold to consider a peak
    
    Returns
    -------
    pd.DataFrame
        The dataframe with the pausing index calculated
    """
    # If y is a list of columns, calculate window_Z_score per column
    if type(y) != list:
        y_list = [y]
    else:
        y_list = y
        # Create a new dataframe to store the pausing index
    pausing_index_df = pd.DataFrame(columns = df.columns)
    n_samples = len(df[sample_col_name].unique())
    # Iterate through each sample
    for sample in df[sample_col_name].unique():
        # Get the data for the current sample
        sample_df = df[df[sample_col_name] == sample]
        for y in y_list:
            # Get the RPM data for the current sample
            rpm_data_pos = sample_df[f'pos_{y}'].values
            rpm_data_neg = sample_df[f'neg_{y}'].values
            rpm_data = sample_df[y].values

            # Calculate the Z-scores for the RPM data
            z_scores_pos = iter_over_windows(rpm_data_pos, window_size=window_size)
            z_scores_neg = iter_over_windows(rpm_data_neg, window_size=window_size)
            z_scores_both = iter_over_windows(rpm_data, window_size=window_size)


            # Add the Z-scores to the dataframe
            sample_df.loc[:, f'pos_{y}_window_Z'] = z_scores_pos
            sample_df.loc[:, f'neg_{y}_window_Z'] = z_scores_neg
            sample_df.loc[:, f'{y}_window_Z'] = z_scores_both

            # Check if the Position qualifies as a peak
            sample_df.loc[:, f'pos_{y}_peak'] = (sample_df[f'pos_{y}_window_Z'] >= z_score_thresh) & (sample_df[f'pos_{cov_y}'] > RPM_thresh)
            sample_df.loc[:, f'neg_{y}_peak'] = (sample_df[f'neg_{y}_window_Z'] >= z_score_thresh) & (sample_df[f'neg_{cov_y}'] > RPM_thresh)
            sample_df.loc[:, f'{y}_peak'] = (sample_df[f'{y}_window_Z'] >= z_score_thresh) & (sample_df[f'{cov_y}'] > RPM_thresh)
        pausing_index_df = pd.concat([pausing_index_df, sample_df], ignore_index=True)

    return pausing_index_df

   
