import scipy.cluster.hierarchy
import scipy.spatial.distance
import numpy as np 
import scipy
import pandas as pd
import matplotlib.pyplot as plt 
import seaborn as sns
import networkx as nx 

def vlr_cluster(vlr_diff_df, p_value_df:pd.DataFrame=None):
    '''Cluster the VLR output values for clearer plotting.'''
    m = scipy.spatial.distance.squareform(vlr_diff_df.values)
    m = scipy.cluster.hierarchy.linkage(m, method='average')
    order = scipy.cluster.hierarchy.leaves_list(m)
    vlr_diff_df = vlr_diff_df.iloc[order, order]

    if p_value_df is not None:
        p_value_df = p_value_df.iloc[order, order]
        return vlr_diff_df, p_value_df
    return vlr_diff_df

def box_cox(values:np.ndarray, alpha:float=0.5):
    return (values ** alpha - 1) / alpha

def log_with_pseudocount(values:np.ndarray):
    return np.log(np.where(values == 0, 1, values))


def vlr_get_p_values(v, method:str='two-sided', n_samples=None):
    F = (n_samples - 2) * (1 - v) / v # Convert to a T statistic. 
    if method == 'two-sided':
        # return 2 * (1 - scipy.stats.t.cdf(np.abs(F), df=n_samples - 2))
        p_values = 2 * (scipy.stats.t.sf(np.abs(F), df=n_samples - 2))
    if method == 'less':
        p_values = scipy.stats.t.cdf(np.abs(F), df=n_samples - 2)
    if method == 'greater':
        # return 1 - scipy.stats.t.cdf(np.abs(F), df=n_samples - 2)
        p_values = scipy.stats.t.sf(np.abs(F), df=n_samples - 2)
    p_values = scipy.stats.false_discovery_control(p_values, method='bh', axis=None)
    return p_values.reshape(F.shape)


# Speed up computation for all genes at once. 
def vlr_get_diff(metat_df, transform=np.log, groups={'2025':['n_top_2025', 'n_middle_2025', 'n_bottom_2025'], '2024':['n_top_2024', 'n_middle_2024', 'n_bottom_2024']}):
    '''The variance in the log-ratios for a group with N samples between gene_i and gene_j is 
    var(log(x_i1/x_j1), ..., log(x_iN/x_jN)).'''
    # metat_df = metat_df.pivot(index=['sample_id', 'year', 'location'], values='read_count', columns='gene_id')
    metat_df = metat_df.pivot(index='gene_id', values='read_count', columns='sample_id') # Shape is (n_genes, n_samples)
    n_genes, n_samples = len(metat_df), len(metat_df.columns)

    groups = {group_id:np.where(np.isin(metat_df.columns.values, group))[0] for group_id, group in groups.items()}

    metat_arr = metat_df.values[None, :, :] # Add a dimension so shape is (1, n_genes, n_samples)
    metat_arr = np.broadcast_to(metat_arr, (n_genes, n_genes, n_samples))
    metat_arr = transform(metat_arr)

    f = lambda arr : np.var(arr - arr.transpose(1, 0, 2), ddof=0, axis=-1) * arr.shape[-1] # Make sure to weight according to the number of samples.
    L_total = f(metat_arr)
    np.fill_diagonal(L_total, 1) # Fill the diagonals to 1 to avoid nan error. Indicates perfect coordination. 
    
    L_group = {group_id:f(metat_arr[:, :, group_idxs]) for group_id, group_idxs in groups.items()}
    v = np.where(L_total != 1, np.sum(list(L_group.values()), axis=0) / L_total, 1)
    p_values = vlr_get_p_values(v, n_samples=n_samples)

    v_df = pd.DataFrame(v, index=metat_df.index, columns=metat_df.index)
    p_values_df = pd.DataFrame(p_values, index=metat_df.index, columns=metat_df.index)

    return v_df, p_values_df # Value closer to zero means that less of the total variance is explained by within-group variances. 


def vlr_plot_diff(vlr_diff_df, p_value_df:pd.DataFrame=None, max_p_value=0.2, title='', annotations:dict=None):
    '''Plot the VLR difference values as a heatmap.'''
    cmap = sns.color_palette('Blues', as_cmap=True)
    cmap.set_bad(color='lightgray')   # Color for NaNs

    figure_df = vlr_diff_df.copy()
    figure_df[figure_df > 1] = 1
    figure_df = 1 - figure_df # Convert from similarity metric (higher values are more similar) to distance metric.

    if p_value_df is not None:
        figure_df, p_value_df = vlr_cluster(figure_df, p_value_df=p_value_df.copy())
    else:
        figure_df = vlr_cluster(figure_df)

    if p_value_df is not None:
        print(f'vlr_plot_diff: {(p_value_df > max_p_value).values.sum()} values do not meet the p-value cutoff of {max_p_value}')
        # figure_df[p_value_df > max_p_value] = np.nan
        figure_df = figure_df.where(p_value_df <= max_p_value)

    if annotations is not None:
        figure_df.index = figure_df.index.map(annotations)
        figure_df.columns = figure_df.columns.map(annotations)

    fig, ax = plt.subplots(figsize=(0.25 * len(figure_df), 0.25 * len(figure_df)))

    sns.heatmap(figure_df, cbar=False, cmap=cmap, vmin=0, vmax=1, lw=0.7, linecolor='black') #, annot=figure_df, fmt='.1f')
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.set_xticks(np.arange(len(figure_df)) + 0.5) # , x_labels, fontsize='xx-small')
    ax.set_yticks(np.arange(len(figure_df)) + 0.5) #, y_labels, fontsize='xx-small')

    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(1)
        spine.set_edgecolor('black')
    ax.tick_params(axis='both', which='both', length=0)

    ax.minorticks_off()
    ax.set_title(title)
    plt.show()


# def _get_vlr(gene_id_i, gene_id_j, metat_df, transform=np.log, shuffle:bool=False):
#     '''The variance in the log-ratios for a group with N samples between gene_i and gene_j is 
#     var(log(x_i1/x_j1), ..., log(x_iN/x_jN)).'''
#     df = metat_df[metat_df.gene_id.isin([gene_id_i, gene_id_j])].copy()
#     df = df.pivot(index=['sample_id', 'year', 'location'], values='read_count', columns='gene_id')
#     df = df.reset_index()
#     n = len(df)

#     if shuffle:
#         df[gene_id_i] = np.random.permutation(df[gene_id_i].values)
#         df[gene_id_j] = np.random.permutation(df[gene_id_j].values)
        
#     L_total = np.var(transform(df[gene_id_i].values) - transform(df[gene_id_j].values), ddof=0)
#     if L_total == 0:
#         return np.nan 
    
#     L_group = dict()
#     # NOTE: Why use ddof=1?
#     for group, df_ in df.groupby('year'):
#         L_group[group] = len(df_) * np.var(transform(df_[gene_id_i].values) - transform(df_[gene_id_j].values), ddof=0)

#     v = np.sum(list(L_group.values())) / (n * L_total)
#     return v # Value closer to zero means that less of the total variance is explained by within-group variances. 

# def get_vlr_matrix(metat_df:pd.DataFrame, n_permutations:int=100, method:str='two-sided'):
#     gene_ids = metat_df.gene_id.unique()

#     matrix = np.empty((len(gene_ids), len(gene_ids)))
#     p_values = np.empty((len(gene_ids), len(gene_ids)))
#     # pbar = tqdm(desc='get_vlr_matrix', total=len(gene_ids)**2)
#     for i, gene_id_i in tqdm(list(enumerate(gene_ids)), desc='get_vlr_matrix'):
#         for j, gene_id_j in enumerate(gene_ids):
#             v, p_value = get_vlr(gene_id_i, gene_id_j, metat_df, n_permutations=n_permutations)
#             matrix[i, j] = v
#             p_values[i, j] = p_value
#     matrix_df = pd.DataFrame(matrix, index=pd.Series(gene_ids, name='gene_id'), columns=gene_ids)
#     p_value_df = pd.DataFrame(p_values, index=pd.Series(gene_ids, name='gene_id'), columns=gene_ids)
#     return matrix_df, p_value_df