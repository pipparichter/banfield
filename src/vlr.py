import scipy.cluster.hierarchy
import scipy.spatial.distance
import numpy as np 
import scipy
import pandas as pd
import matplotlib.pyplot as plt 
import seaborn as sns
import networkx as nx 
import sklearn.neighbors
import igraph
import sklearn
import matplotlib.cm
import matplotlib.colors


# def box_cox(values:np.ndarray, alpha:float=0.5):
#     return (values ** alpha - 1) / alpha

# def log_with_pseudocount(values:np.ndarray):
#     return np.log(np.where(values == 0, 1, values))

GROUPS = {'2025':['n_top_2025', 'n_middle_2025', 'n_bottom_2025'], '2024':['n_top_2024', 'n_middle_2024', 'n_bottom_2024']}
GROUP_ORDER = ['2025', '2024']


def vlr_plot_graph(graph):
    fig, ax = plt.subplots(figsize=(5, 5))

    pos = nx.spring_layout(graph, weight='weight')
    nx.draw_networkx_edges(graph, pos=pos, ax=ax, width=0.7)
    nx.draw_networkx_nodes(graph, pos=pos, node_size=30, node_color='lightgray', ax=ax)

    for spine in ax.spines.values():
        spine.set_visible(False)


# NOTE: A differential proportionality value close to 1 does not necessarily imply coordination, just that the variance between within treatment
# groups largely explains the total variance. 


# Not sure how valid this is given that variance can't be compared across gene pairs, but c'est la vie.
# Unclear if I should set a threshold on a gene-by-gene basis, there does seem to be decent variability.  

# Another way to build the graph would be to treat the ANOVA vectors as distances. 
def vlr_get_graph(df, radius:float=None, n_neighbors:int=10, plot:bool=False, metric:str='precomputed'):
    '''Build a graph using the VLR data, not the ANOVA-like data. ANOVA-lke data captures '''
    if (radius is not None):
        graph = sklearn.neighbors.radius_neighbors_graph(df.values, metric=metric, radius=radius, mode='distance')
    elif (n_neighbors is not None):
        graph = sklearn.neighbors.kneighbors_graph(df.values, metric=metric, mode='distance', n_neighbors=n_neighbors)
    else:
        raise Exception('vlr_get_graph: Either radius or n_neighbors must be specified.')

    graph = nx.from_numpy_array(graph.toarray()) #, edge_attr='distance') # Store the values in the distance attribute.
    for u, v, data in graph.edges(data=True):
        data['distance'] = data['weight']
        data['weight'] = 1.0 / (data['distance'] + 1e-8)
    graph = nx.relabel_nodes(graph, dict(enumerate(df.index.tolist())))

    if plot:
        vlr_plot_graph(graph)
    print(f'vlr_get_graph: Created a graph from the data with {len(list(graph.nodes()))} nodes and {len(list(graph.edges()))} edges.')
    return graph


# # https://medium.com/@swapnil.agashe456/leiden-clustering-for-community-detection-a-step-by-step-guide-with-python-implementation-c883933a1430
# TODO: Read about this algorithm and what the parameters actually do. 
def vlr_cluster_leiden(df, p_value_df:pd.DataFrame=None, radius:float=None, n_neighbors:int=10, resolution:float=1, **kwargs):
    # I think I want the graph to be scale-invariant, so probably best to use cosine. 
    graph = vlr_get_graph(df, radius=radius, n_neighbors=n_neighbors, plot=False, metric='cosine')
    graph = igraph.Graph.from_networkx(graph) # Convert to an igraph object.

    # communities = nx.community.leiden_communities(graph, weight='weight')
    communities = graph.community_leiden(weights='weight', resolution=resolution)
    print(f'vlr_cluster_leiden: Found {len(communities)} communities with Leiden clustering.')
    order = [idx for community in communities for idx in community]

    df = df.iloc[order, order]
    if p_value_df is not None:
        p_value_df = p_value_df.iloc[order, order]
        return df, p_value_df
    
    return df


# Hierarchical clustering doesn't seem to be working supremely well, not totally sure why. Basically seems like groups which should
# be connected are getting sorted into different clusters. Probably a network-based approach is better. 
def vlr_cluster_hierarchical(df, p_value_df:pd.DataFrame=None, **kwargs):
    '''Cluster the VLR output values for clearer plotting.'''
    df[df > 1] = 1
    m = scipy.spatial.distance.squareform((1 - df).values) # Make sure to convert to a distance metric. 
    m = scipy.cluster.hierarchy.linkage(m, method='average')
    order = scipy.cluster.hierarchy.leaves_list(m)
    df = df.iloc[order, order]

    if p_value_df is not None:
        p_value_df = p_value_df.iloc[order, order]
        return df, p_value_df
    return df

def _sort_clusters_kmeans(cluster_centers, labels):
    '''Sort the clusters according to which centers are the most similar for easier visualization.'''
    cluster_labels = np.unique(labels)
    cluster_linkage = scipy.cluster.hierarchy.linkage(scipy.spatial.distance.pdist(cluster_centers, metric='euclidean'), method='average')
    cluster_order = scipy.cluster.hierarchy.leaves_list(cluster_linkage)
    cluster_label_map = {cluster_labels[o]:i for i, o in enumerate(cluster_order)}
    labels = [cluster_label_map[label] for label in labels]
    return labels


def vlr_cluster_kmeans(df, p_value_df:pd.DataFrame=None, n_clusters:int=2, sort_clusters:bool=True, **kwargs):
    '''Cluster the VLR output values for clearer plotting.'''
    kmeans = sklearn.cluster.KMeans(n_clusters=n_clusters, random_state=42)
    kmeans.fit(df.values)
    labels = kmeans.labels_
    if sort_clusters:
        labels = _sort_clusters_kmeans(kmeans.cluster_centers, labels)

    order = np.argsort(labels)
    df = df.iloc[order, order]

    if p_value_df is not None:
        p_value_df = p_value_df.iloc[order, order]
        return df, p_value_df
    return df


# def vlr_variance_of_log_ratios(metat_df, transform=np.log):
#     # metat_df = metat_df.pivot(index=['sample_id', 'year', 'location'], values='read_count', columns='gene_id')
#     metat_df = metat_df.pivot(index='gene_id', values='read_count', columns='sample_id') # Shape is (n_genes, n_samples)
#     n_genes, n_samples = len(metat_df), len(metat_df.columns)

#     metat_arr = metat_df.values[None, :, :] # Add a dimension so shape is (1, n_genes, n_samples)
#     metat_arr = np.broadcast_to(metat_arr, (n_genes, n_genes, n_samples))
#     metat_arr = transform(metat_arr)

#     _variance_of_log_ratios = lambda arr : np.var(arr - arr.transpose(1, 0, 2), ddof=0, axis=-1) * arr.shape[-1] # Make sure to weight according to the number of samples.
#     L = _variance_of_log_ratios(metat_arr)
#     # np.fill_diagonal(L, 1) # Fill the diagonals to 1 to avoid nan error. Indicates perfect coordination. 
    
#     var_df = pd.DataFrame(L, index=metat_df.index, columns=metat_df.index)
#     return var_df


# I don't think a two-sided test makes sense here.
def vlr_anova_get_p_values(v, method:str='greater', n_samples=None):
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

# Disjointed versus emergent proportionality; lower values still indicate higher proportionality, but actual values will be lower for the 
# emergent proportionality. 

# What is the F statistic for the emergent version? Is it the same as before?

# Emergent proportionality is 1 - ((variance in the most variable group) / (total variance)), so the lower the number, the more
# variance is explained by one of the two groups, i.e. demonstrates emergent proportionality. 

# For disjoint proportionality, higher numbers indicate that intra-group variance explains most of the overall variance, indicating a lack 
# of disjoint proportionality; so, lower numbers are indicative of disjoint proportionality. 

# Don't want to use difference in the means directly, as this does not capture consistency across groups of samples, e.g. the mean could
# be negative if there is one outlier

def _mean_of_log_ratios(arr, transform=np.log):
    # L = transform(arr) - transform(arr).transpose(1, 0, 2) # Now the array stores the log ratios.
    L = transform(arr).transpose(1, 0, 2) - transform(arr) # Now the array stores the log ratios.
    E = L.mean(axis=-1, keepdims=True)
    return E.squeeze(axis=-1)

def _variance_of_log_ratios(arr, weighted:bool=False, transform=np.log):
    # Will have per-group weights and per-sample weights for each pair of genes. 
    w = arr * arr.transpose(1, 0, 2) if weighted else np.ones((arr.shape[0], arr.shape[0], arr.shape[-1]))
    L = transform(arr).transpose(1, 0, 2) - transform(arr) # Now the array stores the log ratios.
    omega = w.sum(axis=-1, keepdims=True) # Sum up total counts over all samples.  
    E = (w * L).sum(axis=-1, keepdims=True) / omega # Sum up over the samples. 
    var = (w * (L - E) ** 2).sum(axis=-1, keepdims=True) # / omega Don't divide by omega at the end, because it gets multiplied when computing v. 
    var = var.squeeze(-1)
    np.fill_diagonal(var, omega)
    return var

def vlr_anova(metat_df,  groups=GROUPS, method:str='emergent', weighted:bool=False, signed:bool=True, group_order=GROUP_ORDER):
    '''Computes the equivalent of an ANOVA test, comparing the variance between samples in the same "treatment group" to variance across all samples.'''

    # metat_df = metat_df.pivot(index=['sample_id', 'year', 'location'], values='read_count', columns='gene_id')
    metat_df = metat_df.pivot(index='gene_id', values='read_count', columns='sample_id') # Shape is (n_genes, n_samples)
    n_genes, n_samples = len(metat_df), len(metat_df.columns)

    groups = {group_id:np.where(np.isin(metat_df.columns.values, group))[0] for group_id, group in groups.items()}

    metat_arr = metat_df.values[None, :, :] # Add a dimension so shape is (1, n_genes, n_samples)
    metat_arr = np.broadcast_to(metat_arr, (n_genes, n_genes, n_samples))

    var_total = _variance_of_log_ratios(metat_arr, weighted=weighted)
    var_group = {group_id:_variance_of_log_ratios(metat_arr[:, :, group_idxs], weighted=weighted) for group_id, group_idxs in groups.items()}
    
    if (method == 'disjoint'):
        v = np.where(var_total != 1, np.sum(list(var_group.values()), axis=0) / var_total, 1)
    elif (method == 'emergent'):
        v = np.where(var_total != 1, 1 - np.max(list(var_group.values()), axis=0) / var_total, 1)
    else:
        raise Exception(f'vlr_anova: Specified method must be either emergent or disjoint.')
    
    p_values = vlr_anova_get_p_values(v, n_samples=n_samples)
    p_values_df = pd.DataFrame(p_values, index=metat_df.index, columns=metat_df.index)
    
    v = 1 - v 
    if signed:
        E_group = {group_id:_mean_of_log_ratios(metat_arr[:, :, group_idxs]) for group_id, group_idxs in groups.items()}
        E_diff = E_group[group_order[0]] - E_group[group_order[1]]
        v *= np.sign(E_diff)
    v_df = pd.DataFrame(v, index=metat_df.index, columns=metat_df.index)

    return v_df, p_values_df # Value closer to zero means that less of the total variance is explained by within-group variances. 


def vlr_plot_anova_heatmap(df, p_value_df:pd.DataFrame=None, max_p_value=0.2, title='', annotations:dict=None, scale:float=0.25, ticks:bool=True, cbar:bool=False, **kwargs):
    '''Plot the input values as a heatmap. For the ANOVA-like data, values closer to zero imply that much of the total variance
    is explained by the variance within treatment groups, so lower values imply larger "distance." The opposite is true for the 
    raw VLR values, so these should not be transformed.'''
    # cmap = sns.color_palette('Blues', as_cmap=True)
    # cmap.set_bad(color='black')   # Color for NaNs
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list('palette', ['gray', 'white', 'steelblue'])
    cmap.set_bad(color='white')   # Color for NaNs
    v_min, v_max = -1, 1

    figure_df = df.copy()
    # figure_df = figure_df.where(p_value_df <= max_p_value).fillna(1)

    figure_df, p_value_df = vlr_cluster_kmeans(figure_df, p_value_df=p_value_df, **kwargs)

    if p_value_df is not None:
        print(f'vlr_plot_diff: {(p_value_df > max_p_value).values.sum()} values do not meet the p-value cutoff of {max_p_value}')
        # figure_df[p_value_df > max_p_value] = np.nan
        figure_df = figure_df.where(p_value_df <= max_p_value)

    if annotations is not None:
        figure_df.index = figure_df.index.map(annotations)
        figure_df.columns = figure_df.columns.map(annotations)

    # figure_df = 1 - figure_df # Convert from similarity metric (higher values are more similar) to distance metric.
    fig, ax = plt.subplots(figsize=(scale * len(figure_df), scale * len(figure_df)))

    sns.heatmap(figure_df, cbar=cbar, cmap=cmap, lw=0.7, vmin=v_min, vmax=v_max, linecolor='white', cbar_kws={'shrink': 0.25}) #, annot=figure_df, fmt='.1f')
    ax.set_xlabel('')
    ax.set_ylabel('')
    if ticks:
        ax.set_xticks(np.arange(len(figure_df)) + 0.5) # , x_labels, fontsize='xx-small')
        ax.set_yticks(np.arange(len(figure_df)) + 0.5) #, y_labels, fontsize='xx-small')
    else:
        ax.set_xticks([])
        ax.set_yticks([])
    


    ax.minorticks_off()
    ax.set_title(title)
    plt.show()

    
    # cbar = ax.collections[0].colorbar
    # for spine in list(ax.spines.values()) + (list(cbar.ax.spines.values())):
    #     spine.set_visible(True)
    #     spine.set_linewidth(0.7)
    #     spine.set_edgecolor('black')
    # ax.tick_params(axis='both', which='both', length=0)

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
        
#     var_total = np.var(transform(df[gene_id_i].values) - transform(df[gene_id_j].values), ddof=0)
#     if var_total == 0:
#         return np.nan 
    
#     var_group = dict()
#     # NOTE: Why use ddof=1?
#     for group, df_ in df.groupby('year'):
#         var_group[group] = len(df_) * np.var(transform(df_[gene_id_i].values) - transform(df_[gene_id_j].values), ddof=0)

#     v = np.sum(list(var_group.values())) / (n * var_total)
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