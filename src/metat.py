import os
import skbio
import pandas as pd
import glob 
import numpy as np 
from src.bbduk import bbduk_load
from scipy.stats import gmean
import re


def metat_get_diff(metat_df, genome_id=None, normalization='clr', ref_gene_ids=[], verbose:bool=False, threshold:int=10):

    assert genome_id is not None, 'metat_get_diff: Genome ID must be specified.'
    
    # mask = ((metat_df.genome_id == genome_id) | metat_df.gene_id.isin(ref_gene_ids))
    mask = (metat_df.genome_id == genome_id)
    metat_df = metat_df[mask].copy()
    sample_ids = metat_df.sample_id.unique()
    assert len(sample_ids) == 2, f'metat_get_diff: There should be two samples represented in the input DataFrame, but saw {' '.join(sample_ids)}.'

    n = metat_df.gene_id.nunique()
    metat_df = metat_filter(metat_df, threshold=1, min_samples=1)
    print(f'metat_get_diff: Removed {n - metat_df.gene_id.nunique()} genes with zero expression in both samples.')

    if (normalization == 'alr') and verbose:
        print(f'get_diff: {np.isin(ref_gene_ids, metat_df.gene_id.unique()).sum()} out of {len(ref_gene_ids)} reference genes retained after filtering.')
    metat_df = metat_df.drop_duplicates(['gene_id', 'year'])
    metat_df = metat_normalize(metat_df, method=normalization, ref_gene_ids={genome_id:ref_gene_ids}, add_pseudocount='mzr')
    diff_df = dict()
    diff_df = {year:df.set_index('gene_id').sort_index() for year, df in metat_df.groupby('year')}

    diff_df = diff_df['2025'][['read_count_normalized']] - diff_df['2024'][['read_count_normalized']]
    diff_df = diff_df.reset_index().rename(columns={'read_count_normalized':'diff'})
    diff_df['genome_id'] = genome_id
    return diff_df


def metat_add_library_size(metat_df, bbduk_data_dir='../data/bbduk'):
    bbduk_df = bbduk_load(bbduk_data_dir)
    bbduk_df = bbduk_df[bbduk_df.index.str.contains('metat')].copy()
    metat_df['library_size'] = metat_df.sample_id.map(bbduk_df.library_size)
    return metat_df


# Zero handling 
# (1) Add a pseudocount of one. 
# (2) Multiplicative zero replacement, which preserves the ratios betweek gene expression levels. 
#   For differential expression within an organism, I think this should be done on a per-organism basis
#   (and also per-sample)


def _metat_add_pseudocounts_multiplicative(metat_df:pd.DataFrame):
    n = metat_df.read_count.sum()
    metat_df['read_count'] = skbio.stats.composition.multi_replace(metat_df.read_count.values) * n
    return metat_df

def _metat_add_pseudocounts_impute_median(metat_df:pd.DataFrame):
    median = metat_df.read_count.median()
    metat_df['read_count'] = np.where(metat_df.read_count == 0, median, metat_df.read_count)
    return metat_df  

def _metat_add_pseudocounts_ones(metat_df:pd.DataFrame):
    metat_df['read_count'] = metat_df.read_count.values + 1
    return metat_df

def metat_add_pseudocounts(metat_df:pd.DataFrame, method:str='ones'):
    methods = dict()
    methods['ones'] = _metat_add_pseudocounts_ones
    methods['mzr'] = _metat_add_pseudocounts_multiplicative
    methods['median'] = _metat_add_pseudocounts_impute_median

    # assert 'read_count_original' not in metat_df.columns, 'metat_add_pseudocounts: Seems like the pseudocounts have already been added!'
    if 'read_count_original' in metat_df.columns:
        # print('metat_add_pseudocounts: Skipping pseudocounts, seems like they are already added.')
        return metat_df
    
    metat_df['read_count_original'] = metat_df.read_count.values.copy()
    metat_df_ = list()
    for (sample_id, genome_id), df in metat_df.groupby(['sample_id', 'genome_id'], group_keys=False):
        df = methods[method](df) # MZR should preserve the ratios for a specific organism in a specific sample.
        metat_df_.append(df)
    return pd.concat(metat_df_)

metat_check_single_genome = lambda metat_df : metat_df.genome_id.nunique() == 1
metat_check_single_sample = lambda metat_df : metat_df.sample_id.nunique() == 1

# Reference provided as {genome_id:[(genome_id, gene_id), (genome_id, gene_id)]}

def _metat_normalize_alr(metat_df:pd.DataFrame, ref_gene_ids:dict=None):
    assert metat_check_single_sample(metat_df), '_metat_get_normalization_factors_alr'
    normalization_factors = dict()
    for genome_id, genes in ref_gene_ids.items():
        
        read_counts = metat_df[metat_df.gene_id.isin(genes)]
        read_counts = read_counts.read_count
        normalization_factors[genome_id] = gmean(read_counts)
    metat_df = metat_df.reset_index() # Restore the original index. 
    metat_df['normalization_factor'] = metat_df.genome_id.map(normalization_factors)
    metat_df['read_count_normalized'] = np.log(metat_df.read_count) - np.log(metat_df.normalization_factor)
    return metat_df

def _metat_normalize_clr(metat_df, ref_gene_ids:dict=dict()):
    normalization_factor = metat_df.groupby('genome_id').apply(lambda df : gmean(df.read_count.values), include_groups=False).to_dict()
    metat_df['normalization_factor'] = metat_df.genome_id.map(normalization_factor)
    metat_df['read_count_normalized'] = np.log(metat_df.read_count) - np.log(metat_df.normalization_factor)
    return metat_df


def metat_normalize(metat_df:pd.DataFrame, ref_gene_ids:dict=dict(), add_pseudocount:str='mzr', method:str='clr'):
    methods = dict()
    methods['clr'] = _metat_normalize_clr
    methods['alr'] = _metat_normalize_alr

    if 'read_count_original' not in metat_df.columns:
        metat_df = metat_add_pseudocounts(metat_df.copy(), method=add_pseudocount)

    metat_df_normalized = list()
    for sample_id, df in metat_df.groupby('sample_id', group_keys=True):
        df = methods[method](df, ref_gene_ids=ref_gene_ids)
        metat_df_normalized.append(df)

    return pd.concat(metat_df_normalized)


# Is it better to filter based on CPM or raw counts? Seems like CPM https://combine-australia.github.io/RNAseq-R/slides/RNASeq_filtering_qc.pdf
# Don't think it's worth scaling when looking at a single group of organisms, e.g. all the Methanoperedens. 
# Might be worth scaling when looking at something like the ECE because so many of the genes are so small. 
        
def metat_filter(metat_df, threshold:float=(5 / (1e8 / 1e6)), min_samples:int=8, field:str='read_count'):
    genome_id = metat_df.genome_id.values[0]
    # metat_df['cpm'] = metat_df.read_count / (metat_df.library_size / 1e6)
    mask = metat_df.groupby('gene_id').apply(lambda df : np.sum(df[field] >= threshold) >= min_samples, include_groups=False)
    keep_ids = mask.index[mask.values]
    # print(f'metat_filter: Keeping {len(keep_ids)} out of {metat_df['gene_id'].nunique()} total genes for {genome_id}.')
    metat_df = metat_df[metat_df['gene_id'].isin(keep_ids)].copy()
    return metat_df



def metat_group_genomes(metat_df:pd.DataFrame):
    '''Group the DataFrame values across genomes.'''
    if 'length' in metat_df.columns:
        metat_df['genome_size'] = metat_df.length # For consistency with coverm_df.

    # Aggregate the metat_df by genome ID. 
    agg_funcs = {'read_count':'sum', 'genome_size':'sum', 'library_size':'first', 'year':'first', 'location':'first', 'reactor':'first', 'detected':'mean', 'fraction_reads_in_top_genes':'first'}
    agg_funcs = {col:func for col, func in agg_funcs.items() if col in metat_df.columns} # For if the detected column isn't present. 

    metat_df = metat_df.groupby(['genome_id', 'sample_id']).agg(agg_funcs).reset_index()
    metat_df['rpkm'] = metat_df.read_count / (metat_df.genome_size / 1e3) / (metat_df.library_size / 1e6)
    metat_df['sample_id'] = metat_df.sample_id.str.replace('_metat', '')
    metat_df = metat_df.rename(columns={'detected':'fraction_detected'}) # The fraction of genes detected in the sample. 
    return metat_df



def metat_load(data_dir:str='../data/metat/', read_length:int=150, group_genomes:bool=False):
    metat_df = list()
    for path in glob.glob(os.path.join(data_dir, '*read_counts')):
        if 'summary' in path:
            continue
        df = pd.read_csv(path, sep='\t', comment='#', skiprows=2, header=None, names=['gene_id', 'contig_id', 'start', 'end', 'strand', 'length', 'read_count'])
        df['sample_id'] = os.path.basename(path).replace('_read_counts', '')
        metat_df.append(df)
    
    metat_df = pd.concat(metat_df).reset_index(drop=True)
    metat_df['genome_id'] = [gene_id.split('.')[0] for gene_id in metat_df.gene_id]
    metat_df['coverage'] = read_length * metat_df.read_count / metat_df.length
    metat_df = metat_add_library_size(metat_df) 
    metat_df['detected'] = metat_df.read_count > 0
    metat_df['location'] = [re.search('top|bottom|middle', sample_id).group(0) for sample_id in metat_df.sample_id]
    metat_df['reactor'] = [re.search('(ck|n)_', sample_id).group(1) for sample_id in metat_df.sample_id]
    metat_df['year'] = [re.search('2024|2025', sample_id).group(0) for sample_id in metat_df.sample_id]
    metat_df['sample_id'] = metat_df.sample_id.str.replace('_metat', '') # Redundant because this is all metaT data anyway.

    if group_genomes:
        metat_df = metat_group_genomes(metat_df)
    return metat_df


