import os
import skbio
import pandas as pd
import glob 
import numpy as np 
from src.bbduk import bbduk_load
from scipy.stats import gmean 


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
# (3) 
# TODO: Should I be doing this on a per-sample or overall basis? Seems like I should do it across samples if the goal is to compare that way...

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

    # assert 'read_count_original' not in metat_df.columns, 'metat_add_pseudocounts: Seems like the pseusocounts have already been added!'
    if 'read_count_original' in metat_df.columns:
        # print('metat_add_pseudocounts: Skipping pseudocounts, seems like they are already added.')
        return metat_df
    
    metat_df['read_count_original'] = metat_df.read_count.values.copy()
    metat_df_ = list()
    for _, df in metat_df.groupby('sample_id', group_keys=False):
        df = methods[method](df)
        metat_df_.append(df)
    return pd.concat(metat_df_)

metat_check_single_genome = lambda metat_df : metat_df.target_name.nunique() == 1
metat_check_single_sample = lambda metat_df : metat_df.sample_id.nunique() == 1

# Use the geometric mean (which does not work with zero values, so will need to adjust).
# Should I be using RPKM or read counts?
def _metat_normalize_alr(metat_df, ref_gene_ids:list=[]):
    target_name, sample_id = metat_df.target_name.values[0], metat_df.sample_id.values[0]
    mask = metat_df.gene_id.isin(ref_gene_ids) # & (metat_df.read_count_original > 0)
    # print('_metat_normalize_alr:', mask.sum(), f'nonzero reference genes for {target_name} in {sample_id}.')
    normalization_factor = gmean(metat_df[mask].read_count.values)
    # normalization_factor = gmean(metat_df[metat_df.gene_id.isin(ref_gene_ids)].read_count.values)
    metat_df['normalization_factor_alr'] = normalization_factor
    metat_df['read_count_0_normalized_alr'] = np.log(metat_df.read_count.min()) - np.log(normalization_factor)
    metat_df['read_count_normalized_alr'] = np.log(metat_df.read_count) - np.log(normalization_factor)
    return metat_df

def _metat_normalize_clr(metat_df, ref_gene_ids:list=[]):
    normalization_factor = gmean(metat_df.read_count.values)
    metat_df['normalization_factor_clr'] = normalization_factor
    metat_df['read_count_0_normalized_clr'] = np.log(metat_df.read_count.min()) - np.log(normalization_factor)
    metat_df['read_count_normalized_clr'] = np.log(metat_df.read_count) - np.log(normalization_factor)
    return metat_df

def metat_normalize(metat_df:pd.DataFrame, ref_gene_ids:dict=dict(), add_pseudocount:str='mzr', method:str='alr'):
    methods = dict()
    methods['clr'] = _metat_normalize_clr
    methods['alr'] = _metat_normalize_alr

    metat_df_normalized = list()
    for (sample_id, target_name), df in metat_df.groupby(['sample_id', 'target_name'], group_keys=True):
        df = metat_add_pseudocounts(df.copy(), method=add_pseudocount)
        df = methods[method](df, ref_gene_ids=ref_gene_ids.get(target_name, None))
        metat_df_normalized.append(df)

    return pd.concat(metat_df_normalized)


# Is it better to filter based on CPM or raw counts? Seems like CPM https://combine-australia.github.io/RNAseq-R/slides/RNASeq_filtering_qc.pdf
# Don't think it's worth scaling when looking at a single group of organisms, e.g. all the Methanoperedens. 
# Might be worth scaling when looking at something like the ECE because so many of the genes are so small. 
        
def metat_filter(metat_df, threshold:float=(5 / (1e8 / 1e6)), min_samples:int=8, field:str='cpm'):
    target_name = metat_df.target_name.values[0]
    # metat_df['cpm'] = metat_df.read_count / (metat_df.library_size / 1e6)
    mask = metat_df.groupby('gene_id').apply(lambda df : np.sum(df[field] >= threshold) >= min_samples, include_groups=False)
    keep_ids = mask.index[mask.values]
    print(f'filter_: Keeping {len(keep_ids)} out of {metat_df['gene_id'].nunique()} total genes for {target_name}.')
    metat_df = metat_df[metat_df['gene_id'].isin(keep_ids)].copy()
    return metat_df


def metat_load(metat_dir:str='../data/metat', read_length:int=150, add_pseudocounts:str=None):
    metat_df = list()

    for target_name in os.listdir(metat_dir): # Each of the subdirectories is the organism name. 
        assert os.path.isdir(os.path.join(metat_dir, target_name)), f'load_metat: Directory {os.path.join(metat_dir, target_name)} does not exist.'
        for path in glob.glob(os.path.join(metat_dir, target_name, '*read_counts')):
            if 'summary' in path:
                continue 
            df = pd.read_csv(path, sep='\t', comment='#', skiprows=2, header=None, names=['gene_id', 'contig', 'start', 'end', 'strand', 'length', 'read_count'])
            df['sample_id'] = os.path.basename(path).replace('_read_counts', '')
            df['target_name'] = target_name  
            metat_df.append(df)
    
    metat_df = pd.concat(metat_df).reset_index(drop=True)
    metat_df['coverage'] = read_length * metat_df.read_count / metat_df.length
    if add_pseudocounts is not None:
        metat_df = metat_add_pseudocounts(metat_df, method=add_pseudocounts)
    metat_df = metat_add_library_size(metat_df) 
    return metat_df



def _metat_load_summary(path:str):
    # The summary file only includes the mapped reads in the BAM file, so the total reads mapped to the provided reference genome. 
    # Therefore, these cannot be used to compute RPKM. 
    cols = dict()
    cols['Assigned'] = 'n_assigned'
    cols['Unassigned_Unmapped'] = 'n_unassigned_unmapped' # Could not be mapped to the reference. 
    cols['Unassigned_NoFeatures'] = 'n_unassigned_no_features' # Mapped, but could not be assigned to a feature. 
    cols['Unassigned_Ambiguity'] = 'n_unassigned_ambiguity'
    df = pd.read_csv(path, sep=r'\s+', index_col=0).T.reset_index(drop=True)
    df.columns.name = ''
    df = df[list(cols.keys())].copy()
    df = df.rename(columns=cols)
    return df 


def metat_load_summary(data_dir='../data/metat/'):
    metat_summary_df = list()
    for target_name in os.listdir(data_dir): # Each of the subdirectories is the organism name. 
        assert os.path.isdir(os.path.join(data_dir, target_name)), f'load_metat: Directory {os.path.join(data_dir, target_name)} does not exist.'
        for path in glob.glob(os.path.join(data_dir, target_name, '*summary')):
            df = _metat_load_summary(path)
            df['sample_id'] = os.path.basename(path).replace('_read_counts.summary', '')
            df['target_name'] = target_name  
            metat_summary_df.append(df)
    metat_summary_df = pd.concat(metat_summary_df)
    metat_summary_df['total'] = metat_summary_df.n_unassigned_ambiguity + metat_summary_df.n_assigned + metat_summary_df.n_unassigned_no_features
    metat_summary_df = metat_add_library_size(metat_summary_df)
    return metat_summary_df.reset_index(drop=True)

