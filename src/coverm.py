import glob 
import pandas as pd 
import numpy as np 
import os 
import re 
from src.bbduk import bbduk_load
import importlib.resources
import io 


sample_metadata_path = importlib.resources.files('src').joinpath('data/sample_metadata.csv')
sample_metadata_df = pd.read_csv(io.StringIO(sample_metadata_path.read_text()))

ggkbase_id_to_sample_id_map = sample_metadata_df.set_index('ggkbase_id').sample_id.to_dict()
sample_id_to_ggkbase_id_map = sample_metadata_df.set_index('sample_id').ggkbase_id.to_dict()


def coverm_add_library_size(coverm_df, bbduk_data_dir='../data/bbduk'):
    bbduk_df = bbduk_load(bbduk_data_dir)
    bbduk_df = bbduk_df[~bbduk_df.index.str.contains('metat')].copy()
    coverm_df['library_size'] = coverm_df.sample_id.map(bbduk_df.library_size)
    return coverm_df


def coverm_group_genomes(coverm_df:pd.DataFrame):
    sample_df = coverm_df[coverm_df.sample_id == coverm_df.sample_id.values[0]].copy() # Get the coverage for a single sample. 
    coverm_df = coverm_df.copy()

    coverm_df['genome_size'] = coverm_df.genome_id.map(sample_df.groupby('genome_id').contig_size.sum())
    coverm_df['contig_weight'] = coverm_df.contig_size / coverm_df.genome_size 
    # Weight each of the columns by contig size. 
    coverm_df['variance'] = coverm_df['variance'] * coverm_df.contig_weight 
    coverm_df['trimmed_mean'] = coverm_df['trimmed_mean'] * coverm_df.contig_weight 
    coverm_df['mean'] = coverm_df['mean'] * coverm_df.contig_weight 

    agg_funcs = dict()
    agg_funcs['variance'] = 'mean'
    agg_funcs['mean'] = 'mean'
    agg_funcs['trimmed_mean'] = 'mean'
    agg_funcs['read_count'] = 'sum'
    agg_funcs['genome_size'] = 'first'
    agg_funcs['library_size'] = 'first'
    
    coverm_df = coverm_df.groupby(['sample_id', 'genome_id']).agg(agg_funcs)
    return coverm_df.reset_index()


def _coverm_load(df:pd.DataFrame):
    sample_id_pattern = r'/(.+)_trim_clean.PE.\d.fastq.gz'
    sample_ids = [re.search(sample_id_pattern, col).group(1) for col in df.columns if (re.search(sample_id_pattern, col) is not None)]
    sample_ids = np.unique(sample_ids)

    df = [df[['Contig'] + [col for col in df.columns if sample_id in col]].copy().assign(sample_id=sample_id) for sample_id in sample_ids]

    def rename_columns(df):
        col_map = {col:'_'.join(col.split(' ')[1:]).lower() for col in df.columns}
        col_map.update({'sample_id':'sample_id', 'Contig':'contig_id'})
        return df.rename(columns=col_map)

    df = [rename_columns(df_) for df_ in df]
    return pd.concat(df)

# variance: Coverage variance, which is a statistical measure of how uneven the read depth is across a contig or genome.
# trimmed_mean: Trimmed mean coverage per contig, which is essentially the mean depth after removing the lowest and highest 10% of per-base coverage values.
#   Read depth is defined per nucleotide base. 
# rpkm: Reads per kb per million mapped reads.
# tpm: TPM-normalized coverage, which is the RPKM for a specific contig divided by the sum of all RPKM. This controls for both contig length and sequencing depth. 

def coverm_load(path:str, bbduk_data_dir='../data/bbduk', contig_sizes:dict=None):
    coverm_df = _coverm_load(pd.read_csv(path, sep='\t'))
    coverm_df['genome_id'] = [contig_id.split('.')[0] for contig_id in coverm_df.contig_id]
    coverm_df['sample_id'] = coverm_df.sample_id.map(ggkbase_id_to_sample_id_map)

    if contig_sizes is not None:
        coverm_df['contig_size'] = coverm_df.contig_id.map(contig_sizes)
        coverm_df = coverm_add_library_size(coverm_df, bbduk_data_dir=bbduk_data_dir)
        coverm_df = coverm_group_genomes(coverm_df)
        coverm_df['rpkm'] = (coverm_df['read_count'] / (coverm_df['genome_size'] / 1e3)) / (coverm_df.library_size / 1e6) # RPKM is reads per kilobase of transcript per million mapped reads.
    
    sample_ids = coverm_df.sample_id.values
    coverm_df['location'] = [re.search('top|bottom|middle', sample_id).group(0)  if (re.search('top|bottom|middle', sample_id) is not None) else 'none' for sample_id in sample_ids]
    coverm_df['reactor'] = [re.search('(ck|n)_', sample_id).group(1)  if (re.search('(ck|n)_', sample_id) is not None) else 'none' for sample_id in sample_ids]
    coverm_df['year'] = [re.search('2024|2025', sample_id).group(0) if (re.search('2024|2025', sample_id) is not None) else 'none' for sample_id in sample_ids]

    return coverm_df