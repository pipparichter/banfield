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

ggkbase_name_to_sample_id_map = sample_metadata_df.set_index('ggkbase_name').sample_id.to_dict()
sample_id_to_ggkbase_name_map = sample_metadata_df.set_index('sample_id').ggkbase_name.to_dict()


def coverm_add_library_size(coverm_df, bbduk_data_dir='../data/bbduk'):
    bbduk_df = bbduk_load(bbduk_data_dir)
    bbduk_df = bbduk_df[~bbduk_df.index.str.contains('metat')].copy()
    coverm_df['library_size'] = coverm_df.sample_id.map(bbduk_df.library_size)
    return coverm_df


def coverm_group_targets(coverm_df:pd.DataFrame):
    sample_df = coverm_df[coverm_df.sample_id == coverm_df.sample_id.values[0]].copy() # Get the coverage for a single sample. 
    coverm_df = coverm_df.copy()

    coverm_df['genome_size'] = coverm_df.target_name.map(sample_df.groupby('target_name').contig_size.sum())
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
    
    coverm_df = coverm_df.groupby(['sample_id', 'target_name']).agg(agg_funcs)
    return coverm_df.reset_index()

def coverm_reformat_columns(df:pd.DataFrame):
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

def coverm_load(data_dir='../data/coverm/', bbduk_data_dir='../data/bbduk', contig_size_df:pd.DataFrame=None, exclude_target_names=['mp_1', 'mp_2', 'mp_3', 'mp_4', 'mp_5']):
    coverm_df = list()
    for path in glob.glob(os.path.join(data_dir, '*')):
        file_name = os.path.basename(path).replace('.tsv', '')
        target_name = file_name.replace('.tsv', '')
        if target_name in exclude_target_names:
            continue

        df = coverm_reformat_columns(pd.read_csv(path, sep='\t'))
        df['target_name'] = target_name

        if target_name == 'mp':
            df['target_name'] = [contig_id.split('.')[0] for contig_id in df.contig_id]
            df['contig_id'] = [contig_id.split('.')[-1] for contig_id in df.contig_id]

        if len(df) > 0:
            coverm_df.append(df)
    coverm_df = pd.concat(coverm_df)
    coverm_df['sample_id'] = coverm_df.sample_id.map(ggkbase_name_to_sample_id_map)

    if contig_size_df is not None:
        coverm_df = coverm_df.merge(contig_size_df, on=['contig_id', 'target_name'], how='left')
        coverm_df = coverm_add_library_size(coverm_df, bbduk_data_dir=bbduk_data_dir)
        coverm_df = coverm_group_targets(coverm_df)
        coverm_df['rpkm'] = (coverm_df['read_count'] / (coverm_df['genome_size'] / 1e3)) / (coverm_df.library_size / 1e6) # RPKM is reads per kilobase of transcript per million mapped reads.

    return coverm_df