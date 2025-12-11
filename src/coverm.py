import glob 
import pandas as pd 
import numpy as np 
import os 
import re 
from src.bbduk import bbduk_load



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
    # agg_funcs['rpkm'] = 'sum'
    agg_funcs['variance'] = 'mean'
    agg_funcs['mean'] = 'mean'
    agg_funcs['trimmed_mean'] = 'mean'
    agg_funcs['count'] = 'sum'
    agg_funcs['genome_size'] = 'first'
    agg_funcs['library_size'] = 'first'
    
    coverm_df = coverm_df.groupby(['sample_id', 'target_name']).agg(agg_funcs)
    return coverm_df.reset_index()


# variance: Coverage variance, which is a statistical measure of how uneven the read depth is across a contig or genome.
# trimmed_mean: Trimmed mean coverage per contig, which is essentially the mean depth after removing the lowest and highest 10% of per-base coverage values.
#   Read depth is defined per nucleotide base. 
# rpkm: Reads per kb per million mapped reads.
# tpm: TPM-normalized coverage, which is the RPKM for a specific contig divided by the sum of all RPKM. This controls for both contig length and sequencing depth. 

def coverm_load(data_dir='../data/coverm/', bbduk_data_dir='../data/bbduk'):
    fields = 'mean trimmed_mean covered_bases variance count rpkm tpm'
    cols = fields.split()
    coverm_df = list()
    for path in glob.glob(os.path.join(data_dir, '*')):
        file_name = os.path.basename(path).replace('.tsv', '')
        sample_id, target_name = file_name.split('-')
        df = pd.read_csv(path, sep='\t', names=cols, skiprows=1)
        df['sample_id'], df['target_name'] = sample_id, target_name
        if len(df) > 0:
            coverm_df.append(df)
    coverm_df = pd.concat(coverm_df)
    coverm_df['contig_size'] = coverm_df.index.map(contig_sizes)

    coverm_df = coverm_add_library_size(coverm_df, bbduk_data_dir=bbduk_data_dir)
    coverm_df = coverm_group_targets(coverm_df)
    coverm_df['rpkm'] = (coverm_df['count'] / (coverm_df['genome_size'] / 1e3)) / (coverm_df.library_size / 1e6) # RPKM is reads per kilobase of transcript per million mapped reads. 
    return coverm_df