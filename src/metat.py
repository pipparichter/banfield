import os
import skbio
import pandas as pd
import glob 
import numpy as np 

# Generate a script to submit the transcript-mapping jobs on Biotite. 


def metat_add_pseudocounts(metat_df:pd.DataFrame):
    '''Want to make sure to add pseudocounts based on everything in a sample, not just a particular organism!'''
    modified_df = list()
    for _, df in metat_df.groupby('sample_name', group_keys=False):
        df['read_count_original'] = df.read_count.values # Mark the read counts which were corrected. 
        n = df.read_count.sum()
        df['read_count'] = skbio.stats.composition.multi_replace(df.read_count.values) * n
        modified_df.append(df)
    return pd.concat(modified_df)


def metat_load(metat_dir:str='../data/metat', read_length:int=150):

    metat_df = list()
    for path in glob.glob(os.path.join(metat_dir, '*read_counts')):
        file_name = os.path.basename(path)
        sample_name, target_name = file_name.replace('_read_counts', '').split('-')
        df = pd.read_csv(path, sep='\t', comment='#')
        columns = df.columns.tolist()
        columns[-1] = 'read_count'
        df.columns = [col.lower() for col in columns]
        df = df.rename(columns={'end':'stop', 'chr':'contig_id', 'geneid':'gene_id'})
        df['sample_name'] = sample_name
        df['target_name'] = target_name
        metat_df.append(df)
    metat_df = pd.concat(metat_df).reset_index(drop=True)
    metat_df = metat_add_pseudocounts(metat_df)
    metat_df['coverage'] = read_length * metat_df.read_count / metat_df.length  
    # There should be one entry for each gene ID per sample. 
    # assert np.all(metat_df.value_counts(['sample_name', 'gene_id']) == 1), 'load_transcriptome_data: The gene IDs are not unique.'
    return metat_df


def _metat_load_summary(path:str):
    # The summary file only includes the mapped reads in the BAM file, so the total reads mapped to the provided reference genome. 
    # Therefore, these cannot be used to compute RPKM. 
    file_name = os.path.basename(path).replace('_read_counts.summary', '')
    sample_name, target_name = file_name.split('-')
    cols = dict()
    # cols['Status'] = 'bam_file'
    cols['Assigned'] = 'n_assigned'
    cols['Unassigned_Unmapped'] = 'n_unassigned_unmapped' # Could not be mapped to the reference. 
    cols['Unassigned_NoFeatures'] = 'n_unassigned_no_features' # Mapped, but could not be assigned to a feature. 
    cols['Unassigned_Ambiguity'] = 'n_unassigned_ambiguity'
    df = pd.read_csv(path, sep=r'\s+', index_col=0).T.reset_index(drop=True)
    df.columns.name = ''
    df = df[list(cols.keys())].copy()
    df = df.rename(columns=cols)
    df['target_name'] = target_name 
    df['sample_name'] = sample_name
    return df 

def metat_load_summary(data_dir='../data/metat/'):
    metat_summary_df = list()
    for path in glob.glob(os.path.join(data_dir, '*summary')):
        metat_summary_df.append(_metat_load_summary(path))
    metat_summary_df = pd.concat(metat_summary_df)
    metat_summary_df['total'] = metat_summary_df.n_unassigned_ambiguity + metat_summary_df.n_assigned + metat_summary_df.n_unassigned_no_features
    return metat_summary_df


