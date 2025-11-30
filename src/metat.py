import os
import skbio
import pandas as pd
import glob 
import numpy as np 

# Generate a script to submit the transcript-mapping jobs on Biotite. 

sample_paths = list()
sample_paths = ['/groups/banfield/sequences/2025/SR-VP_Bioreactor_ck_bot_05_06_2024_metaT']
sample_paths += ['/groups/banfield/sequences/2025/SR-VP_Bioreactor_ck_bot_05_17_2025_metaT']
sample_paths += ['/groups/banfield/sequences/2025/SR-VP_Bioreactor_ck_mid_05_17_2025_metaT']
sample_paths += ['/groups/banfield/sequences/2025/SR-VP_Bioreactor_ck_top_05_17_2025_metaT']
sample_paths += ['/groups/banfield/sequences/2025/SR-VP_Bioreactor_N_bot_05_06_2024_metaT']
sample_paths += ['/groups/banfield/sequences/2025/SR-VP_Bioreactor_N_bot_05_17_2025_metaT']
sample_paths += ['/groups/banfield/sequences/2025/SR-VP_Bioreactor_N_top_05_06_2024_metaT']
sample_paths += ['/groups/banfield/sequences/2025/SR-VP_Bioreactor_N_top_05_17_2025_metaT']
sample_paths += ['/groups/banfield/sequences/2025/SR-VP_Bioreactor_N_mid_05_06_2024_metaT']
sample_paths += ['/groups/banfield/sequences/2025/SR-VP_Bioreactor_N_mid_05_17_2025_metaT']

sample_name_map = dict()
sample_name_map['SR-VP_Bioreactor_ck_bot_05_06_2024_metaT'] = 'ck_bottom_2024'
sample_name_map['SR-VP_Bioreactor_ck_bot_05_17_2025_metaT'] = 'ck_bottom_2025'
sample_name_map['SR-VP_Bioreactor_ck_mid_05_17_2025_metaT'] = 'ck_middle_2025'
sample_name_map['SR-VP_Bioreactor_ck_top_05_17_2025_metaT'] = 'ck_top_2025'
sample_name_map['SR-VP_Bioreactor_N_bot_05_06_2024_metaT'] = 'n_bottom_2024'
sample_name_map['SR-VP_Bioreactor_N_bot_05_17_2025_metaT'] = 'n_bottom_2025'
sample_name_map['SR-VP_Bioreactor_N_mid_05_06_2024_metaT'] = 'n_middle_2024'
sample_name_map['SR-VP_Bioreactor_N_mid_05_17_2025_metaT'] = 'n_middle_2025'
sample_name_map['SR-VP_Bioreactor_N_top_05_06_2024_metaT'] = 'n_top_2024'
sample_name_map['SR-VP_Bioreactor_N_top_05_17_2025_metaT'] = 'n_top_2025'

ece_id = 'ece_26_1334'
ref_paths = ['../data/methanoperedens_1.fn', '../data/methanoperedens_2.fn', f'../data/{ece_id}.fn']
output_dir = '../data/metat/'

def get_mapping_command(sample_path:str, ref_path:str=None, output_dir:str=output_dir):
    # TODO: Should look into what paired-end reads are and how that works experimentally. 

    target_name = os.path.basename(ref_path).replace('.fn', '')
    sample_name = os.path.basename(sample_path)
    input_path_1 = os.path.join(sample_path, 'raw.d', f'{sample_name}_trim_clean.PE.1.fastq.gz')
    input_path_2 = os.path.join(sample_path, 'raw.d', f'{sample_name}_trim_clean.PE.2.fastq.gz')
    output_path = os.path.join(output_dir, f'{sample_name_map[sample_name]}-{target_name}.bam')

    params = 'pigz=t unpigz=t ambiguous=random minid=0.96 idfilter=0.97 threads=64 out=stdout.sam editfilter=5 out=stdout.sam'
    cmd = f'bbmap.sh {params} in1={input_path_1} in2={input_path_2} ref={ref_path} nodisk | shrinksam | sambam > {output_path}'
    return cmd, output_path


def get_counting_command(bam_path:str, ref_path:str=None, output_dir:str=output_dir):
    output_file_name = os.path.basename(bam_path).replace('.bam', '')
    output_file_name += '_read_counts'
    output_path = os.path.join(output_dir, output_file_name)
    return f'featureCounts -p -T 64 -g ID -t CDS -a {ref_path} -s 2 -o {output_path} {bam_path}' 


def get_sbatch_command(cmd, job_name:str=None):
    return f'sbatch --wrap "{cmd}" --output ../slurm.out/{job_name}.out'


def metat_add_pseudocounts(metat_df:pd.DataFrame):
    '''Want to make sure to add pseudocounts based on everything in a sample, not just a particular organism!'''
    modified_metat_df = list()
    for _, df in metat_df.groupby('sample_name', group_keys=False):
        df['read_count_original'] = df.read_count.values # Mark the read counts which were corrected. 
        n = df.read_count.sum()
        df['read_count'] = skbio.stats.composition.multi_replace(df.read_count.values) * n
        modified_metat_df.append(df)
    return pd.concat(modified_metat_df)


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

# i = 0 

# mapping_path, counting_path = '../scripts/metat_mapping.sh', '../scripts/metat_counting.sh'
# mapping, counting = list(), list()
# for ref_path in ref_paths:
#     for sample_path in sample_paths:
#         mapping_command, bam_path = get_mapping_command(sample_path, ref_path=ref_path)
#         counting_command = get_counting_command(bam_path, ref_path=ref_path.replace('fn', 'gff'))
#         mapping.append(get_sbatch_command(mapping_command, job_name=i))
#         counting.append(counting_command)

#         i += 1

# with open(mapping_path, 'w') as f:
#     f.write('\n'.join(mapping))

# with open(counting_path, 'w') as f:
#     f.write('\n'.join(counting))