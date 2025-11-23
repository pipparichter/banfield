import os

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