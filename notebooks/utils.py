import re 
import pandas as pd 
import os 
import glob


ggkbase_to_id_map = dict()
ggkbase_to_id_map['SR-VP_05_06_2024_ck_bottom_Methanoperedens_44_24'] = 'methanoperedens_1'
ggkbase_to_id_map['Final_SR-VP_05_06_2024_coassembly_19kb_linear_ECE_26_1334_complete'] = 'methanoperedens_2'
ggkbase_to_id_map['SR-VP_05_06_2024_N_top_Candidatus_Methanoperedens_Black-host_type_44_27'] = 'methanoperedens_2'


get_genus = lambda taxonomy : re.search('g__([^;]+)', taxonomy).group(1) if (re.search('g__([^;]+)', taxonomy) is not None) else 'none'
get_domain = lambda taxonomy : re.search('d__([^;]+)', taxonomy).group(1) if (re.search('d__([^;]+)', taxonomy) is not None) else 'none'
get_species = lambda taxonomy : re.search('s__([^;]+)', taxonomy).group(1) if (re.search('s__([^;]+)', taxonomy) is not None) else 'none'
get_order = lambda taxonomy : re.search('o__([^;]+)', taxonomy).group(1) if (re.search('o__([^;]+)', taxonomy) is not None) else 'none'
get_phylum = lambda taxonomy : re.search('p__([^;]+)', taxonomy).group(1) if (re.search('p__([^;]+)', taxonomy) is not None) else 'none'
get_family = lambda taxonomy : re.search('f__([^;]+)', taxonomy).group(1) if (re.search('f__([^;]+)', taxonomy) is not None) else 'none'
get_class = lambda taxonomy : re.search('c__([^;]+)', taxonomy).group(1) if (re.search('c__([^;]+)', taxonomy) is not None) else 'none'


def load_metadata(path:str='../data/bin_metadata.csv'):
    metadata_df = pd.read_csv(path)
    cols = [col for col in metadata_df.columns if ('relative_abundance' in col) or (col in ['genome_id', 'completeness', 'contamination', 'taxonomy'])]
    metadata_df = metadata_df[~metadata_df.taxonomy.isnull()].copy() # Not sure why some of these end up null. 
    metadata_df = metadata_df[cols].copy()
    metadata_df = metadata_df.melt(id_vars=['genome_id', 'completeness', 'contamination', 'taxonomy'], var_name='column_name', value_name='relative_abundance')
    metadata_df['nitrate'] = metadata_df.column_name.str.startswith('n')
    metadata_df['location'] = [re.search('(top|bottom|middle)', col).group(1) for col in metadata_df.column_name]

    metadata_df = metadata_df.set_index('genome_id')
    return metadata_df


def load_scaffold_to_bin_map(data_dir:str='../data/', cleaned:bool=True):
    scaffold_to_bin_map = list()
    for path in glob.glob(os.path.join(data_dir, f'*_scaffold_to_bin{"_cleaned" if cleaned else ""}.tsv')):
        scaffold_to_bin_map.append(pd.read_csv(path, sep='\t'))
    scaffold_to_bin_map = pd.concat(scaffold_to_bin_map).set_index('scaffold_name')['bin']
    scaffold_to_bin_map = scaffold_to_bin_map.to_dict()
    return scaffold_to_bin_map