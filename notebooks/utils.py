import re 
import pandas as pd 
import os 
import glob
from src.files import *
from src.coverm import * 
import numpy as np 
from src.metat import *
from src.bbduk import bbduk_load
import seaborn as sns 
from scipy.stats import gmean 
from tqdm import tqdm
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from src.kofamscan import *
import requests 
import scipy 
import warnings 
import json
import itertools
import subprocess

id_to_ggkbase_name_map = {}
id_to_ggkbase_name_map['mp_4'] = 'SR-VP_11_27_2022_S1_80cm_Methanoperedens_44_5'
id_to_ggkbase_name_map['mp_1'] = 'SR-VP_05_06_2024_ck_bottom_Methanoperedens_44_47'
id_to_ggkbase_name_map['mp_3'] = 'SR-VP_05_06_2024_N_top_Methanoperedens_44_14'
id_to_ggkbase_name_map['mp_5'] = 'SR-VP_05_06_2024_N_top_Candidatus_Methanoperedens_Black-host_type_44_27' # Confirmed with CRISPR spacer. Same strain as other assembly with CRISPR hit.
id_to_ggkbase_name_map['mp_2'] = 'SR-VP_05_06_2024_ck_bottom_Methanoperedens_41_16'
# All Borgs in the Vernal Pool bioreactor coassembly. 
id_to_ggkbase_name_map['jupiter_mini_borg_1'] = 'SR-VP_05_06_2024_coassembly_Jupiter_mini-Borg_33_21'
id_to_ggkbase_name_map['jupiter_mini_borg_2'] = 'SR-VP_05_06_2024_coassembly_Jupiter_mini-Borg_32_5'
id_to_ggkbase_name_map['jupiter_mini_borg_3'] = 'SR-VP_05_06_2024_coassembly_Jupiter_mini-Borg_34_1246'
id_to_ggkbase_name_map['jupiter_mini_borg_4'] = 'SR-VP_05_06_2024_coassembly_Jupiter_mini-Borg_36_6'
id_to_ggkbase_name_map['jupiter_mini_borg_6'] = 'SR-VP_05_06_2024_coassembly_Jupiter_mini-Borg_35_3'
id_to_ggkbase_name_map['jupiter_mini_borg_7'] = 'SR-VP_05_06_2024_coassembly_Jupiter_mini-Borg_33_6'
id_to_ggkbase_name_map['jupiter_mini_borg_8'] = 'SR-VP_05_06_2024_coassembly_Jupiter_mini-Borg_31_4'
id_to_ggkbase_name_map['jupiter_mini_borg_9'] = 'SR-VP_05_06_2024_coassembly_Jupiter_mini-Borg_35_6'
id_to_ggkbase_name_map['saturn_mini_borg_1'] = 'SR-VP_05_06_2024_coassembly_Saturn_mini-Borg_35_7'
id_to_ggkbase_name_map['saturn_mini_borg_2'] = 'SR-VP_05_06_2024_coassembly_Saturn_mini-Borg_32_200'
id_to_ggkbase_name_map['saturn_mini_borg_3'] = 'SR-VP_05_06_2024_coassembly_Saturn_mini-Borg_33_3'
id_to_ggkbase_name_map['saturn_mini_borg_4'] = 'SR-VP_05_06_2024_coassembly_Saturn_mini-Borg_33_4'
id_to_ggkbase_name_map['unclassified_mini_borg'] = 'SR-VP_05_06_2024_coassembly_mini-Borg_reminscent_42_6-8'
id_to_ggkbase_name_map['unclassified_borg'] = 'SR-VP_05_06_2024_coassembly_new_Borg_34_11'
id_to_ggkbase_name_map['amethyst_borg'] = 'SR-VP_05_06_2024_coassembly_Amethyst_Borg_34_3'
id_to_ggkbase_name_map['oxblood_borg'] = 'SR-VP_05_06_2024_coassembly_Oxblood_Borg_33_20'
id_to_ggkbase_name_map['pink_borg'] = 'SR-VP_05_06_2024_coassembly_Pink_Borg_32_55'
id_to_ggkbase_name_map['purple_borg'] = 'SR-VP_05_06_2024_coassembly_Purple_Borg_33_3'
id_to_ggkbase_name_map['rose_borg'] = 'SR-VP_05_06_2024_coassembly_Rose_Borg_31_2'
id_to_ggkbase_name_map['vermilion_borg'] = 'SR-VP_05_06_2024_coassembly_Vermilion_Borg_34_8'
id_to_ggkbase_name_map['mercury_mini_borg'] = 'SR-VP_05_06_2024_coassembly_Mercury_mini-Borg_37_9'
id_to_ggkbase_name_map['saturn_mini_borg_like'] = 'SR-VP_05_06_2024_coassembly_Saturn_mini-Borg-like_32_7'
id_to_ggkbase_name_map['ruby_borg_related'] = 'SR-VP_05_06_2024_coassembly_Ruby-Borg-related_37_10'
id_to_ggkbase_name_map['black_borg'] = 'SR-VP_05_06_2024_coassembly_Black_Borg_32_272'
id_to_ggkbase_name_map['linear_ece_19kb'] = 'Final_SR-VP_05_06_2024_coassembly_19kb_linear_ECE_26_1334_complete'

is_mp = lambda df : df.target_name.str.startswith('mp_')
is_ece = lambda df : ~df.target_name.str.startswith('mp_')

ece_ggkbase_name = 'Final_SR-VP_05_06_2024_coassembly_19kb_linear_ECE_26_1334_complete'
bb_ggkbase_name = 'SR-VP_05_06_2024_coassembly_Black_Borg_32_272'

ece_id = 'ece_19kb'
mp_id = 'mp_5'

get_genus = lambda taxonomy : re.search('g__([^;]+)', taxonomy).group(1) if (re.search('g__([^;]+)', taxonomy) is not None) else 'none'
get_domain = lambda taxonomy : re.search('d__([^;]+)', taxonomy).group(1) if (re.search('d__([^;]+)', taxonomy) is not None) else 'none'
get_species = lambda taxonomy : re.search('s__([^;]+)', taxonomy).group(1) if (re.search('s__([^;]+)', taxonomy) is not None) else 'none'
get_order = lambda taxonomy : re.search('o__([^;]+)', taxonomy).group(1) if (re.search('o__([^;]+)', taxonomy) is not None) else 'none'
get_phylum = lambda taxonomy : re.search('p__([^;]+)', taxonomy).group(1) if (re.search('p__([^;]+)', taxonomy) is not None) else 'none'
get_family = lambda taxonomy : re.search('f__([^;]+)', taxonomy).group(1) if (re.search('f__([^;]+)', taxonomy) is not None) else 'none'
get_class = lambda taxonomy : re.search('c__([^;]+)', taxonomy).group(1) if (re.search('c__([^;]+)', taxonomy) is not None) else 'none'

level_funcs = dict()
level_funcs['genus'] = get_genus
level_funcs['domain'] = get_domain
level_funcs['species'] = get_species
level_funcs['order'] = get_order
level_funcs['phylum'] = get_phylum
level_funcs['family'] = get_family
level_funcs['class'] = get_class


contig_sizes = dict()
for ggkbase_name in id_to_ggkbase_name_map.values():
    contig_sizes.update(fasta_get_contig_sizes(f'../data/ggkbase/contigs/{ggkbase_name}.contigs.fa'))

genome_sizes = {id_:fasta_get_genome_size(f'../data/ggkbase/contigs/{ggkbase_name}.contigs.fa') for id_, ggkbase_name in id_to_ggkbase_name_map.items()}


def load_interproscan(data_dir:str='../data/interproscan/', max_e_value:float=1e-1):
    interproscan_df = list()
    for path in glob.glob(os.path.join(data_dir, '*')):
        df = InterProScanFileTSV.from_file(path).to_df()
        df['target_name'] = os.path.basename(path).replace('.tsv', '')
        interproscan_df.append(df)
    interproscan_df = pd.concat(interproscan_df)
    interproscan_df = interproscan_df[interproscan_df.e_value < max_e_value].copy()

    analysis_order = ['AntiFam', 'Gene3D','Pfam','SUPERFAMILY', 'TIGRFAM', 'FunFam', 'NCBIfam', 'SMART', 'SFLD', 'ProSitePatterns', 'ProSiteProfiles', 'PANTHER', 'Hamap', 'PRINTS', 'PIRSF', 'CDD']
    interproscan_df =  interproscan_df[interproscan_df.signature_analysis.isin(analysis_order)].copy()
    interproscan_df['signature_analysis'] = pd.Categorical(interproscan_df.signature_analysis, ordered=True, categories=analysis_order)
    interproscan_df = interproscan_df.sort_values('signature_analysis')

    fix_cdd_description = lambda description : description.replace('_', ' ') + ' domain-containing protein'
    interproscan_df['signature_description'] = np.where(interproscan_df.signature_analysis == 'CDD', interproscan_df.signature_description.apply(fix_cdd_description), interproscan_df.signature_description)
    interproscan_df['signature_description'] = interproscan_df.signature_description.str.replace('.', '', regex=False)
    interproscan_df['signature_description'] = np.where(interproscan_df.interpro_description != '-', interproscan_df.interpro_description, interproscan_df.signature_description)
    return interproscan_df


def load_metadata(path:str='../data/bin_metadata.csv'):
    metadata_df = pd.read_csv(path)
    cols = [col for col in metadata_df.columns if ('relative_abundance' in col) or (col in ['genome_id', 'completeness', 'contamination', 'taxonomy'])]
    metadata_df = metadata_df[~metadata_df.taxonomy.isnull()].copy() # Not sure why some of these end up null. 
    metadata_df = metadata_df[cols].copy()
    metadata_df = metadata_df.melt(id_vars=['genome_id', 'completeness', 'contamination', 'taxonomy'], var_name='column_name', value_name='relative_abundance')
    metadata_df['nitrate'] = metadata_df.column_name.str.startswith('n')
    metadata_df['location'] = [re.search('(top|bottom|middle)', col).group(1) for col in metadata_df.column_name]

    metadata_df = metadata_df.set_index('genome_id')

    for level, level_func in level_funcs.items():
        metadata_df[level] = metadata_df.taxonomy.apply(level_func)

    return metadata_df


def load_scaffold_to_bin_map(data_dir:str='../data/', cleaned:bool=True):
    scaffold_to_bin_map = list()
    for path in glob.glob(os.path.join(data_dir, f'*_scaffold_to_bin{"_cleaned" if cleaned else ""}.tsv')):
        scaffold_to_bin_map.append(pd.read_csv(path, sep='\t'))
    scaffold_to_bin_map = pd.concat(scaffold_to_bin_map).set_index('scaffold_name')['bin']
    scaffold_to_bin_map = scaffold_to_bin_map.to_dict()
    return scaffold_to_bin_map


get_metagenome_id = lambda id_ : re.search(r'(S[12]_\d+cm|N_middle|N_top|N_bottom|ck_bottom)', id_).group(0).lower()


def load_organism_info(path:str):
    metagenome_id = os.path.basename(path).replace('_organism_info.tsv', '')
    cols = ['id','name','taxonomy','bin length','GC%','coverage', 'curation status','completion status']
    df = pd.read_csv(path, sep='\t', usecols=cols)
    df.columns = [col.lower().replace(' ', '_') if (col != 'GC%') else 'gc_percent' for col in df.columns]
    df['metagenome_id'] = metagenome_id
    df = df[df.coverage > 0].copy()
    df = df.drop(columns=['id']).set_index('name')
    df['taxonomy'] = np.where(df.taxonomy.isnull(), 'none', df.taxonomy)
    return df 