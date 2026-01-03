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

genome_ids = pd.read_csv('./data/genome_metadata.csv').genome_id.values
contig_sizes = dict()
for genome_id in genome_ids:
    contig_sizes.update(fasta_get_contig_sizes(f'../data/data/{genome_id}.fn'))

is_mp = lambda df : df.genome_id.str.startswith('mp_')
is_ece = lambda df : ~df.genome_id.str.startswith('mp_')

ece_id = 'linear_ece_19kb'
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


def load_interproscan(data_dir:str='../data/interproscan/', max_e_value:float=1e-1):
    interproscan_df = list()
    for path in glob.glob(os.path.join(data_dir, '*')):
        df = InterProScanFileTSV.from_file(path).to_df()
        df['genome_id'] = os.path.basename(path).replace('.tsv', '')
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