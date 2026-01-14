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

# library_sizes = dict()
# library_sizes['']


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


def plot_diff_distributions(metat_df:pd.DataFrame, locations=['bottom', 'middle', 'top'], threshold=5, normalization:str='alr', genome_id=None, axes=None, color='lightgray', ref_gene_ids:list=[]):

    for ax, location in zip(axes, locations):
        figure_df = metat_get_diff(metat_df, genome_id=genome_id, location=location, normalization=normalization, ref_gene_ids=ref_gene_ids, threshold=threshold)
        
        sns.kdeplot(figure_df, x='diff', common_norm=False, color=color, ax=ax, label=threshold)
        
        ax.set_ylabel('density')
        ax.set_title(location)

        # if quantile is not None:
        #     up_threshold = np.quantile(figure_df['diff'].values, 1 - quantile).item()
        #     down_threshold = np.quantile(figure_df['diff'].values, quantile).item()
        #     ax.axvline(up_threshold, ls='--', color='black', lw=0.7)
        #     ax.axvline(down_threshold, ls='--', color='black', lw=0.7)

def plot_read_counts(metat_df:pd.DataFrame, title='ribosomal proteins', drop_empty:bool=True, figsize=(3, 5), reactor:str='n'):

    # assert 'annotation' in metat_df.columns, 'plot_read_counts: Missing required column "annotation"'
    assert 'read_count' in metat_df.columns, 'plot_read_counts: Missing required column "read_count"'
    assert 'sample_id' in metat_df.columns, 'plot_read_counts: Missing required column "sample_id"'

    fig, ax = plt.subplots(figsize=figsize)
    
    figure_df = metat_df[metat_df.reactor == reactor].copy()
    figure_df = figure_df.pivot(columns='sample_id', values='read_count', index='gene_id')
    figure_df.columns = figure_df.columns.str.replace('_metat', '')
    figure_df = figure_df.fillna(0).astype(int)
    if drop_empty:
        print(f'plot_read_counts: Dropping {(figure_df.sum(axis=1).values.ravel() > 0).sum()} genes with no presence in any sample.')
        figure_df = figure_df[figure_df.sum(axis=1).values.ravel() > 0].copy()

    sns.heatmap(figure_df, annot=True, fmt='d', cmap='Grays', lw=1, cbar=False)
    ax.set_ylabel('')
    ax.set_xlabel('')
    ax.set_title(title)

    if 'annotation' in metat_df.columns:
        annotations = metat_df.set_index('gene_id').annotation.to_dict()
        gene_ids = [label.get_text() for label in ax.get_yticklabels()]
        y_tick_labels = [f'{gene_id} {annotations.get(gene_id, gene_id)}' for gene_id in gene_ids]
    else:
        y_tick_labels = ax.get_yticklabels()
    ax.set_yticks(ax.get_yticks(), y_tick_labels, rotation=0)
    plt.show()


mcr_accessions = dict()
mcr_accessions['IPR015823'] = 'methyl-coenzyme M reductase A/B'
mcr_accessions['IPR008924'] = 'methyl-coenzyme M reductase A/B'
mcr_accessions['IPR015811'] = 'methyl-coenzyme M reductase A/B'
mcr_accessions['IPR022680'] = 'methyl-coenzyme M reductase B' # Beta subunit. 
mcr_accessions['IPR003179'] = 'methyl-coenzyme M reductase B' # Beta subunit. 
mcr_accessions['IPR022679'] = 'methyl-coenzyme M reductase B' # Beta subunit. 
mcr_accessions['IPR009047'] = 'methyl-coenzyme M reductase A' # Alpha subunit.
mcr_accessions['IPR003183'] = 'methyl-coenzyme M reductase A' # Alpha subunit. 
mcr_accessions['IPR016212'] = 'methyl-coenzyme M reductase A' # Alpha subunit.  

met_accessions = mcr_accessions.copy()
met_accessions['IPR015823'] = 'methyl-coenzyme M reductase'
met_accessions['IPR008924'] = 'methyl-coenzyme M reductase'
met_accessions['IPR015811'] = 'methyl-coenzyme M reductase'
met_accessions['IPR022680'] = 'methyl-coenzyme M reductase' # Beta subunit. 
met_accessions['IPR003179'] = 'methyl-coenzyme M reductase' # Beta subunit. 
met_accessions['IPR022679'] = 'methyl-coenzyme M reductase' # Beta subunit. 
met_accessions['IPR009047'] = 'methyl-coenzyme M reductase' # Alpha subunit.
met_accessions['IPR003183'] = 'methyl-coenzyme M reductase' # Alpha subunit. 
met_accessions['IPR016212'] = 'methyl-coenzyme M reductase' # Alpha subunit. 
# Tetrapyrrole methylase has two non-similar subdomains. 
# Includes uroporphyrinogen III methyltransferase and diphthine synthase.
# These enzymes catalyse the methylation of their substrates using S-adenosyl-L-methionine as a methyl source.
met_accessions['IPR014776'] = 'tetrapyrrole methylase' 
met_accessions['IPR014777'] = 'tetrapyrrole methylase'
# met_accessions['IPR035996'] = 'tetrapyrrole methylase' # Superfamily. 
met_accessions['IPR000878'] = 'tetrapyrrole methylase'
met_accessions['IPR015421'] = 'pyridoxal phosphate-dependent transferase'
met_accessions['IPR015422'] = 'pyridoxal phosphate-dependent transferase'
met_accessions['IPR036108'] = 'uroporphyrinogen III synthase' # Involved in tetrapyrrole biosynthesis.
met_accessions['IPR003754'] = 'uroporphyrinogen III synthase'
met_accessions['IPR039793'] = 'uroporphyrinogen III synthase'
met_accessions['IPR006366'] = 'uroporphyrin-III C-methyltransferase'
# met_accessions['IPR036803'] = 'porphobilinogen deaminase' # Superfamily. 
# met_accessions['IPR022417'] = 'Porphobilinogen deaminase' # Superfamily. 
met_accessions['IPR022418'] = 'porphobilinogen deaminase'
met_accessions['IPR000860'] = 'porphobilinogen deaminase'
# met_accessions['IPR036343'] = 'glutamyl-tRNA reductase' # Superfamily. 
# met_accessions['IPR036453'] = 'glutamyl-tRNA reductase' # Superfamily. 
met_accessions['IPR015895'] = 'glutamyl-tRNA reductase' 
met_accessions['IPR000343'] = 'glutamyl-tRNA reductase' 
met_accessions['IPR006151'] = 'glutamyl-tRNA reductase' 
met_accessions['IPR015896'] = 'glutamyl-tRNA reductase' # Involved in tetrapyrrole biosynthesis.
met_accessions['IPR001731'] = 'delta-aminolevulinic acid dehydratase'
met_accessions['IPR004551'] = 'diphthine synthase'
met_accessions['IPR012818'] = 'precorrin methyltransferase'
met_accessions['IPR012382'] = 'precorrin methyltransferase'
met_accessions['IPR051810'] = 'precorrin methyltransferase'
met_accessions['IPR006363'] = 'precorrin methyltransferase' 
met_accessions['IPR004639'] = 'glutamate-1-semialdehyde aminotransferase'

rnap_accessions = dict()
rnap_accessions['IPR037034'] = 'RNA polymerase Rpo2' # Superfamily. 
rnap_accessions['IPR042102'] = 'RNA polymerase Rpo1' # Superfamily. 
rnap_accessions['IPR044893'] = 'RNA polymerase Rpo1' # Superfamily.
rnap_accessions['IPR038120'] = 'RNA polymerase Rpo1' # Superfamily. 
rnap_accessions['IPR036603'] = 'RNA polymerase Rpo11' # Rbp11-like subunit
# rnap_accessions['IPR036643'] = 'RNA polymerase'
rnap_accessions['IPR035913'] = 'RNA polymerase Rpo5' # Rbp5-like, superfamily. 
# rnap_accessions['IPR037033'] = 'RNA polymerase'
rnap_accessions['IPR014724'] = 'RNA polymerase Rpo2,'
rnap_accessions['IPR036898'] = 'RNA polymerase Rpo7' # Rbp7-like, superfamily. 
rnap_accessions['IPR007644'] = 'RNA polymerase Rbo1B' # Beta subunit. 
rnap_accessions['IPR007120'] = 'RNA polymerase Rbp2'
rnap_accessions['IPR007642'] = 'RNA polymerase Rpo2'
rnap_accessions['IPR006110'] = 'RNA polymerase Rpo6'
rnap_accessions['IPR007646'] = 'RNA polymerase Rpo2'
rnap_accessions['IPR007641'] = 'RNA polymerase Rpo2'
rnap_accessions['IPR007645'] = 'RNA polymerase Rpo2'
rnap_accessions['IPR000783'] = 'RNA polymerase Rpo5'
rnap_accessions['IPR011263'] = 'RNA polymerase Rpo3'
# rnap_accessions['IPR011262'] = 'RNA polymerase'
rnap_accessions['IPR009025'] = 'RNA polymerase Rbp11'
rnap_accessions['IPR007081'] = 'RNA polymerase Rpo1'
rnap_accessions['IPR000722'] = 'RNA polymerase Rpo1A' # Alpha subunit.
rnap_accessions['IPR007066'] = 'RNA polymerase Rpo1'
rnap_accessions['IPR007083'] = 'RNA polymerase Rpo1'
rnap_accessions['IPR007080'] = 'RNA polymerase Rpo1'
rnap_accessions['IPR005574'] = 'RNA polymerase Rpo4'
rnap_accessions['IPR005576'] = 'RNA polymerase Rpo7' # Rbp7-like
rnap_accessions['IPR000268'] = 'RNA polymerase Rbp10'
# rnap_accessions['IPR043502'] = 'DNA/RNA polymerase superfamily'
rnap_accessions['IPR029040'] = 'RNA polymerase RpoK'
rnap_accessions['IPR023580'] = 'RNA polymerase Rpo10'
rnap_accessions['IPR019969'] = 'RNA polymerase Rpo2'
rnap_accessions['IPR004519'] = 'RNA polymerase Rpo8'
rnap_accessions['IPR012757'] = 'RNA polymerase Rpo1'
rnap_accessions['IPR012758'] = 'RNA polymerase Rpo1'
# rnap_accessions['IPR006592'] = 'RNA polymerase'
rnap_accessions['IPR001529'] = 'RNA polymerase Rpo10' # Also RpoM.
rnap_accessions['IPR006591'] = 'RNA polymerase Rpo12' # Also RpoP.
rnap_accessions['IPR045867'] = 'RNA polymerase Rpo1B' # Subunit beta-prime.
rnap_accessions['IPR015712'] = 'RNA polymerase Rpo2'
rnap_accessions['IPR010924'] = 'RNA polymerase Rpo4'
rnap_accessions['IPR014381'] = 'RNA polymerase Rpo5'
rnap_accessions['IPR012164'] = 'RNA polymerase Rpo11' # Also RpoS.
rnap_accessions['IPR045113'] = 'RNA polymerase Rpo7'
rnap_accessions['IPR006111'] = 'RNA polymerase Rpo6'