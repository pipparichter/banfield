import re 
import pandas as pd 
import os 
import glob
from src.files import fasta_get_contig_sizes, fasta_get_genome_size

ggkbase_to_id_map = dict()
ggkbase_to_id_map['SR-VP_05_06_2024_ck_bottom_Methanoperedens_44_24'] = 'methanoperedens_1'
ggkbase_to_id_map['Final_SR-VP_05_06_2024_coassembly_19kb_linear_ECE_26_1334_complete'] = 'methanoperedens_2'
ggkbase_to_id_map['SR-VP_05_06_2024_N_top_Candidatus_Methanoperedens_Black-host_type_44_27'] = 'methanoperedens_2'

ece_id = 'ece_26_1334'
mp_id = 'methanoperedens_2'

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


contig_sizes = fasta_get_contig_sizes('../data/methanoperedens_2.fn')
contig_sizes.update(fasta_get_contig_sizes('../data/ece_26_1334.fn'))

genome_sizes = dict()
genome_sizes['methanoperedens_2']= fasta_get_genome_size('../data/methanoperedens_2.fn') 
genome_sizes['ece_26_1334'] = fasta_get_genome_size('../data/ece_26_1334.fn')


def bbduk_load(data_dir='../data/bbduk'):
    bbduk_df = list()
    for path in glob.glob(os.path.join(data_dir, '*')):
        with open(path, 'r') as f:
            content = f.read()
        sample_name = os.path.basename(path).replace('.txt', '')
        total = re.search(r'#Total\s+(\d+)', content, flags=re.MULTILINE).group(1)
        bbduk_df.append({'sample_name':sample_name, 'library_size':int(total)})
    return pd.DataFrame(bbduk_df).set_index('sample_name')


def metat_add_library_size(metat_df, bbduk_data_dir='../data/bbduk'):
    bbduk_df = bbduk_load(bbduk_data_dir)
    bbduk_df = bbduk_df[bbduk_df.index.str.contains('metat')].copy()
    bbduk_df.index = bbduk_df.index.str.replace('_metat', '') # Remove the metat suffix from the sample name. 
    metat_df['library_size'] = metat_df.sample_name.map(bbduk_df.library_size)
    # RPKM is reads per kilobase of transcript per Million mapped reads
    # metat_df['rpkm'] = (metat_df['count'] / (metat_df['contig_size'] / 1e3)) / (metat_df.library_size / 1e6)
    return metat_df


def coverm_add_library_size(coverm_df, bbduk_data_dir='../data/bbduk'):
    bbduk_df = bbduk_load(bbduk_data_dir)
    bbduk_df = bbduk_df[~bbduk_df.index.str.contains('metat')].copy()
    coverm_df['library_size'] = coverm_df.sample_name.map(bbduk_df.library_size)
    return coverm_df


def coverm_group_targets(coverm_df:pd.DataFrame):
    sample_df = coverm_df[coverm_df.sample_name == coverm_df.sample_name.values[0]].copy() # Get the coverage for a single sample. 
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
    
    coverm_df = coverm_df.groupby(['sample_name', 'target_name']).agg(agg_funcs)
    return coverm_df.reset_index()


def coverm_load(data_dir='../data/coverm/', bbduk_data_dir='../data/bbduk'):
    fields = 'mean trimmed_mean covered_bases variance count rpkm tpm'
    cols = fields.split()
    coverm_df = list()
    for path in glob.glob(os.path.join(data_dir, '*')):
        file_name = os.path.basename(path).replace('.tsv', '')
        sample_name, target_name = file_name.split('-')
        df = pd.read_csv(path, sep='\t', names=cols, skiprows=1)
        df['sample_name'], df['target_name'] = sample_name, target_name
        if len(df) > 0:
            coverm_df.append(df)
    coverm_df = pd.concat(coverm_df)
    coverm_df['contig_size'] = coverm_df.index.map(contig_sizes)

    coverm_df = coverm_add_library_size(coverm_df, bbduk_data_dir=bbduk_data_dir)
    coverm_df = coverm_group_targets(coverm_df)
    coverm_df['rpkm'] = (coverm_df['count'] / (coverm_df['genome_size'] / 1e3)) / (coverm_df.library_size / 1e6) # RPKM is reads per kilobase of transcript per million mapped reads. 
    return coverm_df


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