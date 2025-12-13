import re 
import pandas as pd 
import os 
import glob
from src.files import fasta_get_contig_sizes, fasta_get_genome_size
from src.coverm import * 
import numpy as np 


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
# id_to_ggkbase_name_map['jupiter_mini_borg_5'] = 'SR-VP_05_06_2024_coassembly_Jupiter_mini-Borg_38_3'
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

# id_to_ggkbase_name_map['phage_69kb'] = 'SR-VP_05_06_2024_coassembly_unknown_69_26'
# id_to_ggkbase_name_map['blue_borg'] = 'SR-VP_05_06_2024_coassembly_SR_VP_9_9_2021_65_4B_1_25m_Amethyst_Borg_complete_34_9'
# id_to_ggkbase_name_map['huge_phage_40kn'] = 'SR-VP_05_06_2024_coassembly_unknown_40_43'
# id_to_ggkbase_name_map['novel_mge'] = 'SR-VP_05_06_2024_coassembly_unknown_40_7'
# id_to_ggkbase_name_map['vulcan_1'] = 'SR-VP_05_06_2024_coassembly_unknown_42_6'
# id_to_ggkbase_name_map['vulcan_2'] = 'SR-VP_05_06_2024_coassembly_unknown_42_1'
# id_to_ggkbase_name_map['linear_ece_55kb'] = 'SR-VP_05_06_2024_coassembly_unknown_55_20'
# id_to_ggkbase_name_map['plasmid_like_41kb'] = 'SR-VP_05_06_2024_coassembly_Euryarchaeota_41_11'
# id_to_ggkbase_name_map['plasmid_like_39kb'] = 'SR-VP_05_06_2024_coassembly_Archaea_39_10'
# id_to_ggkbase_name_map['huge_phage_1mb'] = 'SR-VP_05_06_2024_coassembly_unknown_40_16'

# id_to_ggkbase_name_map['black_borg'] = 'SR-VP_05_06_2024_coassembly_BLACK_SR_VP_26_10_2019_C_40cm_scaffold_23_FINAL_IR_32_277'
id_to_ggkbase_name_map['black_borg'] = 'SR-VP_05_06_2024_coassembly_Black_Borg_32_272'
id_to_ggkbase_name_map['linear_ece_19kb'] = 'Final_SR-VP_05_06_2024_coassembly_19kb_linear_ECE_26_1334_complete'


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


def bbduk_load(data_dir='../data/bbduk'):
    bbduk_df = list()
    for path in glob.glob(os.path.join(data_dir, '*')):
        with open(path, 'r') as f:
            content = f.read()
        sample_id = os.path.basename(path).replace('.txt', '')
        total = re.search(r'#Total\s+(\d+)', content, flags=re.MULTILINE).group(1)
        bbduk_df.append({'sample_id':sample_id, 'library_size':int(total)})
    return pd.DataFrame(bbduk_df).set_index('sample_id')


def metat_add_library_size(metat_df, bbduk_data_dir='../data/bbduk'):
    bbduk_df = bbduk_load(bbduk_data_dir)
    bbduk_df = bbduk_df[bbduk_df.index.str.contains('metat')].copy()
    bbduk_df.index = bbduk_df.index.str.replace('_metat', '') # Remove the metat suffix from the sample name. 
    metat_df['library_size'] = metat_df.sample_id.map(bbduk_df.library_size)
    # RPKM is reads per kilobase of transcript per Million mapped reads
    # metat_df['rpkm'] = (metat_df['count'] / (metat_df['contig_size'] / 1e3)) / (metat_df.library_size / 1e6)
    return metat_df



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