import re 
import pandas as pd 


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
