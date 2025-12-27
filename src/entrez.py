from Bio import Entrez 
import re 
from tqdm import tqdm
import pandas as pd 

Entrez.email = 'prichter@berkeley.edu'

_entrez_fix_key = lambda key : re.sub(r'\s|-', '_', key).replace('GBSeq_', '').lower()

def _entrez_protein_parse_result(result:dict):

    result = {key:value for key, value in dict(result).items() if ('GBReference' not in key) }# Don't think I need any of this data. 

    parsed_result = {_entrez_fix_key(key):value for key, value in result.items() if (type(value) == Entrez.Parser.StringElement)}

    feature_table = result['GBSeq_feature-table'][0]
    crossrefs = result.get('GBSeq_xrefs', [])
    parsed_result.update({_entrez_fix_key(key):value for key, value in feature_table.items() if (type(value) == Entrez.Parser.StringElement)})

    for qualifier in feature_table['GBFeature_quals']:
        parsed_result[qualifier['GBQualifier_name']] = qualifier.get('GBQualifier_value', True)
    for crossref in crossrefs:
        parsed_result[crossref['GBXref_dbname']] = crossref['GBXref_id']

    return parsed_result

    
def entrez_protein_get_metadata(protein_ids:list, path:str='../data/entrez.csv'):
    query_db = 'protein'
    
    df = list()
    for protein_id in tqdm(protein_ids, desc='entrez_protein_get_metadata'):
        handle = Entrez.efetch(id=protein_id, db=query_db, retmode='xml')
        result = list(Entrez.parse(handle))
        assert len(result) == 1, f'entrez_protein_get_metadata: Expected one result, got {len(result)}.'
        df.append(_entrez_protein_parse_result(result[0]))
    df = pd.DataFrame(df)
    df.to_csv(path)
    return df 

        