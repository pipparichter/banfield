# processes = dict()
# processes['methanotrophy'] = ['methyl+coenzyme', 'coM+methyltransferase', 'coenzyme+F430+synthase', 'methyl+coenzyme+M+reductase']
# processes['nitrate_reduction'] = ['narG', 'narH', 'nari', 'narK', 'narJ'] 
# processes['iron_reduction_pilus'] = ['archaeal+pilin']
# processes['iron_reduction_omc'] = ['omcA', 'omcC', 'omcB', 'omcE', 'omcF', 'mtrA', 'mtrB', 'mtrC']


database_patterns = dict()
database_patterns[r'PF\d+'] = 'pfam'
database_patterns[r'PTHR\d+'] = 'panther'
database_patterns[r'NF\d+'] = 'ncbifam'
database_patterns[r'IPR\d+'] = 'interpro'
database_patterns[r'PIRSF\d+'] = 'pirsf'
database_patterns[r'TIGR\d+'] = 'ncbifam'
database_patterns[r'cd\d+'] = 'cdd'
database_patterns[r'PS\d+'] = 'prosite'
database_patterns[r'G[\d\.]+'] = 'gene3d'
database_patterns[r'SSF\d+'] = 'ssf'
database_patterns[r'SFLD\d+'] = 'sfld'

def _interpro_get_metadata(accession:str):

    database = None
    for database_pattern, database in database_patterns.items():
        if re.match(database_pattern, accession) is not None:
            break

    info = dict()
    info['accession'] = accession
    info['database'] = database
    
    try:
        url = f'https://www.ebi.ac.uk/interpro/api/entry/{database}/{accession}'
        result = requests.get(url).text 
        result = json.loads(result)['metadata']
        info['name'] = result['name']['name']
        info['short_name'] = result['name']['short']
        info['description'] = result['description']['text']
        info['database'] = database
    except:
        pass 
    return info

def interpro_get_metadata(process, queries, n_hits:int=20):
    accessions = list()
    url = 'https://www.ebi.ac.uk/ebisearch/ws/rest/interpro?query={query}&size={n_hits}' # Reliability kind of drops off after the first 20 or so queries.
    for query in queries:
        result = requests.get(url.format(query=query, n_hits=n_hits)).text
        accessions += re.findall(r'id="([a-zA-Z0-9\.]+)"', result)
    accessions = np.unique(accessions) # Make sure there are no duplicates. 
    print(f'Obtained {len(accessions)} accessions for {process}')

    interpro_metadata_df = list()
    for accession in tqdm(accessions, desc='Collecting InterPro metadata.'):
        info = _interpro_get_metadata(accession)
        interpro_metadata_df.append(info)
    interpro_metadata_df = pd.DataFrame(interpro_metadata_df)
    interpro_metadata_df = interpro_metadata_df.drop_duplicates('accession')
    interpro_metadata_df.to_csv(f'../data/interpro_{process}.csv')

# for process, queries in processes.items():
#     interpro_get_metadata(process, queries)
# interpro_get_metadata('iron_reduction_pilus', processes['iron_reduction_pilus'])
# interpro_get_metadata('methanotrophy', processes['methanotrophy'])

genes = dict()
for process in processes:
    accessions = pd.read_csv(f'../data/interpro_{process}.csv', usecols=['accession']).accession.values 
    mask = interproscan_df.signature_accession.isin(accessions)
    genes[process] = interproscan_df[mask].gene_id.unique()

accessions = ['G3DSA:1.20.840.10']
interproscan_df = interproscan_df.sort_values('e_value') # Both alpha and beta subunits are covered by this annotation, only want to get the best hit (subunit A, I think)
mask = interproscan_df.signature_accession.isin(accessions) & ~interproscan_df.duplicated(['genome_id', 'signature_accession'], keep='first')
genes['methanotrophy_mcra'] = interproscan_df[mask].gene_id
