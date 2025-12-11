import requests 
import os 
import pandas as pd 
import numpy as np 
import re
from tqdm import tqdm


def kegg_load(path:str):
    '''Load in a KEGG annotation file.'''
    asterisk = list()
    cols = ['id', 'ko', 'threshold', 'score', 'e_value', 'definition']
    pattern = r'\s+([^\s]+)\s+(K\d{5})\s{1,}([^\s]+)\s{1,}([^\s]+)\s+([^\s]+)\s+(.+)'

    df = list()
    with open(path, 'r') as f:
        for line in f.readlines():
            if line.startswith('#'):
                continue 
            asterisk += [line.startswith('*')]
            line = line[1:]
            match_ = re.search(pattern, line)
            df += [dict([(cols[i], match_.group(i + 1)) for i in range(len(cols))])]
    df = pd.DataFrame(df)
    df['*'] = asterisk 
    df['definition'] = df.definition.str.strip()
    df['e_value'] = pd.to_numeric(df.e_value)
    return df 

def kegg_get_pathways_by_ko(kos, path:str='../data/kegg/pathways.tsv'):
    kos = [kos[i:i + 10] for i in range(0, len(kos), 10)] # Split into chunks to speed things up. 
    if not os.path.exists(path):
        content = ''
        for ko in tqdm(kos, desc='get_ko_pathways'):
            url = f'https://rest.kegg.jp/link/pathway/ko:{'+'.join(ko)}'
            content += requests.get(url).text
        with open(path, 'w') as f:
            f.write(content)
    df = pd.read_csv(path, sep='\t', names=['ko', 'pathway'])
    df['ko'] = df['ko'].str.replace('ko:', '')
    df['pathway'] = df['pathway'].str.replace('(ko:|path:)', '', regex=True)
    ko_to_pathway_map = df.groupby('ko').apply(lambda df : df['pathway'].tolist(), include_groups=False).to_dict()
    return ko_to_pathway_map


def kegg_get_modules_by_ko(kos, path:str='../data/kegg/modules.tsv'):
    kos = [kos[i:i + 10] for i in range(0, len(kos), 10)] # Split into chunks to speed things up. 
    if not os.path.exists(path):
        content = ''
        for ko in tqdm(kos, desc='get_ko_modules'):
            url = f'https://rest.kegg.jp/link/module/ko:{'+'.join(ko)}'
            content += requests.get(url).text
        with open(path, 'w') as f:
            f.write(content)
    df = pd.read_csv(path, sep='\t', names=['ko', 'module'])
    df['ko'] = df['ko'].str.replace('ko:', '')
    df['module'] = df['module'].str.replace('(ko:|md:)', '', regex=True)
    ko_to_module_map = df.groupby('ko').apply(lambda df : df['module'].tolist(), include_groups=False).to_dict()
    module_to_ko_map = df.groupby('module').apply(lambda df : df['ko'].tolist(), include_groups=False).to_dict()
    return ko_to_module_map, module_to_ko_map



def kegg_get_module(module, module_dir:str='../data/kegg/modules'):
    '''Download the summary file for a Kegg module and return the name.'''
    path = os.path.join(module_dir, f'{module}.txt')
    if not os.path.exists(path):
        url = f'https://rest.kegg.jp/get/{module}'
        content = requests.get(url).text
        with open(path, 'w') as f:
            f.write(content)
    with open(path, 'r') as f:
        content = f.read()
        name = re.search(r'^NAME\s+([^\n]+)$', content, flags=re.DOTALL|re.MULTILINE)
    if name is not None:
        return name.group(1)
    else:
        return 'none'