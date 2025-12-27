import requests 
import os 
import pandas as pd 
import numpy as np 
import re
from tqdm import tqdm
from src.files.kofamscan import KofamscanFile
import glob


def kofamscan_load(data_dir:str='../data/kofamscan'):
    kofamscan_df = list()
    for path in glob.glob(os.path.join(data_dir, '*')):
        target_name = os.path.basename(path).replace('.tsv', '')
        df = KofamscanFile.from_file(path).to_df()
        df['target_name'] = target_name
        kofamscan_df.append(df)
    kofamscan_df = pd.concat(kofamscan_df).reset_index(drop=True)
    return kofamscan_df

def kofamscan_get_pathways_by_ko(kos, path:str='../data/kofamscan/pathways.tsv'):
    kos = [kos[i:i + 10] for i in range(0, len(kos), 10)] # Split into chunks to speed things up. 
    if not os.path.exists(path):
        content = ''
        for ko in tqdm(kos, desc='get_ko_pathways'):
            url = f'https://rest.kofamscan.jp/link/pathway/ko:{'+'.join(ko)}'
            content += requests.get(url).text
        with open(path, 'w') as f:
            f.write(content)
    df = pd.read_csv(path, sep='\t', names=['ko', 'pathway'])
    df['ko'] = df['ko'].str.replace('ko:', '')
    df['pathway'] = df['pathway'].str.replace('(ko:|path:)', '', regex=True)
    ko_to_pathway_map = df.groupby('ko').apply(lambda df : df['pathway'].tolist(), include_groups=False).to_dict()
    return ko_to_pathway_map


def kofamscan_get_modules_by_ko(kos, path:str='../data/kofamscan/modules.tsv'):
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



def kofamscan_get_module(module, module_dir:str='../data/kofamscan/modules'):
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