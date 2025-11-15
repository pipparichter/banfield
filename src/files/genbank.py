import os 
import re
import pandas as pd 
from src.files.fasta import FASTAFile

feature_pattern = r'^\s{0,5}(\S+)\s+(\S+)'
seq_pattern = r'/translation="([^"]+)"'
product_pattern = r'/product="([^"]+)"'
# locus_tag_pattern = r'/locus_tag="([^"]+)"'


def genbank_read_features(path):
    content= ''
    with open(path, 'r') as f:
        read = False
        for line in f.readlines():
            if 'FEATURES' in line:
                read = True 
                continue 
            if 'ORIGIN' in line:
                read = False 
                break 
            if read:
                content += line
    return content


def genbank_parse_features(path:str):
    content = genbank_read_features(path)
    # feature_pattern = r'(?=^\s{0,5}(\S+)\s+(\S+))'
    features = re.split(feature_pattern, content, flags=re.MULTILINE)
    features = [feature for feature in features if len(feature) > 0]
    assert (len(features) % 3) == 0, 'GenBankFile: Expected the number of features to be divisible by three.'
    features = [(features[i], features[i + 1], features[i + 2]) for i in range(0, len(features), 3)]
    return features 
    

def genbank_parse_feature(feature):
    # feature = re.sub(r'[\s\n]{2,}', '', feature) # Remove any whitespace bigger than one. 
    info = dict()
    info['seq'] = re.search(seq_pattern, feature, flags=re.MULTILINE).group(1).replace('\n', '').replace(' ', '').replace('*', '')
    info['product'] = re.search(product_pattern, feature, flags=re.MULTILINE).group(1).replace('\n', '')
    # info['locus_tag'] = re.search(locus_tag_pattern, feature, flags=re.MULTILINE).group(1).replace('\n', '')
    return info 


class GenBankFile():

    aas = 'ACDEFGHIKLMNPQRSTVWYX*'
    seq_pattern = r'/transl="([^"]+)"'

    def __int__(self):
        pass 

    @classmethod
    def from_file(cls, path:str):
        
        features = genbank_parse_features(path)

        df = list()
        for feature_type, coordinate, feature in features:
            if feature_type != 'CDS':
                continue 
            row = {'feature_type':feature_type, 'coordinate':coordinate}
            row.update(genbank_parse_feature(feature))
            df.append(row)
        df = pd.DataFrame(df)
        df['id'] = [f'{os.path.basename(path).replace('.gbk', '')}_{i + 1}' for i in range(len(df))]
        df = df.set_index('id')
        
        obj = cls()
        obj.df = df.copy()
        obj.path = path 
        obj.file_name = os.path.basename(path)

        return obj 
    
    def to_fasta(self, path:str):
        fasta_file = FASTAFile.from_df(self.df)
        fasta_file.write(path)
    
    def to_df(self):
        return self.df.copy()

    # def to_fasta(path:str):

        