import glob 
import pandas as pd 
import numpy as np 
import os 
import re 


def bbduk_load(data_dir='../data/bbduk'):
    bbduk_df = list()
    for path in glob.glob(os.path.join(data_dir, '*')):
        with open(path, 'r') as f:
            content = f.read()
        sample_id = os.path.basename(path).replace('.txt', '')
        total = re.search(r'#Total\s+(\d+)', content, flags=re.MULTILINE).group(1)
        bbduk_df.append({'sample_id':sample_id, 'library_size':int(total)})
    return pd.DataFrame(bbduk_df).set_index('sample_id')
