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
    bbduk_df = pd.DataFrame(bbduk_df)
    bbduk_df['location'] = [re.search('top|bottom|middle', sample_id).group(0)  if (re.search('top|bottom|middle', sample_id) is not None) else 'none' for sample_id in bbduk_df.sample_id]
    bbduk_df['reactor'] = [re.search('(ck|n)_', sample_id).group(1)  if (re.search('(ck|n)_', sample_id) is not None) else 'none' for sample_id in bbduk_df.sample_id]
    bbduk_df['year'] = [re.search('2024|2025', sample_id).group(0) if (re.search('2024|2025', sample_id) is not None) else 'none' for sample_id in bbduk_df.sample_id]
    return bbduk_df.set_index('sample_id')
