import pandas as pd 
import re 
import numpy as np 
import os 


class KofamscanFile():


    def __init__(self):
        pass 

    @classmethod
    def from_file(cls, path:str):
        cols = ['id', 'ko', 'threshold', 'score', 'e_value'] # , 'definition']
        # pattern = r'\s+(\d+_\d+)\s+(K\d{5})\s+([^\s]+)\s+([^\s]+)\s+([^\s]+)\s+"(.+)"'
        definition_pattern = r'"(.+)"'

        df = list()
        with open(path, 'r') as f:
            for line in f.readlines():
                if line.startswith('#'):
                    continue 
                definition = re.search(definition_pattern, line).group(1)
                line = line.replace('*', '')
                line = re.sub(definition_pattern, '', line)
                line = line.strip().split()
                line = dict(zip(cols, line))
                line['definition'] = definition
                df += [line]

        df = pd.DataFrame(df)
        df['e_value'] = pd.to_numeric(df.e_value)
        df['score'] = pd.to_numeric(df.score)
        df['threshold'] = pd.to_numeric(df.threshold)

        obj = cls()
        obj.df = df  
        return obj

    def to_df(self):
        return self.df.copy()
