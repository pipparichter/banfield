import os 
import re 
from typing import List, Dict, Tuple, NoReturn
from tqdm import tqdm 
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd 
import numpy as np 


def _parse_prodigal_description(description:str):
    pattern = r'# ([\d]+) # ([\d]+) # ([-1]+) # ID=([^;]+);partial=([^;]+);start_type=([^;]+);rbs_motif=([^;]+);rbs_spacer=([^;]+);gc_cont=([\.\w]+)'
    columns = ['start', 'stop', 'strand', 'ID', 'partial', 'start_type', 'rbs_motif', 'rbs_spacer', 'gc_content']
    match = re.search(pattern, description)
    parsed_header = {col:match.group(i + 1) for i, col in enumerate(columns)}
    parsed_header['rbs_motif'] = 'none' if (parsed_header['rbs_motif'] == 'None') else parsed_header['rbs_motif']
    return parsed_header 

get_contig_id = lambda id_ :  '_'.join(id_.split('_')[:-1])

get_reverse_complement = lambda seq : str(Seq(seq).reverse_complement())


class FASTAFile():
    # prodigal_dtypes = {'start':int, 'stop':int, 'strand':int, 'gc_content':float, 'rbs_motif':str, 'rbs_spacer':str}

    def __init__(self):
        '''Initialize a FASTAFile object.'''

        return


    def __len__(self):
        return len(self.seqs)
    
    @classmethod
    def from_df(cls, df:pd.DataFrame):
        obj = cls()

        obj.seqs = df.seq.values
        obj.ids = df.index.values 
        obj.descriptions = df.description.values if ('description' in df.columns) else [''] * len(obj.ids)
        return obj
        
    @classmethod
    def from_file(cls, path:str):
        obj = cls()
        f = open(path, 'r')
        obj.seqs, obj.ids, obj.descriptions = [], [], []
        for record in SeqIO.parse(path, 'fasta'):
            obj.ids.append(record.id)
            obj.descriptions.append(record.description.replace(record.id, '').strip())
            obj.seqs.append(str(record.seq))
        f.close()

        obj.seqs = np.array(obj.seqs)
        obj.ids = np.array(obj.ids)

        # assert len(obj.seqs) == len(np.unique(obj.seqs)), f'FASTAFile.from_file: Some of the sequence IDs in {path} are not unique.'
        return obj 
    
    def get_type(self):
        residues = ''.join(self.seqs)
        residues = np.unique(list(residues))
        if len(residues) < 20:
            return 'nt'
        else:
            return 'aa'
        
    def get_seq(self, id_:str, start:int=0, stop:int=None):
        assert id_ in self.ids, f'FASTAFile.get_seq: Sequence {id_} is not present.'

        seq = self.seqs[self.ids == id_][0]
        stop = len(seq) if (stop is None) else stop
        assert len(seq) >= stop, f'FASTAFile.get_seq: Specified stop {stop} is out of bounds for sequence of length {len(seq)}.' # Prevents an out-of-bound slice from failing silently.

        return seq[start:stop]
    
        
    def get_gc_content(self, exclude_unknown:bool=False, check:bool=False):
        if check:
            assert self.get_type() == 'nt', 'FASTAFile.get_gc_content: Needs to be a FASTA nucleotide file to compute GC content.'
        residues = ''.join(self.seqs)
        if exclude_unknown:
            residues = residues.replace('N', '')
        n = len(residues) # Total number of residues. 
        num_g = residues.count('G')
        num_c = residues.count('C')
        return (num_g + num_c) / n
            
    def to_df(self, parse_description:bool=False) -> pd.DataFrame:

        df = []
        for id_, seq, description in zip(self.ids, self.seqs, self.descriptions):
            row = {'description':description} 
            row['id'] = id_
            row['seq'] = seq
            if parse_description:
                row.update(_parse_prodigal_description(description))
                row['contig_id'] = get_contig_id(id_) # Extract contig ID. 
                row['start'], row['stop'], row['strand'] = int(row['start']), int(row['stop']), int(row['strand'])
            df.append(row)
            
        df = pd.DataFrame(df).set_index('id')


        return df

    def write(self, path:str, mode:str='w') -> NoReturn:
        f = open(path, mode=mode)
        records = []
        for id_, seq, description in zip(self.ids, self.seqs, self.descriptions):
            record = SeqRecord(Seq(seq), id=str(id_), description=description)
            records.append(record)
        SeqIO.write(records, f, 'fasta')
        f.close()





    # @staticmethod
    # def parse_coordinate(coord:str):
    #     strand = '-' if 'complement' in coord else '+'
    #     pattern = r'[<>]*(\d+)..(\d+)[<>]+'
    #     match_ = re.search(pattern, coord)

    #     if match_ is None:
    #         return None 
    #     else:
    #         return int(match_.group(1)), int(match_.group(2)), strand

    # def from_genbank(cls, path:str, ):

    #     with open(path, 'r') as f:
    #         content = f.read()

    #     contig_ids = [match_.group(1) for match_ in re.finditer('seqhdr="([^\s]+)')]
    #     contigs = re.split(r'^DEFINITION.+$', content, flags=re.MULTILINE)
    #     contigs = [re.sub('^FEATURES.+$', '', contig) for contig in contigs]

    #     for contig_id, contig in zip(contig_ids, contig):
    #         coords = 
        