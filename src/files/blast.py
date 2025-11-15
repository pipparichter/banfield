import re 
import json 
import pandas as pd 
import numpy as np 
import os 



class BLASTFileJSON():

    field_map = dict()
    field_map['accession'] = 'subject_id'
    field_map['query_title'] = 'id'
    field_map['title'] = 'subject_description'
    field_map['sciname'] = 'subject_taxon'
    field_map['taxid'] = 'subject_taxonomy_id'
    field_map['bit_score'] = 'bit_score'
    field_map['evalue'] = 'e_value'
    field_map['identity'] = 'identity'
    field_map['positive'] = 'positive'
    field_map['hit_from'] = 'subject_alignment_start'
    field_map['hit_to'] = 'subject_alignment_stop'
    field_map['query_from'] = 'query_alignment_start'
    field_map['query_to'] = 'query_alignment_stop'
    field_map['gaps'] = 'n_gaps'
    field_map['align_len'] = 'alignment_length'
    field_map['qseq'] = 'query_seq'
    field_map['hseq'] = 'subject_seq'
    field_map['len'] = 'subject_length'
    field_map['query_len'] = 'query_length'
    field_map['midline'] = 'alignment'

    def __init__(self):
        pass 

    @classmethod
    def from_file(cls, path):

        with open(path, 'r') as f:
            content = json.load(f)

        df = []
        for query in content['BlastOutput2']:
            results = query['report']['results']['search']
            query_info = {field:value for field, value in results.items() if (field != 'hits')}

            if (len(results['hits']) == 0):
                continue
            
            # Probably fine to just get the first high-scoring pair. There is generally just one anyway. 
            for hit in results['hits']:
                row = hit['description'][0].copy()
                row['len'] = hit['len']
                row.update(query_info)
                row.update(hit['hsps'][0]) # Get the best high-scoring pair. 
                df.append(row)

        df = pd.DataFrame(df)
        df = df[list(BLASTFileJSON.field_map.keys())].rename(columns=BLASTFileJSON.field_map) # Rename columns according to the field map. 
        
        obj = cls()
        obj.df = df 
        return obj
    

    def to_df(self):
        return self.df.copy()


                    # for hsp in hit['hsps']:
                    # row = query_info.copy()
                    # row.update(hit_info)
                    # row.update(hit['description'][0]) # Only get the description for the first entry. 
                    # row.update(hsp)
                    # df.append(row)
        # # Need to make sure not to count the number of high-scoring pairs, just the hits. 
        # n_hits = df.groupby('id').apply(lambda df : df.subject_id.nunique(dropna=True))
        # n_hits.name = 'n_hits'
        # n_subject_taxa = df.groupby('id').apply(lambda df : df.subject_taxonomy_id.nunique(dropna=True))
        # n_subject_taxa.name = 'n_subject_taxa'
        