import glob 
import os 
import pandas as pd
import numpy as np 
import re 
import subprocess

ref_genome_path = '/home/philippar/mp_x/mp_x.fn'
ref_genome_id = os.path.basename(ref_genome_path).split('.')[0]
metagenome_dir = '/home/philippar/data/ggkbase/metagenomes/'

min_percent_identity = 99
min_alignment_length = 100

def get_fasta_subset(handle, path:str=None, contig_ids:list=None):
    contig_id_pattern = re.compile('|'.join([f'({contig_id})' for contig_id in contig_ids]))
    with open(path, 'r') as f:
        read = False
        for line in f:
            if (line[0] == '>'):
                read = False
                contig_id = line.replace('>', '').split()[0]
                read = (re.fullmatch(contig_id_pattern, contig_id) is not None)
            if read:
                handle.write(line)

for metagenome_path in glob.glob(os.path.join(metagenome_dir, '*')):
    print(f'Running on metagenome {metagenome_path}')
    # metagenome_path = '../data/ggkbase/metagenomes/n_middle_2025.fn'
    metagenome_id = os.path.basename(metagenome_path).split('.')[0]

    blast_database_path = f'/home/philippar/db/{metagenome_id}'
    blast_output_dir = f'/home/philippar/{ref_genome_id}/blast/'
    blast_output_path = os.path.join(blast_output_dir, f'{ref_genome_id}-{metagenome_id}.tsv')
    blast_fields = 'qseqid sseqid pident qlen slen length mismatch gapopen qstart qend sstart send evalue bitscore stitle'

    output_path = os.path.join(os.path.dirname(ref_genome_path), f'{ref_genome_id}-{metagenome_id}.fn')

    if not os.path.exists(blast_database_path + '.ndb'):
        blast_cmd = f'makeblastdb -in {metagenome_path} -dbtype nucl -parse_seqids -out {blast_database_path}'
        subprocess.run(blast_cmd, shell=True, check=True)
    blast_cmd = f"blastn -query {ref_genome_path} -db {blast_database_path} -perc_identity 99.0 -out {blast_output_path} -outfmt \"6 {blast_fields}\""
    subprocess.run(blast_cmd, shell=True, check=True)

    blast_df = pd.read_csv(blast_output_path, sep='\t', header=None, names=blast_fields.split())
    # Want to filter for matches at the ends of scaffolds, either in the query or subject. If the match is in the middle, then the rest of the contig has lower identity. 
    mask = (blast_df.pident > min_percent_identity) & (blast_df.length > min_alignment_length) # (blast_df.evalue < 1e-10)

    delta = 2
    mask = mask & (((blast_df.qstart < delta) | ((blast_df.qlen - blast_df.qend) < delta)) | ((blast_df.sstart == 1) | ((blast_df.slen - blast_df.send) < delta)))
    blast_df = blast_df[mask].copy()
    contig_ids= blast_df.sseqid.unique()

    print(f'Num. matching contigs for {ref_genome_id} in {metagenome_id}:', len(contig_ids))
    print(f"Total genome size: {blast_df.drop_duplicates('sseqid').slen.sum() / 1e6:.3f} Mbp")
    # print(f'Genome GC content: {FASTAFile.from_file(output_path).get_gc_content(check=False) * 100:.2f}%')

    handle = open(output_path, 'w')
    get_fasta_subset(handle, path=metagenome_path, contig_ids=contig_ids)
    handle.close()
    print(f'Wrote new genome to {output_path}', end='\n\n')