import requests
import os 
import pandas as pd 
import numpy as np 
import re 
from src.files import FASTAFile
from tqdm import tqdm 
import io
import json 


esm_valid_tokens = {'C', 'D', 'T', 'X', 'E', 'Z', 'H', 'M', 'N', 'L', 'S', 'A', 'P', 'G', 'V', 'W', 'F', 'R', 'Y', 'B', 'J', 'K', 'Q', 'I'}

def fold_esm(path:str, output_dir:str='../data/structures/esmfold/'):
    '''
    Note that ESMFold server will not fold proteins larger than 400 amino acids. 
    
    :param path: Path to the FASTA file containing sequences to fold. 
    :param output_dir: Path to the directory where the structures will be stored. 
    '''

    url = 'https://api.esmatlas.com//foldSequence/v1/pdb/'

    fasta_file = FASTAFile.from_file(path)
    for id_, sequence in tqdm(list(zip(fasta_file.ids, fasta_file.seqs)), desc='fold_esm'):
        path = os.path.join(output_dir, f'{id_}.pdb')
        if os.path.exists(path):
            continue
        sequence = sequence.replace('*', '').strip()
        # valid_tokens = np.array([(aa in esm_valid_tokens) for aa in sequence])
        # invalid_tokens = np.array(list(sequence))[~valid_tokens].tolist()
        # assert np.all(valid_tokens), f'fold_esm: Invalid tokens in the sequence, {', '.join(invalid_tokens)}.'
        result = requests.post(url, data=sequence)
        try:
            assert result.status_code == 200, f'fold_esm: Error in ESM-folding sequence {id_}'
            with open(path, 'w') as f:
                f.write(result.text)
        except:
            # print(result.text)
            print(f'{id_}: {sequence}')

# https://www.rbvi.ucsf.edu/chimerax/data/pae-apr2022/pae.html
# PAE is predicted aligned error, and is a measure of the relative position of residue i to residue j; this is an assessment of inter-domain relative positioning.
# It is an output the model learns during training by aligning the predicted position of residue i to the true structure, and then getting
# the error in the predicted distances between i and all other residues. The model directly outputs a distribution of aligned errors for each residue.

# pLDDT is a local confidence score, and is a measure of the accuracy of a residue's local position with its immediate contacts. In this case, the model outputs
# a distribution of DISTANCES, not errors, and it converts the predicted distribution into a confidence score based on the variance.

# Because PAE relies on supervised training, ESMfold doesn't really have a good way to directly measure PAE. It seems like it can be approximated by 
# ensemble-based error, which might be how they get the PAE for the predicted structures already in the database. 

def fold_esm_load_pdb(path:str):
    with open(path, 'r') as f:
        lines = f.readlines()
    lines = [re.sub(r'ATOM\s+', '', line) for line in lines if (line.startswith('ATOM'))]

    # Confidence is referred to as the "B-factor," not totally sure why. 
    columns = ['atom_number', 'atom_name', 'residue_name', 'chain_id', 'residue_number', 'x', 'y', 'z', 'occupancy', 'confidence', 'element']
    df = pd.read_csv(io.StringIO(''.join(lines)), sep=r'\s+', names=columns, header=None)

    # Overall residue-level confidence score is typically the average of all confidence scores in the B-factor column. 
    # ESM and AlphaFold overwrite the B-factor with the per-residue confidence score by default. 
    return df

def fold_esm_get_confidence(path:str):
    pdb_df = fold_esm_load_pdb(path)
    plddts = pdb_df.groupby('residue_number').confidence.mean()
    return plddts.mean(), (plddts > 0.5).mean()

# https://app.gitbook.com/o/-LzcB3BNVSNh_20MBLKi/s/-M-S5z_vnCqDzDHcfsmi/protein-folding
# sbatch --partition gpu --gpus 1 -is the -wrap "/shared/software/bin/colabfold_batch --use-gpu-relax --amber --templates --num-recycle 3 --model_type monomer ece_26_1334.fa alphafold" --output ./slurm-colabfold.out

# sbatch --wrap "rosettafold2 -o rosettafold/ ece_26_1334.fa" --gres gpu:1 --partition gpu --output slurm-rosettafold.out

def fold_alphafold_make_input_file(input_path:str, path:str='../data/ece_26_1334_alphafold.json'):
    '''Convert a FASTA file to a JSON file to use as AlphaFold input.'''
    fa_df = FASTAFile.from_file(input_path).to_df(parse_description=False)
    content = list()
    for row_ in fa_df.itertuples():
        row = {'modelSeeds':[], 'version':1, 'dialect':'alphafoldserver'}
        row['name'] = row_.Index 
        row['sequences'] = [{'proteinChain':{'sequence':row_.seq, 'count':1, 'useStructureTemplate': False}}]
        content.append(row)
    with open(path, 'w') as f:
        json.dump(content, f)
    return path

def fold_alphafold_make_input_directory(input_path:str, dir_:str='../data/ece_26_1334/'):
    if not os.path.isdir(dir_):
        os.mkdir(dir_)

    fasta_file = FASTAFile.from_file(input_path)
    for id_, sequence in zip(fasta_file.ids, fasta_file.seqs):
        path = os.path.join(dir_, f'{id_}.fa')
        with open(path, 'w') as f:
            f.write(f'>{id_}\n{sequence}')




def fold_make_rosettafold_script(output_dir:str=None):
    pass 
