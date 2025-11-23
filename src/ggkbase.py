import subprocess 
import re
import os 


# https://ggkbase.berkeley.edu/organisms/668077/generate?type=genbank
# https://ggkbase.berkeley.edu/organisms/668075/download?type=genbank 

def scaffold_to_bin_remove_unknown(path:str=None):
    output_path = path.replace('.tsv', '_cleaned.tsv')
    if os.path.exists(output_path):
        return output_path
    
    cmd = f'grep -v "_UNK" {path} > {output_path}' 
    subprocess.run(cmd, shell=True, check=True)
    return output_path

def scaffold_to_bin_get_scaffold_ids(genome_id:str, path=None):
    genome_id = genome_id.replace('.', '_')
    cmd = f'grep "{genome_id}" {path}' # Get all rows in the file with the genome ID. 
    result = subprocess.run(cmd, shell=True,  capture_output=True)
    if result.returncode == 1:
        print(f'scaffold_to_bin_get_scaffold_ids: Scaffold IDs for genome {genome_id} not found in {path}')
        return []
    
    content = result.stdout.decode()
    scaffold_ids = [line.replace(genome_id, '').strip() for line in content.split('\n')]
    return scaffold_ids


def contigs_remove_unknown(path:str=None):
    output_path = path.replace('.fa', '_cleaned.fa')
    if os.path.exists(output_path):
        return output_path
    # This is actually so cool I didn't know you could do this with bash. Basically said to run 
    # the line of code in curly brackets whenever a FASTA header line (^>) is encountered. 
    # Then, if the input line ($0) does not contain UNK (i.e. does not match the pattern in between the forward slashes), 
    # then set the variable "keep" to 1. By default, awk prints the current line when the pattern (^>) evaluates to true, which would
    # mean that only the header lines are kept. However, here we override this behavior with a condition flag, which prints
    # whenever the condition (keep) is true.
    cmd = "awk '/^>/'" + r"'{keep = ($0 !~ /UNK/)} keep' " + f"{path} > {output_path}"
    subprocess.run(cmd, shell=True, check=True)
    return output_path


def contigs_get_genome(genome_id:str, contigs_path:str=None, genomes_dir='../data/genomes', scaffold_to_bin_path:str=None):
    genome_id = genome_id.replace('.', '_')

    output_path = os.path.join(genomes_dir, f'{genome_id}.fn')
    if os.path.exists(output_path):
        return output_path
    
    scaffold_ids = scaffold_to_bin_get_scaffold_ids(genome_id, path=scaffold_to_bin_path)
    if len(scaffold_ids) == 0:
        return None
    
    content, keep = '', False
    f = open(contigs_path, 'r')
    for line in f.readlines():
        match = re.search(r'>([^\s]+)', line) # Extract the scaffold ID, which is the non-whitespace stuff following the > character.
        if match is not None:
            keep = (match.group(1) in scaffold_ids)
        if keep:
            content += line
    f.close()
    
    with open(output_path, 'w') as f:
        f.write(content)