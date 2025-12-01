
def foldseek_submit_query(path:str):
    url = 'https://search.foldseek.com/api/ticket'
    with open(path, 'rb') as f:
        files = { 'q': f }
        data = [('mode', '3diaa'),
            ('database[]', 'BFVD'),
            ('database[]', 'afdb50'),
            ('database[]', 'afdb-proteome'),
            ('database[]', 'afdb-swissprot'),
            ('database[]', 'bfmd'),
            ('database[]', 'cath50'),
            ('database[]', 'mgnify_esm30'),
            ('database[]', 'pdb100'),
            ('database[]', 'gmgcl_id')]
        result = requests.post(url, files=files, data=data)
        if result.status_code != 200:
            print(f'foldseek_submit_query: Failed on {path}. {result.text}')
            return None
        return result.json()['id']
    
def foldseek_check_complete(job_id):
    url = f'https://search.foldseek.com/api/ticket/{job_id}'
    result = requests.get(url).json()
    return result['status'] == 'COMPLETE'


def foldseek_retrieve_results(job_id:str, output_path:str=None):
    url = f'https://search.foldseek.com/api/result/download/{job_id}'
    cmd = f'curl -o {output_path} {url}'
    subprocess.run(cmd, shell=True, check=True, stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)


for path in tqdm(paths, desc='foldseek_search'):
    output_file_name = os.path.basename(path).replace('.pdb', '.tar.gz')
    output_path = os.path.join(output_dir, output_file_name)
    if os.path.exists(output_path):
        continue

    job_id = foldseek_submit_query(path)
    if job_id is not None:
        while not foldseek_check_complete(job_id):
            time.sleep(10)
        foldseek_retrieve_results(job_id, output_path=output_path)


