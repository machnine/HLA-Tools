'''
Download the latest protein alignment files from Github
'''
import os
from requests import request

GIT_URL = 'https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/alignments/'

HLA_LOCI = ['A','B','C','ClassI','DMA','DMB','DOA','DOB','DPA1','DPB1','DQA1', 
            'DQB1','DRA','DRB','E','F','G','HFE','MICA','MICB','TAP1','TAP2']

def download_latest_prot_alignments(locus_list, output_path='.'):
    prot_files = ['{}_prot.txt'.format(locus) for locus in locus_list]
    
    success_count = 0
    
    total_count = len(prot_files)
    
    for file in prot_files:
        try:
            r = request('GET', GIT_URL + file)
            with open(os.path.join(output_path, file), 'w') as fout:
                fout.write(r.text)
                
            success_count += 1
        except:
            pass

    print('{}/{} files downloaded successfully.'.format(success_count, 
                                                       total_count))
                                                       

if __name__ == '__main__':
    download_latest_prot_alignments(HLA_LOCI)
