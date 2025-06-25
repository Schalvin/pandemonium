import pandas as pd
import re
import argparse
from pathlib import Path

def list_kpax_logs(kpax_folder):
    """
    Lists all .pdb files in the given folder using pathlib.
    
    Args:
        folder_path (str): The path to the folder containing the PDB files.
    
    Returns:
        list: A list of file paths for all .pdb files in the folder.
    """
    folder = Path(kpax_folder)
    # Use pathlib to find all .log files recursively
    return folder.rglob('*.log')

def process(kpax_folder, output_name):
    logs = list_kpax_logs(kpax_folder)
    d_outfile = open(f'{output_name}_all_directed.tsv', "w") # this file will contain all edges between all proteins (limited to 1000 edges per prot)
    a_outfile = open(f'{output_name}_filt_adjacency.tsv', "w") # this file will be used for clustering
    files_treated = 0
    pairs_treated= {}
    for log in logs:
        if str(log).split('/')[-1] == "kpax.log" : 
            # avoid double logs
            continue
        print(f'files treated :{files_treated}')
        files_treated += 1
        q_domain = str(log).split('/')[-1].split('.')[0]
        log_file = open(log, "r")
        for line in log_file:
            if line.startswith('Top'):
                q_length = int(line.split()[6].replace(')', ''))
                if q_length < 40 : # check if query is long enough
                    log_file.close()
                    break
            # check if line has n from 1 to 1000 in the begining (alignment rankings)
            if (re.match(r'^(?: {3}[1-9]| {2}[1-9]\d| {1}[1-9]\d{2}|1000)', line)):
                line = line.split()
                t_score = line[5]
                t_domain = line[14].split('[')[0]
                t_length = int(line[11])
                aligned_aa = int(line[7]) #number of aligned aa
                aligned_identity = int(line[9]) #number of identity aa for aligned seq
                if aligned_aa != 0: 
                    identity = aligned_identity*100/aligned_aa
                else: identity = 0
                # fill directed graph abc file with pair, t score, len and aligned sequence identity %
                # fill adjacency with only pairs having TM score > 0.45 and q and t length > 40
                # check if pair has already been seen and choose which TM score to keep (highest)
                if (t_domain, q_domain) in pairs_treated.keys() and pairs_treated[(t_domain, q_domain)][0] < t_score:
                    if float(t_score) > 0.45 and t_length > 40:
                        a_outfile.write(f'{q_domain} {t_domain} {t_score}\n')
                    d_outfile.write(f'{q_domain} {t_domain} {t_score} {identity}\n')
                    pairs_treated.pop((t_domain, q_domain))
                elif (t_domain, q_domain) in pairs_treated.keys():
                    if float(t_score) > 0.45 and t_length > 40:
                        a_outfile.write(f'{t_domain} {q_domain} {pairs_treated[(t_domain, q_domain)][0]}\n')
                    d_outfile.write(f'{t_domain} {q_domain} {t_score} {identity}\n')
                    pairs_treated.pop((t_domain, q_domain))
                else: pairs_treated[(q_domain, t_domain)] = (t_score, identity)
            else : continue
    for (q_domain, t_domain), (t_score, identity) in pairs_treated.items():
        if float(t_score) > 0.45 and t_length > 40:
            a_outfile.write(f'{q_domain} {t_domain} {t_score}\n')
        d_outfile.write(f'{q_domain} {t_domain} {t_score} {identity}\n')
        log_file.close()
    print(f'files treated :{files_treated}')
    a_outfile.close()
    d_outfile.close()
    return

            



if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--kpax_folder", type=str)
    parser.add_argument("-o", "--output_base_name", type=str)

    args = parser.parse_args()
    process(args.kpax_folder, args.output_base_name)