import subprocess
import os
import argparse
from cut_pdb import list_pdb_files
from multiprocessing import Pool


def run_kpax_pairs(pivot_pdb, query_pdb, output_dir):
    #for docstring : pivot pdb is str, and query pdb is list
    # set log file name based on pivot name
    pivot_name = str(pivot_pdb).split('/')[-1].split('.')[0] + '.log'
    log_file = os.path.join(output_dir, pivot_name)
    n_prot = len(query_pdb)
    if os.path.exists(log_file) and os.stat(log_file).st_size > 370000:
        print('skip', flush = True)
        return
    command = ["kpax", "-nopdb", f"-top={n_prot}", f"-show={n_prot}", "-sort=T"] + [pivot_pdb] + query_pdb

    with open(log_file, 'w') as log:
        try:
            subprocess.run(command, check=True, stdout=log)
        except subprocess.CalledProcessError as e:
            print(f"Error while executing command: {e}")

def process(pdb_folder, output_dir):
    # given a folder, align all pdb files in folder pair by pair with kpax
    files = [str(pdb_posix) for pdb_posix in list_pdb_files(pdb_folder)]
    
    if not files:
        print(f"No files found in the folder: {pdb_folder}")
        return
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    pool_args = []

    # generate list of arguments (query, targets, output_dir) for multiprocessing of all pair-pair alignments
    for query in files:
        # copy list of all files and remove given pdb
        targets = files[:]
        targets.remove(query)
        pool_args.append((query, targets, output_dir))
    with Pool(processes=8) as pool:
        pool.starmap(run_kpax_pairs, pool_args)

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--pdb_folder", type=str)
    parser.add_argument("-o", "--output_dir", type=str)

    args = parser.parse_args()
    process(args.pdb_folder, args.output_dir)