import pandas as pd
import re
import os
import argparse

from pathlib import Path
from Bio import PDB

# PDB_FOLDER = "/home/chili/chalvin/projet_long/ATLAS/prot_pdb"
# ATLAS_TB = "/home/chili/chalvin/projet_long/ATLAS/2023_03_09_ATLAS_info.tsv"
# OUTPUT_DIR = "/home/chili/chalvin/projet_long/ATLAS/SWORDdomains_pdb"


def list_pdb_files(folder_path):
    """
    Lists all .pdb files in the given folder using pathlib.
    
    Args:
        folder_path (str): The path to the folder containing the PDB files.
    
    Returns:
        list: A list of file paths for all .pdb files in the folder.
    """
    folder = Path(folder_path)
    # Use pathlib to find all .pdb files recursively
    return folder.rglob('*.pdb')

def calc_overlap(start1, end1, start2, end2):
    overlap = max(0, min(end1, end2) - max(start1, start2))
    return overlap

def domain_alignment(domain1, domain2, threshold):
    len1 = sum([(end - start) for start, end in domain1])
    len2 = sum([(end - start) for start, end in domain2])
    aligned = 0
    not_aligned = 0
    overlaps = []
    for start1, end1 in domain1:
        for start2, end2 in domain2:
            overlaps.append(calc_overlap(start1, end1, start2, end2))
    overlap = sum(overlaps)
    if overlap/len1 > threshold and overlap/len2 > threshold:
        return True
    else : return False

def update_corr_tb(corr_table, domains_dict, sword_domain, classification, treshold):
    found = False
    for id, dom in domains_dict.items():
        if domain_alignment(sword_domain, dom, treshold):
            corr_table[classification].append(id)
            corr_table[f'{classification}_renum_delim'].append(dom)
            return corr_table
    corr_table[classification].append('-')
    corr_table[f'{classification}_renum_delim'].append('-')
    return corr_table

def domain_corr(info_table, treshold, output):
    table = pd.read_csv(info_table, sep = '\t', usecols=["PDB", "SWORD_renum_delineation",  'ECOD_domain_ID', "ECOD_renum_delineation", 'CATH_ID', "CATH_renum_delineation", "SCOP_renum_delineation",  'SCOP_class', 'SCOP_fold', 'SCOP_supfamily', 'SCOP_family'])
    correspondence_table = {'SWORD2':[], 'SWORD2_renum_delim':[], 'ECOD':[], 'ECOD_renum_delim':[], 'SCOP':[], 'SCOP_renum_delim':[], 'CATH':[], 'CATH_renum_delim':[]}
    for prot in table.itertuples():
        SWORD2_domains = [[tuple(map(int, re.split(r'(?<=\d)-', items))) for items in item.split(',')] for item in prot.SWORD_renum_delineation.split(', ')]
        ECOD_domains = {domain_id: [tuple(map(int, re.split(r'(?<=\d)-', items))) for items in item.split(',')] for domain_id, item in zip(prot.ECOD_domain_ID.split(', '), prot.ECOD_renum_delineation.split(', '))}
        if prot.CATH_ID != '-':
            CATH_domains = {domain_id: [tuple(map(int, re.split(r'(?<=\d)-', items))) for items in item.split(',')] for domain_id, item in zip(prot.CATH_ID.split(', '), prot.CATH_renum_delineation.split(', '))}
            print(CATH_domains)
        if prot.SCOP_renum_delineation != '-':
            SCOP_domains = {f'{dom_class}.{fold}.{supfamily}.{family}': [tuple(map(int, re.split(r'(?<=\d)-', items))) for items in item.split(',')] for dom_class, fold, supfamily, family, item in zip(prot.SCOP_class.split(', '), prot.SCOP_fold.split(', '), prot.SCOP_supfamily.split(', '), prot.SCOP_family.split(', '), prot.SCOP_renum_delineation.split(', '))}
            print(SCOP_domains)
        for idx, sword_domain in enumerate(SWORD2_domains):
            domain_name = f"{prot.PDB}_SWORD2d{idx + 1}"
            correspondence_table['SWORD2'].append(domain_name)
            correspondence_table['SWORD2_renum_delim'].append(sword_domain)
            #ECOD
            correspondence_table = update_corr_tb(correspondence_table, ECOD_domains, sword_domain, 'ECOD', treshold)
            correspondence_table = update_corr_tb(correspondence_table, SCOP_domains, sword_domain, 'SCOP', treshold)
            correspondence_table = update_corr_tb(correspondence_table, CATH_domains, sword_domain, 'CATH', treshold)
            
    df = pd.DataFrame(correspondence_table)
    # get nbr of domains for correspondences
    ecod_strict = sum(df["SWORD2_renum_delim"] == df["ECOD_renum_delim"])
    scop_strict = sum(df["SWORD2_renum_delim"] == df["SCOP_renum_delim"])
    cath_strict = sum(df["SWORD2_renum_delim"] == df["CATH_renum_delim"])
    ecod_loose = sum(df["ECOD_renum_delim"] != '-')
    scop_loose = sum(df["SCOP_renum_delim"] != '-')
    cath_loose = sum(df["CATH_renum_delim"] != '-')
    print(df["CATH_renum_delim"] != '-')
    print(f'total number of domains: {len(df.index)}, loose overlap treshold : {treshold}')
    print(f'       ECOD   SCOP   CATH')
    print(f'strict {str(ecod_strict):<7s}{str(scop_strict):<7s}{str(cath_strict):<7s}')
    print(f'loose  {str(ecod_loose):<7s}{str(scop_loose):<7s}{str(cath_loose):<7s}')
    df.to_csv(output, sep = "\t", index = False)

def dict_prot_domain_pos(info_table, relate_ecod_scop_cath = False):
    table = pd.read_csv(info_table, sep = '\t', usecols=["PDB", "SWORD_renum_delineation", "ECOD_renum_delineation", "CATH_renum_delineation", "SCOP_renum_delineation"])
    table  = table.set_index("PDB")
    tb_dict = table.to_dict(orient = 'index')
    prot_domain_pos = {key :
    [[tuple(map(int, re.split(r'(?<=\d)-', items))) for items in item.split(',')] for item in sublist['SWORD_renum_delineation'].split(', ')]
    for key, sublist in tb_dict.items()}
    # if relate_ecod_scop_cath:
    #     prot_domain_pos_ECOD = {key :
    #     [[tuple(map(int, re.split(r'(?<=\d)-', items))) for items in item.split(',')] for item in sublist['ECOD_renum_delineation'].split(', ')]
    #     for key, sublist in tb_dict.items()}
    #     prot_domain_pos_SCOP = {key :
    #     [[tuple(map(int, re.split(r'(?<=\d)-', items))) for items in item.split(',')] for item in sublist['SCOP_renum_delineation'].split(', ')]
    #     for key, sublist in tb_dict.items()}
    #     prot_domain_pos_CATH = {key :
    #     [[tuple(map(int, re.split(r'(?<=\d)-', items))) for items in item.split(',')] for item in sublist['CATH_renum_delineation'].split(', ')]
    #     for key, sublist in tb_dict.items()}
    #     for prot, pos:
    #         domain_corr = domain_corr(prot_domain_pos[prot], prot_domain_pos_ECOD[prot], prot_domain_pos_SCOP[prot], prot_domain_pos_CATH[prot])

    return prot_domain_pos

def slice_pdb(pdb_name, input_pdb, ranges, output_dir):
    """
    Slices a PDB file into multiple PDB files based on given start-end positions.
    
    Args:
        input_pdb (str): The path to the input PDB file.
        ranges (list of tuples): List of tuples (start, end) specifying the atom range to slice.
        output_dir (str): Directory where the sliced PDB files will be saved.
    
    Returns:
        None
    """
    # Create a PDB parser object
    parser = PDB.PDBParser(QUIET=True)
    
    # Parse the input PDB file
    structure = parser.get_structure("structure", input_pdb)
    
    # Since there's only one chain, we can safely access the first chain
    model = structure[0]  # Usually there's only one model
    chain = model.get_list()[0]  # The first chain in the model
    
    # Loop through each range and create a sliced PDB file
    for idx, positions in enumerate(ranges):
        # Create a new structure containing only the atoms in the range
        new_structure = PDB.Structure.Structure("structure")
        new_model = PDB.Model.Model(0)
        new_structure.add(new_model)
        
        # Add a new chain to the model
        new_chain = PDB.Chain.Chain(chain.get_id())
        new_model.add(new_chain)
        print(idx, positions)
        # Add residues to the new chain within the specified range
        for (start, end) in positions:

            for residue in chain:
                residue_id = residue.get_id()[1]  # Residue number
                if start <= residue_id <= end:
                    new_chain.add(residue.copy())
        
        # Define the output PDB filename
        output_pdb = f"{output_dir}/{pdb_name}_SWORD2d{idx + 1}.pdb"
        
        # Write the sliced structure to a new PDB file
        io = PDB.PDBIO()
        io.set_structure(new_structure)
        io.save(output_pdb)
        
        print(f"Written sliced PDB file: {output_pdb}")


def process(atlas_info, pdb_folder, output):
    prot_domain_pos = dict_prot_domain_pos(atlas_info)
    pdb_files = list_pdb_files(pdb_folder)
    if not os.path.exists(output):
        os.mkdir(output)
    for pdb in pdb_files:
        pdb_name  = '_'.join(str(pdb).split('/')[-1].split('_')[0:2])
        slice_pdb(pdb_name, pdb, prot_domain_pos[pdb_name], output)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--programm", type=str, help="'cut', 'relate_to_ECOD_CATH_SCOP'")
    parser.add_argument("-d", "--pdb_folder", type=str, help="path to folder containing pdb files")
    parser.add_argument("-a", "--atlas_info", type=str)
    parser.add_argument("-o", "--output", type=str)
    parser.add_argument("-t", "--treshold_overlap", type=float)

    args = parser.parse_args()
    if args.programm == "cut":
        process(args.atlas_info, args.pdb_folder, args.output)
    if args.programm == "relate_to_ECOD_CATH_SCOP":
        domain_corr(args.atlas_info, args.treshold_overlap, args.output)