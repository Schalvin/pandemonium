import pandas as pd
import re
import os

from pathlib import Path
from Bio import PDB


PDB_FOLDER = "/home/chili/chalvin/projet_long/ATLAS/prot_pdb"
ATLAS_TB = "/home/chili/chalvin/projet_long/ATLAS/2023_03_09_ATLAS_info.tsv"
OUTPUT_DIR = "/home/chili/chalvin/projet_long/ATLAS/SWORDdomains_pdb"

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

def dict_prot_domain_pos(info_table):
    table = pd.read_csv(info_table, sep = '\t', usecols=["PDB", "SWORD_renum_delineation", "ECOD_renum_delineation", "CATH_renum_delineation", "SCOP_renum_delineation"])
    table  = table.set_index("PDB")
    prot_domain_pos = table.to_dict(orient = 'index')
    prot_domain_pos = {key :
    [[tuple(map(int, re.split(r'(?<=\d)-', items))) for items in item.split(',')] for item in sublist['SWORD_renum_delineation'].split(', ')]
    for key, sublist in prot_domain_pos.items()}
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


def process():
    prot_domain_pos = dict_prot_domain_pos(ATLAS_TB)
    pdb_files = list_pdb_files(PDB_FOLDER)
    if not os.path.exists(OUTPUT_DIR):
        os.mkdir(OUTPUT_DIR)
    for pdb in pdb_files:
        pdb_name  = '_'.join(str(pdb).split('/')[-1].split('_')[0:2])
        slice_pdb(pdb_name, pdb, prot_domain_pos[pdb_name], OUTPUT_DIR)


if __name__ == '__main__':
    process()