import Bio.PDB, sys, os
import numpy as np
from pathlib import Path

# usage: <pdb directory path> <contact matrix output path>
pdb_dir = sys.argv[1]
output_dir = sys.argv[2]

def calc_dist_matrix(structure, chain_id=None):
    if chain_id is not None:
        chains = [chain for chain in structure.get_chains() if chain.id == chain_id]
        if len(chains) != 1:
            print(f"{len(chains)} chains matching {chain_id}, 1 is required")
            exit()
        else:
            atoms = [a for a in chains[0].get_atoms() if a.get_name() == "CA"] # alpha carbons
    else:
        atoms = [a for a in structure.get_atoms() if a.get_name() == "CA"] # alpha carbons
    n_atoms = len(atoms)
    dist_matrix = np.zeros((n_atoms, n_atoms))
    for i in range(n_atoms):
        for j in range(i + 1, n_atoms):
            if abs(i-j) > 7:
                dist = np.inf # disregard closeness of contacts that are within 7 amino acids of each other in primary sequence
            else:
                dist = atoms[i] - atoms[j]
            dist_matrix[i, j] = dist
            dist_matrix[j, i] = dist
    return dist_matrix

def contact_matrix(dist_matrix, threshold=8): #8 angstrom threshold
    return dist_matrix < threshold

pdb_parser = Bio.PDB.PDBParser()


pdb_files = []
for root,dirs,files in os.walk(pdb_dir):
    for name in files:
        if name.endswith(('pdb')): 
            pdb_files.append(os.path.join(root,name))

Path(output_dir).mkdir(parents=True, exist_ok=True)
for pdb_file in pdb_files:
    id = str(Path(os.path.basename(pdb_file)).stem).split("_")[0]
    structure = pdb_parser.get_structure(id, pdb_file)
    # Calculate distance matrix
    dist_matrix = calc_dist_matrix(structure)
    # Create contact matrix
    result_matrix = contact_matrix(dist_matrix).astype(int)
    matrix_path = os.path.join(output_dir,str(Path(os.path.basename(pdb_file)).stem) + ".txt") 
    print(matrix_path)
    np.savetxt(matrix_path, result_matrix, header=id, fmt = '%d')
