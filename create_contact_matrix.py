import Bio.PDB, sys
import numpy as np

# usage: <PDB code> <pdb file path>"
pdb_code = sys.argv[1]
pdb_file = sys.argv[2]

def calc_dist_matrix(structure):
    atoms = [a for a in structure.get_atoms() if a.get_name() == "CA"] # carbon alpha atoms
    n_atoms = len(atoms)
    dist_matrix = np.zeros((n_atoms, n_atoms))
    for i in range(n_atoms):
        for j in range(i + 1, n_atoms):
            dist = atoms[i] - atoms[j]
            dist_matrix[i, j] = dist
            dist_matrix[j, i] = dist
    return dist_matrix

def contact_matrix(dist_matrix, threshold=12.0):
    return dist_matrix < threshold

pdb_parser = Bio.PDB.PDBParser()
structure = pdb_parser.get_structure(pdb_code, pdb_file)

# Calculate distance matrix
dist_matrix = calc_dist_matrix(structure)

# Create contact matrix
contact_matrix = contact_matrix(dist_matrix).astype(int)


print(contact_matrix)