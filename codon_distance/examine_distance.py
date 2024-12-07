import numpy as np
import glob
import pandas as pd
import os
import scipy.stats as st
import matplotlib.pyplot as plt
import seaborn as sns


rate_matrix_file = "../data/learned_rate_matrix_complexI_non_membrane.txt"

Q = np.zeros((20, 20))

amino_acids = ["A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V"]

with open(rate_matrix_file) as f:
    lines = f.readlines()
        

    for i in range(1,len(lines)):
        line = lines[i].split()[1:]
        for j in range(20):
            Q[i-1][j] = float(line[j])

# Create a heatmap
plt.figure(figsize=(10, 8))  # Adjust figure size for better readability
sns.heatmap(Q, annot=False, cmap="viridis", xticklabels=amino_acids, yticklabels=amino_acids)  # Change annot to True to display values
plt.title("Rate matrix for non membrane proteins")
plt.show()



amino_acid_to_codons = {
    'A': {'GCU', 'GCC', 'GCA', 'GCG'},  # Alanine
    'R': {'CGU', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'},  # Arginine
    'N': {'AAU', 'AAC'},  # Asparagine
    'D': {'GAU', 'GAC'},  # Aspartic Acid
    'C': {'UGU', 'UGC'},  # Cysteine
    'Q': {'CAA', 'CAG'},  # Glutamine
    'E': {'GAA', 'GAG'},  # Glutamic Acid
    'G': {'GGU', 'GGC', 'GGA', 'GGG'},  # Glycine
    'H': {'CAU', 'CAC'},  # Histidine
    'I': {'AUU', 'AUC', 'AUA'},  # Isoleucine
    'L': {'UUA', 'UUG', 'CUU', 'CUC', 'CUA', 'CUG'},  # Leucine
    'K': {'AAA', 'AAG'},  # Lysine
    'M': {'AUG'},  # Methionine (start codon)
    'F': {'UUU', 'UUC'},  # Phenylalanine
    'P': {'CCU', 'CCC', 'CCA', 'CCG'},  # Proline
    'S': {'UCU', 'UCC', 'UCA', 'UCG', 'AGU', 'AGC'},  # Serine
    'T': {'ACU', 'ACC', 'ACA', 'ACG'},  # Threonine
    'W': {'UGG'},  # Tryptophan
    'Y': {'UAU', 'UAC'},  # Tyrosine
    'V': {'GUU', 'GUC', 'GUA', 'GUG'},  # Valine
}





def average_hamming_distance(set1, set2):
    def hamming_distance(s1, s2):
        """Compute the Hamming distance between two strings of equal length."""
        if len(s1) != len(s2):
            raise ValueError("Strings must be of the same length to calculate Hamming distance.")
        return sum(c1 != c2 for c1, c2 in zip(s1, s2))
    
    total_distance = 0
    pair_count = 0

    for s1 in set1:
        for s2 in set2:
            if len(s1) == len(s2):  # Only compare strings of the same length
                total_distance += hamming_distance(s1, s2)
                pair_count += 1

    return total_distance / pair_count if pair_count > 0 else 0.0





#Compute hamming distance matrix
hamming_matrix = np.zeros((20, 20))

for i, aa1 in enumerate(amino_acids):
    for j, aa2 in enumerate(amino_acids):
        hamming_matrix[i][j] = average_hamming_distance(amino_acid_to_codons[aa1], amino_acid_to_codons[aa2])




# Create a heatmap
plt.figure(figsize=(10, 8))  # Adjust figure size for better readability
sns.heatmap(hamming_matrix, annot=False, cmap="viridis", xticklabels=amino_acids, yticklabels=amino_acids)  # Change annot to True to display values
plt.title("Average Hamming distance matrix")
plt.show()






rate_list = []
hamming_list = []

for i in range(20):
    for j in range(20):
        if(i != j):
            rate_list.append(Q[i][j])
            hamming_list.append(hamming_matrix[i][j])


sns.set_theme(style="whitegrid", context="talk")
plt.figure(figsize=(10, 8))
sns.regplot(x=hamming_list,y=rate_list, scatter=True, line_kws={"color": "darkblue", "linewidth": 2}, 
            scatter_kws={"s": 100, "color": "skyblue"}, ci=95)
plt.title("Non membrane protein transition probability to Hamming distance comparison")
plt.xlabel("Average Hamming distance between codons")
plt.ylabel("Transition Probability")
plt.show()






"""
#Treating it as nucleotide transitions
nucleotide_transitions = ["AC","AG","AU","CG","CU","GA","GC","GU","UA","UC","UG"]






def transition_count(set1, set2, transition):
    total_transition_count = 0
    pair_count = 0

    for s1 in set1:
        for s2 in set2:
            if len(s1) == len(s2):  # Only compare strings of the same length
                for i in range(len(s1)):
                    if s1[i] == transition[0] and s2[i] == transition[1]:
                        total_transition_count += 1
                    pair_count += 1

    return total_transition_count / pair_count if pair_count > 0 else 0.0


amino_acid_mutations = []
for i in range(len(amino_acids)):
    for j in range(len(amino_acids)):
        if i != j:
            amino_acid_mutations.append(amino_acids[i]+amino_acids[j])


#This matrix maps amino acid transitions to the proportional nucleotide transitions that could cause it
amino_muts_to_nucl_muts = np.zeros((len(amino_acid_mutations), len(nucleotide_transitions)))

for i in range(len(amino_acid_mutations)):
    for j in range(len(nucleotide_transitions)):
        aa1 = amino_acid_mutations[i][0]
        aa2 = amino_acid_mutations[i][1]

        amino_muts_to_nucl_muts[i][j] = transition_count(amino_acid_to_codons[aa1],amino_acid_to_codons[aa2], nucleotide_transitions[j])

# Create a heatmap
sns.heatmap(amino_muts_to_nucl_muts, annot=False, cmap="viridis", xticklabels=nucleotide_transitions, yticklabels=amino_acid_mutations)  # Change annot to True to display values
plt.title("Hamming matrix")
plt.show()
"""