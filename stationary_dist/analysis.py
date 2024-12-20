import numpy as np
import glob
import pandas as pd
import os
import scipy.stats as st
import matplotlib.pyplot as plt




def calc_stationary_dist_from_rate_matrix_file(rate_matrix_file):
    Q = np.zeros((20, 20))

    with open(rate_matrix_file) as f:
        lines = f.readlines()
        

        for i in range(1,len(lines)):
            line = lines[i].split()[1:]
            for j in range(20):
                Q[i-1][j] = float(line[j])
    
    x_labels = ["A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V"]
    y_labels = ["A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V"]
    Q_T = Q.T
    
    fig, axes = plt.subplots(1,1)
    im = axes.imshow(Q.T, cmap = 'Blues')
    cbar = fig.colorbar(im)
    cbar.set_label('Generator Value')
    axes.set_xticks(np.arange(len(x_labels)))  # Position of ticks
    axes.set_xticklabels(x_labels)  # Labels for ticks
    axes.set_yticks(np.arange(len(y_labels)))
    axes.set_yticklabels(y_labels)
    axes.set_xlabel('Amino Acid')
    axes.set_ylabel('Amino Acid')
    fig.savefig('../figures/figure1.tiff', dpi = 300)
    plt.close()

    # Add the normalization condition
    A = np.vstack([Q_T[:-1], np.ones(Q.shape[1])])  # Replace one equation with sum(pi) = 1
    b = np.zeros(Q.shape[0])
    b[-1] = 1  # Corresponds to sum(pi) = 1

    # Solve the system
    pi = np.linalg.lstsq(A, b, rcond=None)[0]

    return pi

def getCurrentDistributions():
    paths = glob.glob('../data/msa_ETC_fasta_files/*.txt')

    outputDataFrame = pd.DataFrame()
    familyNames = []
    for path in paths:
        with open(path, "r") as file:
        # Read the contents of the file
            content = file.readlines()
    
            familyNames.append(os.path.split(path)[1][:-4])
    
            countDictionary = getDistributionPerFamily(content)
    
            newRow = pd.DataFrame(countDictionary) 
            newRow /= newRow.iloc[0,:].sum()
           
         
            outputDataFrame = pd.concat([outputDataFrame, newRow], axis = 0)
        
    outputDataFrame.index = familyNames
    return outputDataFrame

def getDistributionPerFamily(list):

    aminoAcidDictionary = {'A':[0],'R':[0],'N':[0],'D':[0],'C':[0],'Q':[0],
                           'K':[0],'L':[0],'I':[0],'H':[0],'G':[0],'E':[0], 
                           'M':[0],'F':[0],'P':[0],'S':[0],'T':[0],'W':[0],
                           'Y':[0],'V':[0]}
    for element in list:
        if element[0] != '>':
            for character in element[:-1]:
                if character != '-':
                    if character in aminoAcidDictionary.keys():
                        
                        aminoAcidDictionary[character][0] += 1
    return aminoAcidDictionary

if __name__ == '__main__':


    Kagan_matrix_file = "../data/learned_rate_matrix_complexI.txt"

    pi = calc_stationary_dist_from_rate_matrix_file(Kagan_matrix_file)
    
    piDF = pd.DataFrame(pi)
    piDF.columns = ['value']
    piDF.index = ["A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V"]
    
    piDF.sort_values(by = ['value'], axis = 0, inplace = True, ascending= False)
    
    
    
    fig, axes = plt.subplots(1,1)
    axes.bar(piDF.index, piDF['value'], color = 'gray')
    axes.set_xlabel('Amino Acid')
    axes.set_ylabel('Steady State Frequency')
    fig.savefig('../figures/figure2.tiff', dpi = 300)
    plt.close()

    pi = calc_stationary_dist_from_rate_matrix_file(Kagan_matrix_file)
    pi2 = pd.DataFrame(pi)
    pi2.index = ["A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V"]

    currentDistributions = getCurrentDistributions()

    KL_divergences = {}
    for row in range(currentDistributions.shape[0]):
        rowValues = np.array(currentDistributions.iloc[row,:])

        
        KL_divergences[str(currentDistributions.index[row])] = [st.entropy(rowValues, pi)]
        
    klDF = pd.DataFrame(KL_divergences).T
    klDF.columns = ['KL']
    klDF.sort_values(by = 'KL', axis = 0, inplace = True, ascending = False)
    
    fig, axes = plt.subplots(1,1, figsize = (10,8))
    axes.bar(klDF.index, klDF['KL'], color = 'gray')
    axes.set_xlabel('Protein Family')
    axes.set_ylabel('KL Divergence')
    axes.tick_params(axis='x', rotation=75, size = 2)
    fig.savefig('../figures/figure3.tiff', dpi = 300)
    plt.close()

    x_labels = ["A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V"]
    y_labels = ["A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V"]
    fig, axes = plt.subplots(1,1)

    highCorr = st.pearsonr(pi, np.array(currentDistributions.loc['PF01059'])).statistic 
    lowCorr = st.pearsonr(pi, np.array(currentDistributions.loc['PF01370'])).statistic 
    
    axes.scatter(pi, currentDistributions.loc['PF01059'], s = 5, c = 'r', label = 'Highest (PF01059) KL - corr = '+ str(round(highCorr, 2)))
    axes.scatter( pi, currentDistributions.loc['PF01370'], s = 5, c = 'g', label = 'Lowest (PF01370) KL - corr = '+str(round(lowCorr, 2)))
    axes.plot(np.linspace(0, max(currentDistributions.loc['PF01059']), num = 100), 
    np.linspace(0, max(currentDistributions.loc['PF01059']), num = 100), c = 'k', linewidth = 1)
    axes.set_ylabel('Observed Distribution')
    axes.set_xlabel('Steady State Distribution')
    axes.legend(bbox_to_anchor=(0, 1.15), loc='upper left')
    #axes.set_xticks(np.arange(len(x_labels)))  # Position of ticks
    #axes.set_xticklabels(x_labels)  # Labels for ticks
    #axes.set_yticks(np.arange(len(y_labels)))
    #axes.set_yticklabels(y_labels)

    fig.savefig('../figures/figure4.tiff', dpi = 300)
#matrix_file = "../data/learned_rate_matrix_complexI.txt"

#pi = calc_stationary_dist_from_rate_matrix_file(matrix_file)
#print(pi)



# States corresponding to the stationary probabilities
#states = ["A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V"]

# Plot the stationary distribution
#plt.figure(figsize=(8, 5))
#plt.bar(states, pi, color='skyblue', edgecolor='black')

# Add labels and title
#plt.xlabel('States', fontsize=14)
#plt.ylabel('Probability', fontsize=14)
#plt.title('Stationary Distribution', fontsize=16)
#plt.ylim(0, 1)  # Set y-axis range to [0, 1]
#plt.grid(axis='y', linestyle='--', alpha=0.7)

# Display values on top of bars
#for i, prob in enumerate(pi):
#    plt.text(i, prob + 0.02, f"{prob:.2f}", ha='center', fontsize=12)

#plt.tight_layout()
#plt.show()

