import glob
import pandas as pd

import os
def main():
    paths = glob.glob('/Users/mikemoubarak/cherryproj/data/msa_ETC_fasta_files/*.txt')

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
    main()