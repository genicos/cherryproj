import numpy as np





def calc_stationary_dist_from_rate_matrix_file(rate_matrix_file):
    matrix = np.zeros((20, 20))

    with open(rate_matrix_file) as f:
        lines = f.readlines()
        

        for i in range(1,len(lines)):
            line = lines[i].split()[1:]
            for j in range(20):
                matrix[i-1][j] = float(line[j])
    
    print(matrix)






Kagan_matrix_file = "../data/learned_rate_matrix_complexI.txt"

calc_stationary_dist_from_rate_matrix_file(Kagan_matrix_file)