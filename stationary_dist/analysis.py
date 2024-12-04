import numpy as np





def calc_stationary_dist_from_rate_matrix_file(rate_matrix_file):
    Q = np.zeros((20, 20))

    with open(rate_matrix_file) as f:
        lines = f.readlines()
        

        for i in range(1,len(lines)):
            line = lines[i].split()[1:]
            for j in range(20):
                Q[i-1][j] = float(line[j])
    
    print(Q)

    Q_T = Q.T

    # Add the normalization condition
    A = np.vstack([Q_T[:-1], np.ones(Q.shape[1])])  # Replace one equation with sum(pi) = 1
    b = np.zeros(Q.shape[0])
    b[-1] = 1  # Corresponds to sum(pi) = 1

    # Solve the system
    pi = np.linalg.lstsq(A, b, rcond=None)[0]

    print("Stationary distribution:", pi)






Kagan_matrix_file = "../data/learned_rate_matrix_complexI.txt"

calc_stationary_dist_from_rate_matrix_file(Kagan_matrix_file)