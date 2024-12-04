import numpy as np
import matplotlib.pyplot as plt




def calc_stationary_dist_from_rate_matrix_file(rate_matrix_file):
    Q = np.zeros((20, 20))

    with open(rate_matrix_file) as f:
        lines = f.readlines()
        

        for i in range(1,len(lines)):
            line = lines[i].split()[1:]
            for j in range(20):
                Q[i-1][j] = float(line[j])
    

    Q_T = Q.T

    # Add the normalization condition
    A = np.vstack([Q_T[:-1], np.ones(Q.shape[1])])  # Replace one equation with sum(pi) = 1
    b = np.zeros(Q.shape[0])
    b[-1] = 1  # Corresponds to sum(pi) = 1

    # Solve the system
    pi = np.linalg.lstsq(A, b, rcond=None)[0]

    return pi









matrix_file = "../data/learned_rate_matrix_complexI.txt"

pi = calc_stationary_dist_from_rate_matrix_file(matrix_file)
print(pi)



# States corresponding to the stationary probabilities
states = ["A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V"]

# Plot the stationary distribution
plt.figure(figsize=(8, 5))
plt.bar(states, pi, color='skyblue', edgecolor='black')

# Add labels and title
plt.xlabel('States', fontsize=14)
plt.ylabel('Probability', fontsize=14)
plt.title('Stationary Distribution', fontsize=16)
plt.ylim(0, 1)  # Set y-axis range to [0, 1]
plt.grid(axis='y', linestyle='--', alpha=0.7)

# Display values on top of bars
for i, prob in enumerate(pi):
    plt.text(i, prob + 0.02, f"{prob:.2f}", ha='center', fontsize=12)

plt.tight_layout()
plt.show()

