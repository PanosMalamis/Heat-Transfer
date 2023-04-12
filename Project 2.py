# Panos Malamis
# 4/6/2023
# ME 321 - Heat Transfer


'''Assign nodes (FROM PROFESSOR)'''

import numpy as np
for m in range(nx):  # Loop over the nodes in x direction
    for n in range(ny):  # Loop over the nodes in the y direction
        # print(m,n)
        index = i + j*ny  # This is the rule for how a node gets labeled in A, T, or C matrix

        if (m >= 1) and (m <= nx-2) and (n >= 1) and (n <= ny-2):  # Interior node rule…
            # print('This is an interior node')
            # Follow Rule 1 for Table 4.2 (Eq. 4.29)
            A[index, index] = -4  # Set the A[m,n] to -4
            A[index, index+1] = 1  # Set the A[m+1,n] to 1
            A[index, index-1] = 1  # Set the A[m-1,n] to 1
            A[index, index-ny] = 1  # Set the A[m,n-1] to 1
            A[index, index+ny] = 1  # Set the A[m,n+1] to 1
# Note that there are no terms on the RHS for Rule 1, so you don’t have to change anything in [C].

'''Gauss-Seidel Solver'''


def gauss_seidel(A, b, x0, epsilon, max_iterations):
    n = len(A)
    x = x0.copy()
 # Gauss-Seidal Method [By Bottom Science]
 for i in range(max_iterations):
    x_new = np.zeros(n)
    for j in range(n):
        s1 = np.dot(A[j, :j], x_new[:j])
        s2 = np.dot(A[j, j + 1:], x[j + 1:])
        x_new[j] = (b[j] - s1 - s2) / A[j, j]
 if np.allclose(x, x_new, rtol=epsilon):
        return x_new
    x = x_new
        return x