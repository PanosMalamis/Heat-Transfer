import numpy as np
import matplotlib.pyplot as plt

nx = 20
ny = 20
lx = 1
ly = 1
dx = lx/nx
dy = lx/ny

# case 1

T_right = 50
T_bot = 20
T_inf = 30

k = 60
h = 0.6

T0 = np.ones(nx*ny, dtype=np.float64)*20
B = np.zeros_like(T0)
A = np.zeros((nx*ny, nx*ny), dtype=np.float64)

for i in range(nx):
    for j in range(ny):
        # print (i,j)
        index = i + j*ny

        if (i >= 1) and (i <= nx-2) and (j >= 1) and (j <= ny-2):
            # print ('interior node')
            A[index, index] = -4  # (i, j)
            A[index, index+1] = 1  # (i+1, j)
            A[index, index-1] = 1  # (i-1, j)
            A[index, index+ny] = 1  # (i, j+1)
            A[index, index-ny] = 1  # (i, j-1)

        elif (i >= 1) and (i <= nx-2) and (j == ny-1):
            # print ('top node') AAA
            A[index, index] = -4  # (i, j)
            A[index, index+1] = 1  # (i+1, j)
            A[index, index-ny] = 2  # (i, j-1)
            A[index, index-1] = 1  # (i-1, j)

           # A[index, index] = 1 # (i, j)
           # B[index] = T_top # (i, j)

        elif (i >= 1) and (i <= nx-2) and (j == 0):
            # print ('bottom node')
            A[index, index] = 1  # (i, j)
            B[index] = T_bot

        elif (i == 0) and (j >= 1) and (j <= ny-2):
            # print ('left node') AAAAA
            A[index, index] = -4  # (i, j)
            A[index, index+1] = 2  # (i+1, j)
            A[index, index+ny] = 1  # (i, j+1)
            A[index, index-ny] = 1  # (i, j-1)

        elif (i == nx-1) and (j >= 1) and (j <= ny-2):
            # print ('right node')

            A[index, index] = 1  # (i, j)
            B[index] = T_right  # (i, j)

            # A[index, index] = -2 * (h * dx/k + 2) # (i, j)
            # A[index, index-1] = 2 # (i-1, j)
            # A[index, index+ny] = 1 # (i, j+1)
            # A[index, index-ny] = 1 # (i, j-1)
            # B[index] = -2 * h * dx/k * T_inf # (i, j)

        elif (i == 0) and (j == 0):
            # print ('bottom left')

            A[index, index] = -2  # (i, j)
            A[index, index+1] = 1  # (i+1, j)
            A[index, index+ny] = 1  # (i, j+1)

        elif (i == nx-1) and (j == 0):
            # print ('bottom right')

            A[index, index] = -2  # (i, j)
            A[index, index-1] = 1  # (i-1, j)
            A[index, index+ny] = 1  # (i, j+1)

        elif (i == 0) and (j == ny-1):
            # print ('top left')

            A[index, index] = -2  # (i, j)
            A[index, index+1] = 1  # (i+1, j)
            A[index, index-ny] = 1  # (i, j-1)

        elif (i == nx-1) and (j == ny-1):
            # print ('top right')

            A[index, index] = -2  # (i, j)
            A[index, index-1] = 1  # (i-1, j)
            A[index, index-ny] = 1  # (i, j-1)


def gauss_seidel(A, B, x0, epsilon, max_iterations):
    n = len(A)
    x = x0.copy()

    # gauss-seidal method (by bottom science)

    for i in range(max_iterations):
        x_new = np.zeros(n)
        for j in range(n):
            s1 = np.dot(A[j, :j], x_new[:j])
            s2 = np.dot(A[j, j+1:], x[j+1:])
            x_new[j] = (B[j] - s1 - s2) / A[j, j]
        if np.allclose(x, x_new, rtol=epsilon):
            return x_new
        x = x_new
    return x


def unwrap(theta, nx, ny):
    T = np.zeros((nx, ny), dtype=np.float64)
    for i in range(nx):
        for j in range(ny):
            T[i, j] = theta[i + j*ny]
    return np.transpose(T)


x = gauss_seidel(A, B, T0, 0.0001, 1000)

T = unwrap(x, nx, ny)
fig, ax = plt.subplots()
CS = ax.contour(T)
ax.clabel(CS, inline=True, fontsize=10)
plt.show()
