# Panos Malamis
# Brandon Wong
# 04/09/2023

# HEAT TRANSFER PROJECT

import numpy as np
import matplotlib.pyplot as plt
import time as t

# ------------------------Functions ---------------------------#


def gauss_seidel(A, b, x0, epsilon, max_iterations):
    """This function will solve the system of linear equations using the Gauss Seidel method. From the Professors Rubric."""
    n = len(A)
    x = x0.copy()

    # gauss-seidal method (by bottom science)

    for i in range(max_iterations):
        x_new = np.zeros(n)
        for j in range(n):
            s1 = np.dot(A[j, :j], x_new[:j])
            s2 = np.dot(A[j, j+1:], x[j+1:])
            x_new[j] = (b[j] - s1 - s2) / A[j, j]
        if np.allclose(x, x_new, rtol=epsilon):
            return x_new
        x = x_new
    return x


def unwrap(theta, nx, ny):
    """ Transposes the matrix [T] """
    T = np.zeros((nx, ny), dtype=np.float64)
    for i in range(nx):
        for j in range(ny):
            T[i, j] = theta[i + j*ny]
    return np.transpose(T)


# -------- VARIABLES --------#
nx = 20
ny = 20
Lx = Ly = 1
dx = Lx/nx
dy = Ly/ny

# ---------- MAIN ----------#

# We made it easier for the grader to pull up each graph.
print("Hello, and welcome to our Heat Transfer Project!")
t.sleep(2)
print("This Program should be able to solve most problems...")
t.sleep(2)
print("but for your convenience, please select a case study.")
t.sleep(2)

user_input = input(
    "\nSelect either 1 , 2 , or 3.\nEnter '4' for your custom values\n \nEnter your choice: ")
# User input above switches and displays the tables that you ask for. CASE 1, 2, or 3.

if int(user_input) == 1:  # CASE 1

    T_top = 50  # side 2
    T_bot = 20  # side 4
    T_inf = 30  # side 3
# side 1 is adiabatic
    k = 60  # conduction constant
    h = 100  # convection constant

    T0 = np.ones(nx*ny, dtype=np.float64)*20
    b = np.zeros_like(T0)
    A = np.zeros((nx*ny, nx*ny), dtype=np.float64)
# Graph Customizations
    cases = "Case #1 - "


if int(user_input) == 2:  # CASE 2

    T_top = 50
    T_bot = 20
    T_inf = 30

    k = 60
    h = 0.6

    T0 = np.ones(nx*ny, dtype=np.float64)*20
    b = np.zeros_like(T0)
    A = np.zeros((nx*ny, nx*ny), dtype=np.float64)

    cases = "Case #2 - "


if int(user_input) == 3:  # CASE 3

    T_right = 50
    T_bot = 20

# for case 3, we remove h, k, and T-Infinity. They will no longer be used.
    k = 60
    h = 1

    T0 = np.ones(nx*ny, dtype=np.float64)*20
    b = np.zeros_like(T0)
    A = np.zeros((nx*ny, nx*ny), dtype=np.float64)

    cases = "Case #3 - "

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

            if (i >= 1) and (i <= nx-2) and (j == ny-1):
                # print ('top node') AAA
                A[index, index] = -4  # (i, j)
                A[index, index+1] = 1  # (i+1, j)
                A[index, index-ny] = 2  # (i, j-1)
                A[index, index-1] = 1  # (i-1, j)

            # A[index, index] = 1 # (i, j)
            # B[index] = T_top # (i, j)

            if (i >= 1) and (i <= nx-2) and (j == 0):
                # print ('bottom node')
                A[index, index] = 1  # (i, j)
                b[index] = T_bot

            if (i == 0) and (j >= 1) and (j <= ny-2):
                # print ('left node') AAAAA
                A[index, index] = -4  # (i, j)
                A[index, index+1] = 2  # (i+1, j)
                A[index, index+ny] = 1  # (i, j+1)
                A[index, index-ny] = 1  # (i, j-1)

            if (i == nx-1) and (j >= 1) and (j <= ny-2):
                # print ('right node')

                A[index, index] = 1  # (i, j)
                b[index] = T_right  # (i, j)

            if (i == 0) and (j == 0):
                # print ('bottom left')

                A[index, index] = -2  # (i, j)
                A[index, index+1] = 1  # (i+1, j)
                A[index, index+ny] = 1  # (i, j+1)

            if (i == nx-1) and (j == 0):
                # print ('bottom right')

                A[index, index] = -2  # (i, j)
                A[index, index-1] = 1  # (i-1, j)
                A[index, index+ny] = 1  # (i, j+1)

            if (i == 0) and (j == ny-1):
                # print ('top left')

                A[index, index] = -2  # (i, j)
                A[index, index+1] = 1  # (i+1, j)
                A[index, index-ny] = 1  # (i, j-1)

            if (i == nx-1) and (j == ny-1):
                # print ('top right')

                A[index, index] = -2  # (i, j)
                A[index, index-1] = 1  # (i-1, j)
                A[index, index-ny] = 1  # (i, j-1)

    x = gauss_seidel(A, b, T0, 0.0001, 1000)

    T = unwrap(x, nx, ny)
    fig, ax = plt.subplots()
    CS = ax.contour(T)
    ax.clabel(CS, inline=True, fontsize=10)

    print("\n   Loading Graph now, Please wait...")
    t.sleep(2)
    ax.set_title(cases + "Temperature Gradient")
    ax.set_xlabel("Nodes")
    ax.set_ylabel("Nodes")
    plt.show()

    quit()


# --------------CHANGE THE CODE BELOW FOR CUSTOM RESULTS ----------------------#
if int(user_input) == 4:

    T_top = 50  # side 2
    T_bot = 20  # side 4
    T_inf = 30  # side 3

    k = 60  # conduction constant
    h = 100  # convection constant

    T0 = np.ones(nx*ny, dtype=np.float64)*20
    b = np.zeros_like(T0)
    A = np.zeros((nx*ny, nx*ny), dtype=np.float64)

    cases = "Custom Case -"


# ------------------------------------------------------------------------------#

"""From Professors Rubric"""

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

        if (j == ny-1):
            # print ('top node')
            A[index, index] = 1  # (i, j)
            b[index] = T_top  # (i, j)

        if (j == 0):
            # print ('bottom node')
            A[index, index] = 1  # (i, j)
            b[index] = T_bot

        if (i == 0) and (j >= 1) and (j <= ny-2):
            # print ('left node')
            A[index, index] = -4  # (i, j)
            A[index, index+1] = 2  # (i+1, j)
            A[index, index+ny] = 1  # (i, j+1)
            A[index, index-ny] = 1  # (i, j-1)

        if (i == nx-1) and (j >= 1) and (j <= ny-2):
            # print ('right node')
            A[index, index] = -2 * (h * dx/k + 2)  # (i, j)
            A[index, index-1] = 2  # (i-1, j)
            A[index, index+ny] = 1  # (i, j+1)
            A[index, index-ny] = 1  # (i, j-1)
            b[index] = -2 * h * dx/k * T_inf  # (i, j)

x = gauss_seidel(A, b, T0, 0.0001, 1000)

T = unwrap(x, nx, ny)
fig, ax = plt.subplots()
CS = ax.contour(T)
ax.clabel(CS, inline=True, fontsize=10)

print("\n   Loading Graph now, Please wait...")
t.sleep(2)
ax.set_title(cases + " Temperature Gradient")
ax.set_xlabel("Nodes")
ax.set_ylabel("Nodes")
plt.show()
