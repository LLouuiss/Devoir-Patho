"""Modélisation de la diffusion 2D de chlore dans une section de béton avec schéma implicite"""

# Importation
import matplotlib.pyplot as plt
import scipy.linalg
import sys
import numpy as np
import time

# Paramètres
theta = 0.5 # Choix du theta
X = 100 # Longueur du béton selon x (mm)
Y = 100 # Longueur du béton selon y (mm)
T = 100 # Durée totale d'analyse (années)
D_app = 0.000001 # Coefficient de diffusion supposé constant (mm^2/s)
C_s = 1 # Concentration à la surface (%)

# Partie numérique : différences finies (implicite)

# Discrétisation
Nx = 100 # Nombre de points spatiaux selon x
Ny = 100 # Nombre de points spatiaux selon y
Deltax = X/(Nx) # Pas de discrétisation selon x
Deltay = Y/(Ny) # Pas de discrétisation selon y
Deltat = 1 # Incrément temporel
Nt = int((T/Deltat) + 1) # Nombre de points temporels

# Listes
x = np.linspace(0, X, Nx) # Liste des valeurs de profondeur x (mm)
y = np.linspace(0, Y, Ny) # Liste des valeurs de profondeur y (mm)
t = np.linspace(0, T, Nt) # Liste des valeurs temporelles (années)
c = np.zeros((Nx+1, Ny)) # Liste des valeurs de concentrations initiales C(x) initié à 0
c_n = np.zeros((Nx+1, Ny)) # Liste des valeurs de concentrations finales C(x) initié à 0

t0 = time.time() # Temps initial processus

# Mapping
m = lambda i, j: j*(Nx) + i

# Creation de la matrice A et du vecteur b
A = np.zeros(((Nx)*(Ny), (Nx) *(Ny)))
b = np.zeros((Nx)*(Ny))
Fx = D_app*Deltat*365*24*3600/Deltax**2
Fy = D_app*Deltat*365*24*3600/Deltay**2
for j in range(Ny):
    for i in range(Nx):
        p = m(i, j)
        if j == 0 or j == Ny-1 or i == 0 or i == Nx-1:
            A[p, p] = 1 # Conditions frontières sur les 4 bords
        else:
            A[p, m(i, j - 1)] = -theta * Fy
            A[p, m(i - 1, j)] = -theta * Fx
            A[p, p] = 1 + 2 * theta * (Fx + Fy)
            A[p, m(i + 1, j)] = -theta * Fx
            A[p, m(i, j + 1)] = -theta * Fy

# Mise à jour de la liste concentration C[x, y, t]
C_num = []
for n in range(0, Nt):
    # Iteration dans b
    for j in range(Ny):
        for i in range(Nx):
            p = m(i, j)
            # Conditions de bord
            if j == 0 or j == Ny-1 or i == 0 or i == Nx-1:
                b[p] = C_s # Condition en C[0, y, t], C[1000, y, t], C[x, 0, t], C[x, 1000, t]
            else:
                b[p] = c_n[i, j] + (1 - theta) * (
                    Fx * (c_n[i + 1, j] - 2 * c_n[i, j] + c_n[i - 1, j]) +
                    Fy * (c_n[i, j + 1] - 2 * c_n[i, j] + c_n[i, j - 1])
                )

    # Resolution du systeme
    c_sys = np.linalg.solve(A, b)
    # Remplissage de c avec c_sys
    for i in range(Nx):
        for j in range(Ny):
            c[i,j] = c_sys[m(i,j)]
    C_num.append(c.copy())
    c_n = c
C_num = np.array(C_num)
t1 = time.time() # Temps final processus

plt.imshow(C_num[-1].T, extent=[0, X, 0, Y], origin='lower', cmap='Blues')
plt.colorbar(label='Concentration (%)')
plt.title("Concentration de chlore après {} années".format(T))
plt.xlabel('x (mm)')
plt.ylabel('y (mm)')
plt.show()
