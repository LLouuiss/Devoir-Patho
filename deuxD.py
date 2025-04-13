""" Implémentation des automates cellulaires de la diffusion 2D
    
Le code a été repris et adapté à notre géométrie à partir du mémoire de fin
d’étude "Méthode d’évaluation numérique de la diffusion du chlore à l’intérieur du béton armé" de Arthur Herbiet.

"""
#Import
import numpy as np
import math
import matplotlib.pyplot as plt
import time


def calcul_2D() :
    # Paramètres
    X = 1400  # Longueur de béton selon x (mm)
    Y = 2000  # Longueur de béton selon y (mm)
    T = 150  # Temps total (année)
    D_app = 5 * (10**(-7))  # Coefficient de diffusion constant (mm^2/s)
    C_s = 0.15  # Concentration de chlorure en surface (%)

    # Coefficients d'évolution
    phi_0 = 1/2
    phi_i = (1 - phi_0) / 4

    # Discrétisation
    Nx = 140  # Nombre de points spatiaux selon x
    Ny = 200  # Nombre de points spatiaux selon y
    Deltax = X / Nx  # Pas de discrétisation selon x
    Deltay = Y / Ny  # Pas de discrétisation selon y
    Deltat = (1 / (D_app * 4)) * ((1 - phi_0) * Deltax**2) * (1 / (365 * 24 * 3600))  # Critère de régulation
    #Deltat = 0.39  # Incrément temporel choisi
    Nt = int((T / Deltat) + 1)  # Nombre de points temporels
    C_num = np.zeros((Nx, Ny, Nt))  # Matrice des concentrations C(x,t) initié à 0

    # Listes
    x = np.linspace(0, X, Nx)  # Liste des valeurs de profondeur x (mm)
    y = np.linspace(0, Y, Ny)  # Liste des valeurs de pronfondeur y (mm)
    t = np.linspace(0, T, Nt)  # Liste des valeurs de temps (année)

    # Condition de bord

    C_num[:, -1, :] = C_s  # Condition en C[x, 2000, t]
    C_num[:, 0, :] = C_s  # Condition en C[x, 0, t]

    # Poutr en en I semelle de 30 cm
    for i in range(Nx):  # On parcourt les lignes (x)
        for j in range(Ny):  # On parcourt les colonnes (y)
            if (y[j] >= 300 and y[j] <= 1700) and (500 >= x[i] or x[i] >= 900):  # Condition pour l'âme
                C_num[i, j, :] = C_s  
            if (y[j] <= 300) and (300 >= x[i] or  x[i] >= 1100):  # Condition semelle inferieure
                C_num[i, j, :] = C_s
            if (y[j] >= 1700) and (x[i] == 0 or  x[i] == 1400):  # Condition semelle superieure
                C_num[i, j, :] = C_s
            

    # Mise à jour de la liste concentration C[x, y, t]
    t0 = time.time()  # Temps initial processus
    for n in range(Nt - 1):
        C_num[1:Nx - 1, 1:Ny - 1, n + 1] = (
            phi_0 * C_num[1:Nx - 1, 1:Ny - 1, n]
            + (phi_i * C_num[:Nx - 2, 1:Ny - 1, n])
            + (phi_i * C_num[2:Nx, 1:Ny - 1, n])
            + (phi_i * C_num[1:Nx - 1, :Ny - 2, n])
            + (phi_i * C_num[1:Nx - 1, 2:Ny, n])
        )
        C_num[:, -1, :] = C_s  # Condition en C[x, 2000, t]
        C_num[:, 0, :] = C_s  # Condition en C[x, 0, t]

        # Poutr en en I semelle de 30 cm
        for i in range(Nx):  # On parcourt les lignes (x)
            for j in range(Ny):  # On parcourt les colonnes (y)
                if (y[j] >= 300 and y[j] <= 1700) and (500 >= x[i] or x[i] >= 900):  # Condition pour l'âme
                    C_num[i, j, :] = C_s  
                if (y[j] <= 300) and (300 >= x[i] or  x[i] >= 1100):  # Condition semelle inferieure
                    C_num[i, j, :] = C_s
                if (y[j] >= 1700) and (x[i] == 0 or  x[i] == 1400):  # Condition semelle superieure
                    C_num[i, j, :] = C_s
                
        t1 = time.time()  # Temps final processus

    ### Visualisation des résultats
    # 0 ans
    #plt.figure(figsize=(8, 6))
    plt.imshow(np.transpose(C_num[:, :, 0]), cmap='Blues', interpolation='nearest', origin='lower')
    plt.title("Concentration de chlorure après {} années".format(0))
    plt.colorbar(label="Concentration (%)")
    plt.xlabel("Position en x (cm)")
    plt.ylabel("Position en y (cm)")
    # Ajout des contours de I
    plt.plot([50,30, 30, 110, 110,90], [30,30,0, 0, 30, 30], color='red')  # Semelle inférieure
    plt.plot([50,0,0,140, 140,90], [170,170,200,200,170,170], color='red')  # Semelle supérieure
    plt.plot([50, 50], [30, 170], color='red')  # Âme
    plt.plot([90, 90], [30, 170], color='red')  # Âme
    plt.savefig("figures/2D_C_0.pdf")
    plt.show()

    # Moteir temps
    #plt.figure(figsize=(8, 6))
    plt.imshow(np.transpose(C_num[:, :, int(Nt/2)]), cmap='Blues', interpolation='nearest',origin='lower')
    plt.title("Concentration de chlorure après {} années".format(T/2))
    plt.colorbar(label="Concentration (%)")
    plt.xlabel("Position en x (mm)")
    plt.ylabel("Position en y (mm)")
    # Ajout des contours de I
    plt.plot([50,30, 30, 110, 110,90], [30,30,0, 0, 30, 30], color='red')  # Semelle inférieure
    plt.plot([50,0,0,140, 140,90], [170,170,200,200,170,170], color='red')  # Semelle supérieure
    plt.plot([50, 50], [30, 170], color='red')  # Âme
    plt.plot([90, 90], [30, 170], color='red')  # Âme
    plt.savefig("figures/2D_C_{}.pdf".format(int(T/2)))
    plt.show()

    # T ans
    #plt.figure(figsize=(8, 6))
    plt.imshow(np.transpose(C_num[:, :, -1]), cmap='Blues', interpolation='nearest',origin='lower')
    plt.title("Concentration de chlorure après {} années".format(T))
    plt.colorbar(label="Concentration (%)")
    plt.xlabel("Position en x (mm)")
    plt.ylabel("Position en y (mm)")
    # Ajout des contours de I
    plt.plot([50,30, 30, 110, 110,90], [30,30,0, 0, 30, 30], color='red')  # Semelle inférieure
    plt.plot([50,0,0,140, 140,90], [170,170,200,200,170,170], color='red')  # Semelle supérieure
    plt.plot([50, 50], [30, 170], color='red')  # Âme
    plt.plot([90, 90], [30, 170], color='red')  # Âme
    plt.savefig("figures/2D_C_{}.pdf".format(T))
    plt.show()


    # Plot évolution armature coin 
    # 5 cm 5cm car pour C c'est en cm vu que Nx = 140 et Ny = 200
    plt.plot(t, C_num[135,195, :], label="2D : Coin supérieur droit")
    plt.plot(t, C_num[70,195, :], label="2D : Milieu supérieur")


#calcul_2D()