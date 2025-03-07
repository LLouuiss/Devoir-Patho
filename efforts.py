import matplotlib.pyplot as plt
import numpy as np

def M_rd(As1,As2, fyd,z):
    
    F_s = (As1 + As2) * fyd*10e-5
    M_rd = F_s * z
    return M_rd

def V_rd(As_s,fyd,z,cot_theta):
    
    V_rd = As_s*fyd*z*cot_theta*10e-5
    return V_rd
A_s1 = np.array([213.71,200,186.29,171.43,150,128.57,100,50])
A_s2 = np.array([42.83,40,37.26,34.29,30,25.71,20,10])
As_s = np.array([20.1,20,19,18,15,12,8,4])
t = np.array([0,5,10,15,20,25,30,35])
V_ed = 1275
M_ed = 15500
cot_theta = 2
fyd = 435e6
z_M = 1.7
z_v = 0.731

# Calcul des valeurs de M_rd et V_rd
M_rd_values = M_rd(A_s1,A_s2, fyd, z_M)*10e-4
V_rd_values = V_rd(As_s, fyd, z_v, cot_theta)*10e-4
print("M_rd = ",M_rd_values)
print("V_rd =",V_rd_values)
# Tracé du graphe de M_rd en fonction de t
plt.figure(figsize=(10, 5))
plt.plot(t, M_rd_values, label='M_rd', marker='o')
plt.axhline(y=M_ed, color='r', linestyle='--', label='M_ed')
plt.xlabel('t (années)')
plt.ylabel('M_rd (kN.m)')
plt.title('M_rd en fonction du temps')
plt.legend()
plt.grid()
plt.show()

# Tracé du graphe de V_rd en fonction de t
plt.figure(figsize=(10, 5))
plt.plot(t, V_rd_values, label='V_rd', marker='o')
plt.axhline(y=V_ed, color='r', linestyle='--', label='V_ed')
plt.xlabel('t (années)')
plt.ylabel('V_rd (kN)')
plt.title('V_rd en fonction du temps ')
plt.legend()
plt.grid()
plt.show()