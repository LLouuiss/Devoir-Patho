def calcul_caracteristiques_section_non_fissuree(b1, h1,b2,h2,b3,h3, d, d_prime, A_s, A_s_prime, E_cm, phi_eff, E_s):
    """
    Calcule les caractéristiques mécaniques de la section non fissurée (y_G,I et I_I)
    
    :param b1: Largeur de la section en fibre inférieure  (m)
    :param h1: Hauteur totale de la section (m)
    :param d: Hauteur utile (m)
    :param d_prime: Hauteur des armatures comprimées (m)
    :param A_s: Aire des armatures tendues (m^2)
    :param A_s_prime: Aire des armatures comprimées (m^2)
    :param E_cm: Module d'élasticité du béton (GPa)
    :param phi_eff: Coefficient de fluage effectif
    :param E_s: Module d'élasticité de l'acier (GPa)
    :return: y_G,I (m), I_I (m^4), alpha_e ()
    """
    # Calcul du module d'élasticité effectif du béton
    E_c_eff = E_cm / (1 + phi_eff)  # GPa
    
    # Calcul du coefficient d'équivalence
    alpha_e = E_s / E_c_eff
   
    # Calcul du centre de gravité de la section non fissurée
    
    y_G_I = (b1*h1*h1/2+b2*h2*(h1+h2/2)+b3*h3*(h1+h2+h3/2) + alpha_e * A_s * d + alpha_e * A_s_prime * d_prime) / (b1 * h1+b2*h2+b3*h3 + alpha_e * A_s + alpha_e * A_s_prime)
    
    # Calcul de l'inertie de la section non fissurée
    I_I = (b1 * h1**3 / 12) +(b2 * h2**3 / 12) +(b3 * h3**3 / 12) + b1 * h1 * (h1/2 - y_G_I)**2+b2 * h2 * (h2/2+h1 - y_G_I)**2+b3 * h3 * (h3/2+h2+h1 - y_G_I)**2 + alpha_e * A_s * (y_G_I - d)**2 + alpha_e * A_s_prime * (y_G_I - d_prime)**2
    
    return y_G_I, I_I,alpha_e

def calcul_caracteristiques_section_fissuree(b1, h1,b2,h2,x, d, d_prime, A_s, A_s_prime, alpha_e):
    """
    Calcule les caractéristiques mécaniques de la section fissurée (y_G,II et I_II)
    
    :param b: Largeur de la section (m)
    :param x: Position de l'axe neutre (m)
    :param d: Hauteur utile (m)
    :param d_prime: Hauteur des armatures comprimées (m)
    :param A_s: Aire des armatures tendues (mm^2)
    :param A_s_prime: Aire des armatures comprimées (mm^2)
    :param alpha_e: Coefficient d'équivalence acier/béton
    :return: y_G_II (m), I_II (m^4)
    """
    
    for i in range(10):
        
        x =  (b1 * x**2/2 + alpha_e * A_s * d + alpha_e * A_s_prime * d_prime) / (b1 *x + alpha_e * A_s + alpha_e * A_s_prime)
    # Calcul de l'inertie de la section fissurée
    x = 0.88
    
    I_II = (b1 * h1**3 / 12) +(b2 * (x-h1)**3 / 12) + b1 * h1 * (h1/2 - x)**2+b2 * (x-h1) * (h2/2+h1 - x)**2 + alpha_e * A_s * (x - d)**2 + alpha_e * A_s_prime * (x - d_prime)**2
    return x, I_II
def calcul_ouverture_fissures(M_II, d, x, I_II, A_s, b, h, E_s, alpha_e, k_t, f_ctm,c,k1,k2,phi) :
    """Calcule l'ouverture des fissures à l'ELS quasi-permanent
    :param M_II: Moment fléchissant (kNm)
    :param d: Hauteur utile (m)
    :param x: Position de la fibre neutre (m)
    :param I_II: Moment d'inertie (m^4)
    :param A_s: Aire des armatures (mm^2)
    :param b: Largeur de la section (mm)
    :param h: Hauteur totale de la section (mm)
    :param E_s: Module d'élasticité de l'acier (MPa)
    :param alpha_e: Facteur de correction
    :param k_t: Facteur de durée de charge
    :param f_ctm: Résistance en traction du béton (MPa)
    :param c: enrobage des armatures (m)
    :param k1: Paramètre de fissuration
    :param k2: Paramètre de fissuration
    :param phi: Diamètre des armatures (mm)
    :return: Ouverture des fissures W_k (mm)
    """
    
    
    # Calcul de la contrainte dans les armatures tendues (sigma_s)
    sigma_s = alpha_e * (M_II / I_II) * (d - x)  # MPa
    print("d = ",d)
    print("x = ",x)
    print("sigma_s = ",sigma_s)
    # Calcul de la hauteur efficace du béton autour des armatures
    h_ceff = 2.5 * (h - d)  # mm
    print("h_ceff = ",h_ceff)
    # Calcul du taux d'armatures efficaces
    rho_seff = A_s / (h_ceff * b)
    print("rho_seff = ",rho_seff)
    # Calcul de la déformation epsilon_sm - epsilon_cm
    terme1 = sigma_s / E_s
    print("terme1 = ",0.6*terme1)

    terme2 = k_t * (f_ctm / (E_s * rho_seff)) * (1 + alpha_e * rho_seff)
   
    epsilon_sm_cm = max(terme1 * 0.6, terme1 - terme2)
    print("epsilon_sm_cm = ",epsilon_sm_cm)
    
    # Espacement des fissures (approximé selon des formules usuelles)
    s_r_max = 3.4*c + 0.425*k1*k2*phi/rho_seff  # Approche simplifiée
    print("s_r_max = ",s_r_max)
    # Ouverture des fissures
    W_k = s_r_max * epsilon_sm_cm  # mm
    
    return W_k*1000


#Paramètres du problème
b1 = 0.8  # m
h1 = 0.3  # m
b2 = 0.4
h2 = 1.4
b3 = 1.4
h3 = 0.3
d = 1.85  # m
d_prime = 0.15  # m
A_s = 225.12e-4  # cm^2 28phi32
A_s_prime = 44.19e-4  # m^2 9phi25
E_cm = 35.2e9  # Pa béton C40/50
phi_eff = 2
E_s = 200e9  # Pa
M = 15300000  # Nm
k_t = 0.4 #dans le cas d un chargement de longue durée
f_ctm = 3.51e6 #Pa beton C40/50
c = 0.05 #m enrobage des armatures
k1 = 0.8 #pour des barres à haute adhérence
k2 = 0.5 #en flexion simple
phi = 0.032 #m diamètre des armatures
y_G_I, I_I,alpha_e = calcul_caracteristiques_section_non_fissuree(b1, h1,b2,h2,b3,h3, d, d_prime, A_s, A_s_prime, E_cm, phi_eff, E_s)
y_G_II, I_II = calcul_caracteristiques_section_fissuree(b1, h1,b2,h2,0.2, d, d_prime, A_s, A_s_prime, alpha_e)
print("y_G_I = ",y_G_I)
print("I_I = ",I_I)
print("alpha_e = ",alpha_e)
print("y_G_II = ",y_G_II)
print("I_II = ",I_II)
w = calcul_ouverture_fissures(M, d, y_G_II, I_II, A_s, b3, h1+h2+h3, E_s, alpha_e,k_t,f_ctm,c,k1,k2,phi)
print("w [mm], si <0.3 -> OK = ",w)