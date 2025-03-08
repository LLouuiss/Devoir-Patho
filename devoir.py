import numpy as np
from scipy.special import erf,erfinv
import matplotlib.pyplot as plt


def x_compact(K,t) :
    """ Fonction qui calcule la pénétration de la carbonatation du béton compact en fonction du temps"""
    return K * np.sqrt(t)

def x_fissure(w,t) : 
    """ Fonction qui calcule la pénétration de la carbonatation du béton fissuré en fonction du temps"""
    return 50 * np.sqrt(w) * np.sqrt(np.sqrt(t))

def C (x,Cs,Dce,t) : 
    """ Fonction qui calcule la concentration totale des chlorures C
        x : position armature (mm)
        Cs : concentration en surface (%)
        Dce : coefficient de diffusion des chlorures (mm^2/ans)
        t : temps (ans)
    """
    return Cs*(1-erf(x/(2*np.sqrt(Dce*t))))

def D_fissure (Dce,w,Smtheta) :
    """ Fonction qui calcule le coefficient de diffusion des chlorures dans le béton fissuré
        Dce : coefficient de diffusion des chlorures du béton non fissuré (mm^2/ans)
        w : largeur de la fissure (mm)
        Smtheta : distance entre les fissures (mm)
        Dcr : coefficient de diffusion des chlorures entre les fissures (m^2/ans)
    """
    Dcr = 5e-10 # [m^2/s], suposistion (slide du cours)
    Dcr = Dcr * 1e6 # [mm^2/s]
    Dcr = Dcr * 365*24*60*60 # [mm^2/ans]
    return Dce + (w/Smtheta)*Dcr
         
def As_reduction (Ns,d0,tcarbo,tchlor,t,vcorr1,vcorr2,plot = False) :
    """ Fonction qui calcule la réduction de la section d'armature
        Ns : nombre de barres d'armature
        d0 : diamètre des barres d'armature (mm)
        tcarbo : temps de carbonatation (ans)
        tchlor : temps d'initiation chlorure (ans)
        t : temps (ans)
        vcorr1 : vitesse de corrosion pour la carbonatation (mm/ans)
        vcorr2 : vitesse de corrosion pour les chlorures (mm/ans)
    """
    if tcarbo < tchlor :
        t0 = tcarbo
        print("Carbonatation Frist")
    else :
        t0 = tchlor
        print("Chlorures First")

    As = np.zeros(len(t))      
    for i in range(len(t)) : 
        if t[i] <= t0 :
            As[i] = Ns * np.pi * d0**2 / 4
        if t[i] > t0 and t[i] < tchlor:
            As[i] =  Ns * np.pi * (d0 - 2*vcorr1*(t[i]-t0))**2 / 4
        if t[i] > t0 and t[i] >= tchlor :
            As[i] =  Ns * np.pi * (d0 -2*vcorr1*(tchlor-t0) - 2*vcorr2*(t[i]-tchlor))**2 / 4
    
    if plot :
        plt.plot(t,As,label="Section d'armature")
        plt.vlines(tcarbo,0,Ns * np.pi * d0**2 / 4,'red',linestyles='dashed',label="Temps de carbonatation")
        plt.text(tcarbo,-1000,str(int(round(t0,0))),color='red')
        plt.vlines(tchlor,0,Ns * np.pi * (d0 - 2*vcorr1*(tchlor-t0))**2 / 4,'brown',linestyles='dashed',label="Chlorures + Carbonatation")
        plt.text(tchlor,-1000,str(int(round(tchlor,0))),color='brown')
        plt.xlabel("Temps [ans]")
        plt.ylabel("Section d'armature [mm^2]")
        plt.title("Réduction de la section d'armature en fonction du temps")
        plt.legend()
        plt.grid()
        plt.show()
    return As

def find_ti_carbonatation (K,e,w,plot = False) :
    """ Fonction qui trouve temps d'initiation de la carbonatation
    K : coefficient de diffusion de la carbonatation (mm/sqrt(ans))
    e : enrobage (mm)
    w : largeur de la fissure (mm)
    """
    # Béton compact
    t_compact = (e/K) ** 2 # [ans]
    # Béton fissuré
    t_fissure = (e/(50 * np.sqrt(w)))**4 # [ans]
    # Plot
    if plot : 
        plt.plot(t_compact,x_compact(K,t_compact),'or',markersize = 5,label="Temps d'initiation béton compacte")
        plt.hlines(e,0,max(t),'black',linestyles='dashed')
        plt.vlines(t_compact,0,e,'r',linestyles='dashed')
        plt.text(t_compact,-5,str(int(round(t_compact,0))),color='r')
        plt.plot(t_fissure,x_fissure(w,t_fissure),'ob',markersize = 5,label="Temps d'initiation béton fissuré")
        plt.vlines(t_fissure,0,e,'b',linestyles='dashed')
        plt.text(t_fissure,-5,str(int(round(t_fissure,0))),color='b')
        plt.plot(t,x_compact(K,t),label="Béton compacte")
        plt.plot(t,x_fissure(w,t),label="Béton fissuré")
        plt.xlabel("Temps [ans]")
        plt.ylabel("Pénétration de la carbonatation[mm]")
        plt.title("Pénétration de la carbonatation en fonction du temps")
        plt.grid()
        plt.legend()
        plt.show()
    
    return t_compact, t_fissure
    
def find_ti_corrosion (e,Dce,Clim,Cs,w,Sm0,plot =  False) :
    """ Fonction qui trouve temps d'initiation de la corrosion
    e : enrobage (mm)
    Dce : coefficient de diffusion des chlorures (mm^2/ans)
    Clim : concentration limite en chlorures (%)
    Cs : concentration en surface (%)
    w : largeur de la fissure (mm)
    Sm0 : distance entre les fissures (mm)
    """
    t_compact = (e/(2*np.sqrt(Dce)* erfinv((Cs-Clim)/Cs)))**2
    t_fissure = (e/(2*np.sqrt(D_fissure(Dce,w,Sm0))* erfinv((Cs-Clim)/Cs)))**2
    # Plot 
    if plot : 
        plt.plot(t,C(x,Cs,Dce,t),label="Béton compacte")
        plt.plot(t,C(x,Cs,D_fissure(Dce,w,Sm0),t),label="Béton fissuré")
        plt.plot(t,Clim*np.ones(len(t)),'--',color = 'black',label="Limite de concentration du ciment")
        plt.plot(t_compact,C(x,Cs,Dce,t_compact),'or',markersize = 5,label="Temps d'initiation béton compacte")
        plt.vlines(t_compact,0,Clim,'r',linestyles='dashed')
        plt.text(t_compact,-0.003,str(int(round(t_compact,0))),color='r')
        plt.plot(t_fissure,C(x,Cs,D_fissure(Dce,w,Sm0),t_fissure),'ob',markersize = 5,label="Temps d'initiation béton fissuré")
        plt.vlines(t_fissure,0,Clim,'b',linestyles='dashed')
        plt.text(t_fissure,-0.003,str(int(round(t_fissure,0))),color='b')
        plt.xlabel("Temps [ans]")
        plt.ylabel("Concentration en chlorures [%]")
        plt.title("Concentration en chlorures en fonction du temps")
        plt.grid()
        plt.legend()
        plt.show()
    
    
    return t_compact, t_fissure

def V_rd(As_s,fywd,d,cot_theta):
    """ Fonction qui calcule la résistance à l'effort tranchant
    As_s : section d'armature transversale par longueur [mm^2/m]
    fywd : résistance acier avec coefficient sécurité [MPa] fy/1,15
    cot_theta : inclinaison bielle comprimée 
    d : distance entre armature béton exterieur (h-e) [m]
    """
    z = 0.9 * d # Méthode simplifiée distance entre armature et résultante force béton comprimé [m]
    V_rd = As_s*fywd*z*cot_theta # [N]
    return V_rd

def M_rd(As, fyd,d):
    """ Fonction qui calcule le moment résistant de la section
    d : distance entre armature béton exterieur (h-e) [m]
    As : aire d'armature [mm^2/m]
    fyd : limite élasticité acier avec coefficent sécurité [MPa]
    """
    z = 0.9 * d # [m] méthode simplifiée
    F_s = As * fyd # [N]
    M_rd = F_s * z #[Nm]
    return M_rd

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
    
    return W_k*1000, s_r_max*1000


if __name__ == "__main__" :
    """ Données du problème """
    K = 7 #[mm/sqrt(ans)] coefficient de diffusion de la carbonatation
    t = np.linspace(0,200,201) # [ans]
    #w = 0.3 # [mm]
    Clim = 0.4 # [%] concentration limite en chlorures
    Clim = Clim  * 0.15 # [%] car dans le béton C40 y a 15% de ciment
    Dce = 5e-13 # [m^2/s] coefficient de diffusion des chlorures dans le béton non fissuré
    Dce = Dce * 1e6 # [mm^2/s]
    Dce = Dce * 365*24*60*60 # [mm^2/ans]
    Cs = 0.15 # [%] concentration en surface
    
    """Géométrie de la section"""
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
    M = 7920000  # Nm
    k_t = 0.4 #dans le cas d un chargement de longue durée
    f_ctm = 3.51e6 #Pa beton C40/50
    e = 0.05 #m enrobage des armatures
    k1 = 0.8 #pour des barres à haute adhérence
    k2 = 0.5 #en flexion simple
    phi = 0.032 #m diamètre des armatures
        
    y_G_I, I_I,alpha_e = calcul_caracteristiques_section_non_fissuree(b1, h1,b2,h2,b3,h3, d, d_prime, A_s, A_s_prime, E_cm, phi_eff, E_s)
    y_G_II, I_II = calcul_caracteristiques_section_fissuree(b1, h1,b2,h2,0.2, d, d_prime, A_s, A_s_prime, alpha_e)          
    w, Sm0 = calcul_ouverture_fissures(M, d, y_G_II, I_II, A_s, b3, h1+h2+h3, E_s, alpha_e,k_t,f_ctm,e,k1,k2,phi)
    print("w [mm], si <0.3 -> OK = ",w)
    
    e = 50 # [mm] enrobage
    t_compact1, t_fissure1 = find_ti_carbonatation(K,e,w)
    t_compact2, t_fissure2 = find_ti_corrosion(e,Dce,Clim,Cs,w,Sm0)
    
    print("Temps d'initiation de la corrosion : ",min(t_compact1,t_compact2,t_fissure1,t_fissure2)," ans")
    
    vcorr_c = 2 # [microm/ans] si seulement carbonatation
    vcorr_c = vcorr_c * 1e-3 # [mm/ans]
    vcorr_cc = 40 # [micron/ans] si carbonatation et chlorures
    vcorr_cc = vcorr_cc * 1e-3 # [mm/ans]
    t0 = min(t_compact1,t_compact2,t_fissure1,t_fissure2)
    
    # Plot de la réduction de la section d'armature pour ferraillage logitudinal et ferraillage étrier
    Ns_L = 17 # nombre de barres d'armature longitudinales
    d0_L = 40 # [mm] diamètre des barres d'armature longitudinales
    As_L = As_reduction(Ns_L,d0_L,min(t_compact1,t_fissure1),min(t_compact2,t_fissure2),t,vcorr_c,vcorr_cc)
    
    #### ATTENTION POUR ETRIER ON VA DEVOIR DIVISER AS EN DEUX CAR FCT PREND EN COMPTE SEULEMENT UNE SECTION D'ARMATURE ####
    e= 50 # [mm] enrobage
    fywd = 435 # [MPa] pour acier BE500S
    cot_theta = 2 # inclinaison bielle comprimée
    d = 2 - e*(10**-3) # [m] distance entre armature béton extérieur ZONE 2 => 2 m
    Ns_E = 2*2 # 2 etriers de 2 barres
    d0_E = 11 # [mm] diamètre des barres d'armature étriers
    As_E = As_reduction(Ns_E,d0_E,min(t_compact1,t_fissure1),min(t_compact2,t_fissure2),t,vcorr_c,vcorr_cc)
    As_Es = As_E / 0.11 # [mm^2/m] 2 etriers tout les 11 cm
    V_rd_values = V_rd(As_Es,fywd,d,cot_theta)
    plot = plt.plot(t,V_rd_values,label="Résistance à l'effort tranchant")
    plt.xlabel("Temps [ans]")
    plt.ylabel("Résistance à l'effort tranchant [N]")
    plt.title("Résistance à l'effort tranchant en fonction du temps")
    plt.legend()
    plt.grid()
    plt.show()
    
    
    