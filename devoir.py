import numpy as np
from scipy.special import erf,erfinv
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.optimize import fsolve

from deuxD import calcul_2D

def find_intersection(x):
    """ Fonction qui trouve l'intersection de deux fonctions"""
    return f1(x) - f2(x)

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
         
def As_reduction (Ns,d0,tcarbo,tchlor,t,vcorr1,vcorr2,plot = False,titre = None) :
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
            if (d0 -2*vcorr1*(tchlor-t0) - 2*vcorr2*(t[i]-tchlor)) < 0: 
                As[i] = 0
            else:
                As[i] =  Ns * np.pi * (d0 -2*vcorr1*(tchlor-t0) - 2*vcorr2*(t[i]-tchlor))**2 / 4
    
    if plot :
        plt.plot(t,As,label="Section d'armature")
        plt.vlines(tcarbo,0,Ns * np.pi * d0**2 / 4,'red',linestyles='dashed',label="Temps de carbonatation")
        plt.text(tcarbo,-15,str(int(round(t0,0))),color='red')
        plt.vlines(tchlor,0,Ns * np.pi * (d0 - 2*vcorr1*(tchlor-t0))**2 / 4,'brown',linestyles='dashed',label="Chlorures + Carbonatation")
        plt.text(tchlor,-15,str(int(round(tchlor,0))),color='brown')
        plt.xlabel("Temps [ans]")
        plt.ylabel("Section d'armature [mm^2]")
        plt.title("Réduction de la section d'armature en fonction du temps")
        plt.legend()
        plt.grid()
        if titre != None :
            plt.savefig(titre)
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
    #Béton compact avec CO2
    t_compact_CO2 = (e/(K*2**0.5)) ** 2 # [ans]
    # Béton fissuré
    t_fissure = (e/(50 * np.sqrt(w)))**4 # [ans]
    # Plot
    if plot : 
        plt.plot(t_compact,x_compact(K,t_compact),'or',markersize = 5,label="Temps d'initiation béton compact")
        plt.hlines(e,0,max(t),'black',linestyles='dashed')
        plt.vlines(t_compact,0,e,'r',linestyles='dashed')
        plt.text(t_compact,-5,str(int(round(t_compact,0))),color='r')
        plt.plot(t_fissure,x_fissure(w,t_fissure),'ob',markersize = 5,label="Temps d'initiation béton fissuré")
        plt.vlines(t_fissure,0,e,'b',linestyles='dashed')
        plt.text(t_fissure,-5,str(int(round(t_fissure,0))),color='b')
        plt.plot(t,x_compact(K,t),label="Béton compact")
        plt.plot(t,x_fissure(w,t),label="Béton fissuré")
        plt.xlabel("Temps [ans]")
        plt.ylabel("Pénétration de la carbonatation[mm]")
        plt.title("Pénétration de la carbonatation en fonction du temps")
        plt.grid()
        plt.legend()
        plt.savefig("figures/carbonatation.pdf")
        plt.show()
        
        #plot influence Co2
        plt.plot(t_compact,x_compact(K,t_compact),'ob',markersize = 5,label="Temps d'initiation béton compact : K = "  + str(K))
        plt.hlines(e,0,max(t),'black',linestyles='dashed',label="Enrobage")
        plt.vlines(t_compact,0,e,'b',linestyles='dashed')
        plt.text(t_compact,-5,str(int(round(t_compact,0))),color='b')
        plt.plot(t_compact_CO2,x_compact(K*2**0.5,t_compact_CO2),'o',markersize = 5, color = 'orange',label="Temps d'initiation béton compact : K = "  + str(9.9))
        plt.hlines(e,0,max(t),'black',linestyles='dashed')
        plt.vlines(t_compact_CO2,0,e,'orange',linestyles='dashed')
        plt.text(t_compact_CO2,-5,str(int(round(t_compact_CO2,0))),color='orange')
        plt.plot(t_fissure,x_fissure(w,t_fissure),'og',markersize = 5,label="Temps d'initiation béton fissuré")
        plt.vlines(t_fissure,0,e,'g',linestyles='dashed')
        plt.text(t_fissure,-5,str(int(round(t_fissure,0))),color='g')
        plt.plot(t,x_compact(K,t),label="Béton compact : K = "  + str(K))
        plt.plot(t,x_compact(K*2**0.5,t),label="Béton compact : K = "  + str(9.9))
        plt.plot(t,x_fissure(w,t),label="Béton fissuré")
        plt.xlabel("Temps [ans]")
        plt.ylabel("Pénétration de la carbonatation[mm]")
        plt.title("Pénétration de la carbonatation en fonction du temps")
        plt.grid()
        plt.legend()
        plt.savefig("figures/carbonatation_CO2.pdf")
        plt.show()
    
    return t_compact, t_fissure
    
def find_ti_chlorure (e,Dce,Clim,Cs,w,Sm0,plot =  False,deux_D = False) :
    """ Fonction qui trouve temps d'initiation de la corrosion par les chlorures
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
        if deux_D : 
            calcul_2D()
        plt.plot(t,C(e,Cs,Dce,t),label="1D : Béton compact")
        plt.plot(t,C(e,Cs,D_fissure(Dce,w,Sm0),t),label="1D : Béton fissuré")
        plt.plot(t,Clim*np.ones(len(t)),'--',color = 'black',label="Limite de concentration du ciment")
        plt.plot(t_compact,C(e,Cs,Dce,t_compact),'or',markersize = 5,label="Temps d'initiation béton compact")
        plt.vlines(t_compact,0,Clim,'r',linestyles='dashed')
        plt.text(t_compact,-0.003,str(int(round(t_compact,0))),color='r')
        plt.plot(t_fissure,C(e,Cs,D_fissure(Dce,w,Sm0),t_fissure),'ob',markersize = 5,label="Temps d'initiation béton fissuré")
        plt.vlines(t_fissure,0,Clim,'b',linestyles='dashed')
        plt.text(t_fissure,-0.003,str(int(round(t_fissure,0))),color='b')
        plt.xlabel("Temps [ans]")
        plt.ylabel("Concentration en chlorures [%]")
        plt.title("Concentration en chlorures en fonction du temps")
        plt.grid()
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=10)
        plt.savefig("figures/chlorure_2D.pdf",bbox_inches='tight')
        plt.show()
    
    
    
    return t_compact, t_fissure

def V_rd(As_s,fywd,d,cot_theta):
    """ Fonction qui calcule la résistance à l'effort tranchant
    As_s : section d'armature transversale par longueur [mm^2/m]
    fywd : résistance acier avec coefficient sécurité [MPa] fy/1,15
    cot_theta : inclinaison bielle comprimée 
    d : distance entre armature béton exterieur (h-e) [m]
    """
    z = 0.9 * d  # Méthode simplifiée distance entre armature et résultante force béton comprimé [m]
    V_rd = As_s*fywd*z*cot_theta # [N]
    return V_rd

def M_rd(As, fyd,d):
    """ Fonction qui calcule le moment résistant de la section
    d : distance entre armature béton exterieur (h-e) [mm]
    As : aire d'armature [mm^2]
    fyd : limite élasticité acier avec coefficent sécurité [MPa]
    """
    z = 0.9 * d # [mm] méthode simplifiée
    F_s = As * fyd # [N]
    M_rd = F_s * z #[Nmm]
    return M_rd

def calcul_caracteristiques_section_non_fissuree(b1, h1,b2,h2,b3,h3, d, d_prime, A_s, A_s_prime, E_cm, phi_eff, E_s):
    """
    Calcule les caractéristiques mécaniques de la section non fissurée (y_G,I et I_I)
    b1,h1,b2,h2,b3,h3 ;  géométrie de la section (mm)
    d : Hauteur armature tendue (mm)
    d_prime : Hauteur des armatures comprimées (mm)
    A_s : Aire des armatures tendues (mm^2)
    A_s_prime : Aire des armatures comprimées (mm^2)
    E_cm : Module d'élasticité du béton (MPa)
    phi_eff : Coefficient de fluage effectif du béton
    E_s : Module d'élasticité de l'acier (MPa)
    return : y_G_I (mm), I_I (mm^4) et alpha_e
    """
    # Calcul du module d'élasticité effectif du béton
    E_c_eff = E_cm / (1 + phi_eff)  # GPa
    
    # Calcul du coefficient d'équivalence
    alpha_e = E_s / E_c_eff
   
    # Calcul du centre de gravité de la section non fissurée
    # Réference = bas de la poutre
    y_G_I = (b1*h1*h1/2+b2*h2*(h1+h2/2)+b3*h3*(h1+h2+h3/2) + alpha_e * A_s * d + alpha_e * A_s_prime * d_prime) / (b1 * h1 + b2*h2 + b3*h3 + alpha_e * A_s + alpha_e * A_s_prime)
    
    # Calcul de l'inertie de la section non fissurée
    I_I = (b1 * h1**3 / 12) +(b2 * h2**3 / 12) +(b3 * h3**3 / 12) + b1 * h1 * (h1/2 - y_G_I)**2+b2 * h2 * (h2/2+h1 - y_G_I)**2+b3 * h3 * (h3/2+h2+h1 - y_G_I)**2 + alpha_e * A_s * (y_G_I - d)**2 + alpha_e * A_s_prime * (y_G_I - d_prime)**2
    
    return y_G_I, I_I,alpha_e

def calcul_caracteristiques_section_fissuree(b1, h1,b2,h2, d, d_prime, A_s, A_s_prime, alpha_e):
    """
    Calcule les caractéristiques mécaniques de la section fissurée (y_G,II et I_II)
    b1,h1,b2,h2 ;  géométrie de la section (mm) là ou c'est pas fissuré
    d : Hauteur armature tendue (mm)
    d_prime : Hauteur des armatures comprimées (mm)
    A_s : Aire des armatures tendues (mm^2)
    A_s_prime : Aire des armatures comprimées (mm^2)
    alpha_e : coefficient d'équivalence
    return : y_G_II (mm) et I_II (mm^4) par rapport au bas/haut de h1
    """
    # Valeur initiale de x
    x = h1
    for i in range(10):
        # Hypothèse position de l'axe neutre est dans la semelle
        #x =  ((b1 * x**2)/2 + alpha_e * A_s * d + alpha_e * A_s_prime * d_prime) / (b1 *x + alpha_e * A_s + alpha_e * A_s_prime)
        # Hypothèse position de l'axe neutre est dans l'âme
        x =  (((b1 * h1**2)/2 + b2 * (x-h1) * (h1 +(x-h1)/2 )) + alpha_e * A_s * d + alpha_e * A_s_prime * d_prime) / (b1 * h1 + b2 * (x-h1) + alpha_e * A_s + alpha_e * A_s_prime)
        #print("x = ",x)
    # Calcul de l'inertie de la section fissurée
    # Hypothèse position de l'axe neutre est dans l'âme
    I_II = (b1 * h1**3 / 12) +(b2 * (x-h1)**3 / 12) + b1 * h1 * (h1/2 - x)**2 + b2 * (x-h1) * ((x-h1)/2 + h1 - x)**2 + alpha_e * A_s * (x - d)**2 + alpha_e * A_s_prime * (x - d_prime)**2
    return x, I_II


def calcul_ouverture_fissures(M_II, d, x, I_II, A_s, b, h, E_s, alpha_e, k_t, f_ctm,c,k1,k2,phi,print_values = False) :
    """Calcule l'ouverture des fissures à l'ELS quasi-permanent
    M_II: Moment fléchissant (Nmm)
    d : hauteur armature tendue (mm)
    x: Position de la fibre neutre (mm)
    I_II: Moment d'inertie section fissurée (mm^4)
    A_s: Aire des armatures tendue (mm^2)
    b: Largeur de la section (mm)
    h: Hauteur totale de la section (mm)
    E_s: Module d'élasticité de l'acier (MPa)
    alpha_e: coefficient d'équivalence
    k_t: Facteur de durée de charge
    f_ctm: Résistance en traction du béton (MPa)
    c: enrobage des armatures (mm)
    k1: Paramètre de fissuration
    k2: Paramètre de fissuration
    phi: Diamètre des armatures (mm)
    
    """
    
    ### d et x doivent avoir meme ref
    
    # Calcul de la contrainte dans les armatures tendues (sigma_s)
    sigma_s = alpha_e * (M_II / I_II) * (d - x)  # MPa
    # Calcul de la hauteur efficace du béton autour des armatures
    h_ceff = min(2.5 * (h - d),(h-x)/3, h/2 ) # mm
    
    ##### TO EDIT #####
    rho_seff = A_s / (h_ceff * b)
    # Calcul de la déformation epsilon_sm - epsilon_cm
    terme1 = sigma_s / E_s

    terme2 = k_t * (f_ctm / (E_s * rho_seff)) * (1 + alpha_e * rho_seff)
   
    epsilon_sm_cm = max(terme1 * 0.6, terme1 - terme2)
    
    # Espacement des fissures (approximé selon des formules usuelles)
    s_r_max = 3.4*c + 0.425*k1*k2*phi/rho_seff  # Approche simplifiée (mm)
    #print("s_r_max = ",s_r_max)
    # Ouverture des fissures
    w_k = s_r_max * epsilon_sm_cm  # mm
    if print_values :
        print("d = ",d)
        print("x = ",x)
        print("sigma_s = ",sigma_s)
        print("h_ceff = ",h_ceff)
        print("terme1 = ",0.6*terme1)
        print("epsilon_sm_cm = ",epsilon_sm_cm)
    
    return w_k, s_r_max


if __name__ == "__main__" :
    """ Données du problème """
    plot = False   
    K = 7 #[mm/sqrt(ans)] coefficient de diffusion de la carbonatation
    t = np.linspace(1,150,152) # [ans]
    Clim = 0.4 # [%] concentration limite en chlorures
    Clim = Clim  * 0.15 # [%] car dans le béton C40 y a 15% de ciment
    Dce = 5e-13 # [m^2/s] coefficient de diffusion des chlorures dans le béton non fissuré
    Dce = Dce * 1e6 # [mm^2/s]
    Dce = Dce * 365*24*60*60 # [mm^2/ans]
    Cs = 0.15 # [%] concentration en surface
    vcorr_c = 2 # [microm/ans] si seulement carbonatation
    vcorr_c = vcorr_c * 1e-3 # [mm/ans]
    vcorr_cc = 40 # [micron/ans] si carbonatation et chlorures
    vcorr_cc = vcorr_cc * 1e-3 # [mm/ans]
    
    # Efforts sollicitants
    M = 7920000  # Nm
    M = M * 1e3 # Nmm
    V =  1320000  # N
    
    # Caractéristiques béton et acier
    fywd = 435 # [MPa] pour acier BE500S
    f_ctm = 3.51e6 # Pa beton C40/50
    f_ctm = f_ctm * 1e-6 # MPa
    E_s = 200e9  # Pa Module d'élasticité de l'acier
    E_s = E_s * 1e-6 # MPa 
    E_cm = 35.2e9  # Pa béton C40/50
    E_cm = E_cm * 1e-6 # MPa
    phi_eff = 2 # coefficient de fluage effectif du béton
    
    #enrobage  = np.linspace(50,5,10) # [mm] enrobage
    #enrobage  = np.linspace(80,5,16) # [mm] enrobage
    #enrobage = np.linspace(150,5,30) # [mm] enrobage
    enrobage = [50]
    solution_M = []
    solution_M_NL = []
    solution_M_compact = []
    solution_V = []
    solution_V_NL = []
    solution_V_compact = []
    temps_cl_NL = []
    temps_cl_L = []
    for i in enrobage :
        print("Enrobage = ",i)
        e = i # [mm] enrobage
        """Géométrie de la section au nu de la colonne """
        b1 = 800 # mm
        h1 = 300 # mm  
        b2 = 400 # mm 
        h2 = 1400 # mm
        b3 = 1400 # mm
        h3 = 300 # mm
        d = 2000 - 150   # mm
        d_prime = 150 # mm
        k_t = 0.4 # dans le cas d un chargement de longue durée
        k1 = 0.8 # pour des barres à haute adhérence
        k2 = 0.5 # en flexion simple


        A_s =  28 * 32**2 * np.pi/4 # mm^2 28phi32 en traction
        A_s_prime = 9 * 25**2 * np.pi/4  # mm^2 9phi25 en compression
        
        y_G_I, I_I, alpha_e = calcul_caracteristiques_section_non_fissuree(b1, h1,b2,h2,b3,h3, d, d_prime, A_s, A_s_prime, E_cm, phi_eff, E_s)
        y_G_II, I_II = calcul_caracteristiques_section_fissuree(b1, h1,b2,h2, d, d_prime, A_s, A_s_prime, alpha_e)
        phi = 32 # mm diamètre des armatures en traction        
        w, Sm0 = calcul_ouverture_fissures(M, d, y_G_II, I_II, A_s, b3, h1+h2+h3, E_s, alpha_e,k_t,f_ctm,e,k1,k2,phi)
        
        t_compact_carbo, t_fissure_carbo = find_ti_carbonatation(K,e,w,False)
        t_compact_chl, t_fissure_chl = find_ti_chlorure(e,Dce,Clim,Cs,w,Sm0,True,True)
        
        print("Temps d'initiation de la corrosion : ",min(t_compact_carbo,t_compact_chl,t_fissure_carbo,t_fissure_chl)," ans")
        
        
        t0 = min(t_compact_carbo,t_compact_chl,t_fissure_carbo,t_fissure_chl)
        
        # Plot de la réduction de la section d'armature pour ferraillage logitudinal et ferraillage étrier
        Ns_L1 = 28 # nombre de barres d'armature longitudinales 28
        d0_L1 = 32 # [mm] diamètre des barres d'armature longitudinales en traction
        Ns_L2 = 9 # nombre de barres d'armature longitudinales
        d0_L2 = 25 # [mm] diamètre des barres d'armature longitudinales en compression
        
        As_L = As_reduction(Ns_L1,d0_L1,min(t_compact_carbo,t_fissure_carbo),min(t_compact_chl,t_fissure_chl),t,vcorr_c,vcorr_cc,False,titre="figures/red_section_long_traction.pdf") # béton compact et fissuré
        As_L_linear = As_reduction(Ns_L1,d0_L1,min(t_compact_carbo,t_fissure_carbo),min(t_compact_chl,t_fissure_chl),t,vcorr_c,vcorr_cc) # béton compact et fissuré
        As_L_compact = As_reduction(Ns_L1,d0_L1,t_compact_carbo,t_compact_chl,t,vcorr_c,vcorr_cc) # béton compact
        #As_L_fissure = As_reduction(Ns_L1,d0_L1,t_fissure_carbo,t_fissure_chl,t,vcorr_c,vcorr_cc) # béton fissuré
        
        As_L2 = As_reduction(Ns_L2,d0_L2,min(t_compact_carbo,t_fissure_carbo),min(t_compact_chl,t_fissure_chl),t,vcorr_c,vcorr_cc,False, titre="figures/red_section_long_compression.pdf") # béton compact et fissuré
        As_L_linear2 = As_reduction(Ns_L2,d0_L2,min(t_compact_carbo,t_fissure_carbo),min(t_compact_chl,t_fissure_chl),t,vcorr_c,vcorr_cc) # béton compact et fissuré
        As_L2_compact = As_reduction(Ns_L2,d0_L2,t_compact_carbo,t_compact_chl,t,vcorr_c,vcorr_cc) # béton compact
        #As_L2_fissure = As_reduction(Ns_L2,d0_L2,t_fissure_carbo,t_fissure_chl,t,vcorr_c,vcorr_cc) # béton fissuré
        
        Ns_E = 2*2 # 2 etriers de 2 barres
        d0_E = 10 # [mm] diamètre des barres d'armature étriers
        As_E = As_reduction(Ns_E,d0_E,min(t_compact_carbo,t_fissure_carbo),min(t_compact_chl,t_fissure_chl),t,vcorr_c,vcorr_cc,False,titre="figures/red_section_etrier.pdf") # béton compact et fissuré
        As_E_linear = As_reduction(Ns_E,d0_E,min(t_compact_carbo,t_fissure_carbo),min(t_compact_chl,t_fissure_chl),t,vcorr_c,vcorr_cc) # béton compact et fissuré
        As_E_compact = As_reduction(Ns_E,d0_E,t_compact_carbo,t_compact_chl,t,vcorr_c,vcorr_cc) # béton compact
        #As_E_fissure = As_reduction(Ns_E,d0_E,t_fissure_carbo,t_fissure_chl,t,vcorr_c,vcorr_cc) # béton fissuré

        
        print("-----Sans non linéarité de la réduction de la section d'armature----")
        print("t chlorure compact = ",t_compact_chl)
        print("t chlorure fissure = ",t_fissure_chl)
        print("s_r_max = ",Sm0)
        print("w [mm], si <0.3 -> OK = ",w)
        temps_cl_L.append(min(t_compact_chl,t_fissure_chl))
        
        # Non linéarité de la réduction de la section d'armature
        NL = False 
        if NL :
            for i in range(len(t)) :
                if t[i] <= t0 :
                    As_L[i] = As_L[i]
                if t[i] > t0 and t[i] < min(t_compact_chl,t_fissure_chl) : 
                    y_G_II, I_II = calcul_caracteristiques_section_fissuree(b1, h1,b2,h2, d, d_prime, As_L[i], As_L2[i], alpha_e)
                    phi = np.sqrt(4*As_L[i]/(np.pi*Ns_L1)) # [mm] diamètre des armatures en traction
                    #print("phi = ",phi)
                    w, Sm0 = calcul_ouverture_fissures(M, d, y_G_II, I_II, As_L[i], b3, h1+h2+h3, E_s, alpha_e,k_t,f_ctm,e,k1,k2,phi)
                    t_compact_chl, t_fissure_chl = find_ti_chlorure(e,Dce,Clim,Cs,w,Sm0)
                    As_L[i:] = As_reduction(Ns_L1,d0_L1,min(t_compact_carbo,t_fissure_carbo),min(t_fissure_chl,t_compact_chl),t,vcorr_c,vcorr_cc)[i:]
                    As_L2[i:] = As_reduction(Ns_L2,d0_L2,min(t_compact_carbo,t_fissure_carbo),min(t_fissure_chl,t_compact_chl),t,vcorr_c,vcorr_cc)[i:]
                    As_E[i:] = As_reduction(Ns_E,d0_E,min(t_compact_carbo,t_fissure_carbo),min(t_compact_chl,t_fissure_chl),t,vcorr_c,vcorr_cc)[i:]
                if t[i] > t0 and t[i] >= min(t_compact_chl,t_fissure_chl) :
                    break
            
        print("-----Avec non linéarité de la réduction de la section d'armature----")
        print("t chlorure compact = ",t_compact_chl)
        print("t chlorure fissure = ",t_fissure_chl)
        print("w [mm], si <0.3 -> OK = ",w)
        print("s_r_max = ",Sm0)
        temps_cl_NL.append(min(t_compact_chl,t_fissure_chl))

        
        
        M_rd_values = M_rd(As_L,fywd,d) # [Nmm]
        M_rd_values_compact = M_rd(As_L_compact,fywd,d) # [Nmm]
        if plot :
            plt.plot(t,M_rd_values*10**(-6),label="Moment résistant")
            plt.plot(t,M*np.ones(len(t))*10**(-6),'--',color = 'black',label="Moment sollicitant")
            plt.xlabel("Temps [ans]")
            plt.ylabel("Moment résistant [kNm]")
            plt.title("Moment résistant en fonction du temps pour un enrobage de " + str(e) + " mm")
            plt.legend()
            plt.grid()
            plt.savefig("figures/moment_resistant_5cm.pdf",bbox_inches='tight')
            plt.show()
        # Créer des interpolations linéaires
        f1 = interp1d(t, M_rd_values, kind='linear', fill_value="extrapolate")
        f2 = interp1d(t, M*np.ones(len(t)), kind='linear', fill_value="extrapolate")
        # Trouver l'intersection numériquement
        x_guess = min(t_fissure_chl,t_compact_chl) + 10  # Estimation initiale de l'intersection
        x_intersection = fsolve(find_intersection, x_guess)[0]
        y_intersection = f1(x_intersection)
        # Afficher les coordonnées d'intersection
        print(f"Rupture M  = {x_intersection:.4f} ans ")
        solution_M.append(x_intersection)
        
        ### ----- Compact ----
        f1 = interp1d(t, M_rd_values_compact, kind='linear', fill_value="extrapolate")
        f2 = interp1d(t, M*np.ones(len(t)), kind='linear', fill_value="extrapolate")
        x_guess = t_compact_chl + 10  # Estimation initiale de l'intersection
        x_intersection = fsolve(find_intersection, x_guess)[0]
        y_intersection = f1(x_intersection)
        print(f"Rupture M compact = {x_intersection:.4f} ans")
        solution_M_compact.append(x_intersection)
        
        ## ---- Non linéarité : cas linéaire----
        if NL : 
            M_NL = M_rd(As_L_linear,fywd,d)
            f1 = interp1d(t, M_NL, kind='linear', fill_value="extrapolate")
            f2 = interp1d(t, M*np.ones(len(t)), kind='linear', fill_value="extrapolate")
            x_guess = min(t_fissure_chl,t_compact_chl) + 10   # Estimation initiale de l'intersection
            x_intersection = fsolve(find_intersection, x_guess)[0]
            y_intersection = f1(x_intersection)
            print(f"Rupture M lineaire = {x_intersection:.4f} ans")
            solution_M_NL.append(x_intersection)
            
   
        

        ## On analyse la resistence des etriers pour la zone au nu de la colonne car effort tranchant plus critique
        d = 1500 - e  # mm
        cot_theta = 2 # inclinaison bielle comprimée
        As_Es = As_E / 0.11 # [mm^2/m] 2 etriers tout les 11 cm
        As_Es_compact = As_E_compact / 0.11 # [mm^2/m] 2 etriers tout les 11 cm
        V_rd_values = V_rd(As_Es,fywd,d*10**-3,cot_theta) # [N]
        V_rd_values_compact = V_rd(As_Es_compact,fywd,d*10**-3,cot_theta) # [N]
        if plot :
            plt.plot(t,V_rd_values*10**(-3),label="Effort tranchant résistant")
            plt.plot(t,V*np.ones(len(t))*10**(-3),'--',color = 'black',label="Effort tranchant sollicitant")
            plt.xlabel("Temps [ans]")
            plt.ylabel("Effort tranchant résistant [kN]")
            plt.title("Effort tranchant résistant en fonction du temps pour un enrobage de " + str(e) + " mm")
            plt.legend()
            plt.grid()
            plt.savefig("figures/effort_tranchant_5cm.pdf",bbox_inches='tight')
            plt.show()
        f1 = interp1d(t, V_rd_values, kind='linear', fill_value="extrapolate")
        f2 = interp1d(t, V*np.ones(len(t)), kind='linear', fill_value="extrapolate")
        x_guess = min(t_fissure_chl,t_compact_chl) + 30  # Estimation initiale de l'intersection
        x_intersection = fsolve(find_intersection, x_guess)[0]
        y_intersection = f1(x_intersection)
        print(f"Rupture V = {x_intersection:.4f} ans")
        solution_V.append(x_intersection)
        
        ### ----- Compact ----
        f1 = interp1d(t, V_rd_values_compact, kind='linear', fill_value="extrapolate")
        f2 = interp1d(t, V*np.ones(len(t)), kind='linear', fill_value="extrapolate")
        x_guess = t_compact_chl + 30  # Estimation initiale de l'intersection
        x_intersection = fsolve(find_intersection, x_guess)[0]
        y_intersection = f1(x_intersection)
        print(f"Rupture V compact= {x_intersection:.4f} ans")
        solution_V_compact.append(x_intersection)
        
        ### ---- Non linéarité : cas linéaire----
        if NL :
            V_NL = V_rd(As_E_linear/0.11,fywd,d*10**-3,cot_theta)
            f1 = interp1d(t, V_NL, kind='linear', fill_value="extrapolate")
            f2 = interp1d(t, V*np.ones(len(t)), kind='linear', fill_value="extrapolate")
            x_guess = min(t_fissure_chl,t_compact_chl) + 30  # Estimation initiale de l'intersection
            x_intersection = fsolve(find_intersection, x_guess)[0]
            y_intersection = f1(x_intersection)
            print(f"Rupture V linéaire {x_intersection:.4f} ans")
            solution_V_NL.append(x_intersection)
            
        
        
        
    
    
    # Plot évolutions temps avant rupture en fonction de l'enrobage
    plt.plot(enrobage,solution_M,label="Rupture par flexion")
    plt.plot(enrobage,solution_V,label="Rupture par effort tranchant")
    plt.xlabel("Enrobage [mm]")
    plt.ylabel("Temps [ans]")
    plt.title("Temps avant la rupture en fonction de l'enrobage")
    plt.legend()
    plt.grid()
    #plt.savefig("figures/temps_avant_rupture.pdf",bbox_inches='tight')
    plt.show()
    
    # Rupture moment en fonction enrobage béton fissuré et compact
    plt.plot(enrobage,solution_M,label="Rupture par flexion béton fissuré")
    plt.plot(enrobage,solution_M_compact,label="Rupture par flexion béton compact")
    plt.xlabel("Enrobage [mm]")
    plt.ylabel("Temps [ans]")
    plt.title("Temps avant la rupture en fonction de l'enrobage")
    plt.legend()
    plt.grid()
    #plt.savefig("figures/M_fissure_compact.pdf",bbox_inches='tight')
    plt.show()
    
    # Rupture effort tranchant en fonction enrobage béton fissuré et compact
    plt.plot(enrobage,solution_V,label="Rupture par effort tranchant béton fissuré")
    plt.plot(enrobage,solution_V_compact,label="Rupture par effort tranchant béton compact")
    plt.xlabel("Enrobage [mm]")
    plt.ylabel("Temps [ans]")
    plt.title("Temps avant la rupture en fonction de l'enrobage")
    plt.legend()
    plt.grid()
    #plt.savefig("figures/V_fissure_compact.pdf",bbox_inches='tight')
    plt.show()
    
    if NL : 
        # Différence temps de rupture
        difference_M = np.array(solution_M) - np.array(solution_M_NL)
        difference_V = np.array(solution_V) - np.array(solution_V_NL)
        difference_t_cl = np.array(temps_cl_NL) - np.array(temps_cl_L)
        
        
        # Non linéarité impact enrobage - temps de rupture 
        fig, axs = plt.subplots(1, 3, figsize=(18, 6))

        # Plot for Initiation corrosion by chlorides
        axs[0].plot(enrobage, difference_t_cl, label="Initiation corrosion par chlorure")
        axs[0].set_xlabel("Enrobage [mm]")
        axs[0].set_ylabel("Différence temps [ans]")
        axs[0].legend()
        axs[0].grid()
        axs[0].set_title("Initiation corrosion par chlorure")

        # Plot for Rupture by flexion
        axs[1].plot(enrobage, difference_M, label="Rupture par flexion")
        axs[1].set_xlabel("Enrobage [mm]")
        axs[1].set_ylabel("Différence temps rupture [ans]")
        axs[1].legend()
        axs[1].grid()
        axs[1].set_title("Rupture par flexion")

        # Plot for Rupture by shear force
        axs[2].plot(enrobage, difference_V, label="Rupture par effort tranchant")
        axs[2].set_xlabel("Enrobage [mm]")
        axs[2].set_ylabel("Différence temps rupture [ans]")
        axs[2].legend()
        axs[2].grid()
        axs[2].set_title("Rupture par effort tranchant")

        plt.tight_layout()
        #plt.savefig("figures/comparison_NL.pdf", bbox_inches='tight')
        plt.show()
        
    

    
    
    
    
    
   
    
    
    
