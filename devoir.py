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
         
def As_reduction (Ns,d0,t0,t,vcorr) :
    """ Fonction qui calcule la réduction de la section d'armature
        Ns : nombre de barres d'armature
        d0 : diamètre des barres d'armature (mm)
        t0 : temps de début de corrosion (ans)
        t : temps (ans)
        vcorr : vitesse de corrosion (mm/ans)
    """
    if t <= t0 :
        return Ns * np.pi * d0**2 / 4
    else :
        return Ns * np.pi * (d0 - 2*vcorr*(t-t0))**2 / 4

def find_ti_carbonatation (K,e,w) :
    """ Fonction qui trouve temps d'initiation de la carbonatation
    K : coefficient de diffusion de la carbonatation (mm/sqrt(ans))
    x : pénétration de la carbonatation (mm)
    e : enrobage (mm)
    w : largeur de la fissure (mm)
    """
    # Béton compact
    t_compact = (e/K) ** 2
    # Béton fissuré
    t_fissure = (e/(50 * np.sqrt(w)))**4
    # Plot
    plt.plot(t_compact,x_compact(K,t_compact),'or',markersize = 5,label="Temps d'initiation béton compacte")
    plt.hlines(e,0,t_compact,'r',linestyles='dashed')
    plt.vlines(t_compact,0,e,'r',linestyles='dashed')
    plt.plot(t_fissure,x_fissure(w,t_fissure),'ob',markersize = 5,label="Temps d'initiation béton fissuré")
    plt.hlines(e,0,t_fissure,'b',linestyles='dashed')
    plt.vlines(t_fissure,0,e,'b',linestyles='dashed')
    plt.plot(t,x_compact(K,t),label="Béton compacte")
    plt.plot(t,x_fissure(w,t),label="Béton fissuré")
    plt.xlabel("Temps [ans]")
    plt.ylabel("Pénétration de la carbonatation[mm]")
    plt.title("Pénétration de la carbonatation en fonction du temps")
    plt.grid()
    plt.legend()
    plt.show()
    
    return t_compact, t_fissure
    
def find_ti_corrosion (e,Dce,Clim,Cs,w,Sm0) :
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
    plt.plot(t,C(x,Cs,Dce,t),label="Béton compacte")
    plt.plot(t,C(x,Cs,D_fissure(Dce,w,Sm0),t),label="Béton fissuré")
    plt.plot(t,Clim*np.ones(len(t)),'--',label="Limite de concentration du ciment")
    plt.plot(t_compact,C(x,Cs,Dce,t_compact),'or',markersize = 5,label="Temps d'initiation béton compacte")
    plt.hlines(Clim,0,t_compact,'r',linestyles='dashed')
    plt.vlines(t_compact,0,Clim,'r',linestyles='dashed')
    plt.plot(t_fissure,C(x,Cs,D_fissure(Dce,w,Sm0),t_fissure),'ob',markersize = 5,label="Temps d'initiation béton fissuré")
    plt.hlines(Clim,0,t_fissure,'b',linestyles='dashed')
    plt.vlines(t_fissure,0,Clim,'b',linestyles='dashed')
    plt.xlabel("Temps [ans]")
    plt.ylabel("Concentration en chlorures [%]")
    plt.title("Concentration en chlorures en fonction du temps")
    plt.grid()
    plt.legend()
    plt.show()
    
    
    return t_compact, t_fissure
    

if __name__ == "__main__" :
    K1 = 9 #[mm/sqrt(ans)]
    K2 = 7.5 #[mm/sqrt(ans)]
    K3 = 4 #[mm/sqrt(ans)]
    t = np.linspace(0,200,100) # [ans]
    w = 0.3 # [mm]
    t_compact1, t_fissure1 = find_ti_carbonatation(K1,50,0.3)
    
    Clim = 0.4 # [%]
    Clim = Clim  * 0.15 # [%] car dans le béton C40 y a 15% de ciment
    Dce = 5e-13 # [m^2/s]
    Dce = Dce * 1e6 # [mm^2/s]
    Dce = Dce * 365*24*60*60 # [mm^2/ans]
    Cs = 0.15 # [%]
    x = 50 # [mm], pos[ition de l'armature (enrobage)
    Sm0 = 300 # [mm]
    t_compact2, t_fissure2 = find_ti_corrosion(50,Dce,Clim,Cs,w,Sm0)
    
    print("Temps d'initiation de la corrosion : ",min(t_compact1,t_compact2,t_fissure1,t_fissure2)," ans")
    
    vcorr = 2 # [microm/ans] si seulement carbonatation
    vcorr = 40 # [mm/ans] si carbonatation et corrosion
    
    