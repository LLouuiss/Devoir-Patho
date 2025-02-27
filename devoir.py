import numpy as np
from scipy.special import erf
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

if __name__ == "__main__" :
    K1 = 9 #[mm/sqrt(ans)]
    K2 = 7.5 #[mm/sqrt(ans)]
    K3 = 4 #[mm/sqrt(ans)]
    t = np.linspace(0,600,100) # [ans]
    w = 0.3 # [mm]
    plt.plot(t,x_compact(K1,t),label="Béton compacte K = 9")
    plt.plot(t,x_compact(K2,t),label="Béton compacte K = 7,5")
    plt.plot(t,x_compact(K3,t),label="Béton compacte K = 4")
    plt.plot(t, x_fissure(w, t), '--', label="Béton fissuré")
    plt.xlabel("Temps [ans]")
    plt.ylabel("Pénétration de la carbonatation[mm]")
    plt.title("Pénétration de la carbonatation en fonction du temps")
    plt.grid()
    plt.legend()
    plt.show()
    
    Clim = 0.4 # [%]
    Dce = 5e-13 # [m^2/s]
    Dce = Dce * 1e6 # [mm^2/s]
    Dce = Dce * 365*24*60*60 # [mm^2/ans]
    Cs = 0.15 # [%]
    x = 50 # [mm], pos[ition de l'armature (enrobage)
    Sm0 = 300 # [mm]
    plt.plot(t,C(x,Cs,Dce,t),label="Béton compacte")
    plt.plot(t,C(x,Cs,D_fissure(Dce,w,Sm0),t),label="Béton fissuré")
    plt.plot(t,Clim*np.ones(len(t)),'--',label="Limite de concentration du ciment")
    plt.xlabel("Temps [ans]")
    plt.ylabel("Concentration en chlorures [%]")
    plt.title("Concentration en chlorures en fonction du temps")
    plt.grid()
    plt.legend()
    plt.show()
    
    