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
         
def As_reduction (Ns,d0,tcarbo,tchlor,t,vcorr1,vcorr2) :
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

def find_ti_carbonatation (K,e,w) :
    """ Fonction qui trouve temps d'initiation de la carbonatation
    K : coefficient de diffusion de la carbonatation (mm/sqrt(ans))
    x : pénétration de la carbonatation (mm)
    e : enrobage (mm)
    w : largeur de la fissure (mm)
    """
    # Béton compact
    t_compact = (e/K) ** 2 # [ans]
    # Béton fissuré
    t_fissure = (e/(50 * np.sqrt(w)))**4 # [ans]
    # Plot
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
    

if __name__ == "__main__" :
    K = 7 #[mm/sqrt(ans)]
    t = np.linspace(0,200,201) # [ans]
    w = 0.3 # [mm]
    t_compact1, t_fissure1 = find_ti_carbonatation(K,50,0.3)
    
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
    