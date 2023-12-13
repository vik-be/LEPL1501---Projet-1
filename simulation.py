#####################################
## Simulation physique groupe 1136 ##
#####################################

import numpy as np
import matplotlib.pyplot as plt

"""Simulation informatique de l'évolution de l'angle de la barge sur le temps
---- Formules du modèle
dv = Somme des forces /m
dx = x0 + v.dt
---- Principes du modèle
Energie (inertie) constante
"""
# ----- Constantes initiales ------ #
g = 9.81

# ----Paramètres de la barge ------ #
m1 = 4 #masse de la plateforme
m2 = 1 #masse de la grue seule
mtot = 5    #masse totale [kg]
l = 0.6  #Longueur de la plateforme

h1 = 0.08    #hauteur de la plateforme [m]
h2=0.1  

# ----- Variables ----- #
dist = 1  #distance totale parcourue [m]
m_charge = 0.2  #masse de la charge
d = 15 #coefficient d'amortissement






# ----- Initialisation (tableaux) ------ #

end = 20  #Déplacement de la masse de la grue
step = 0.01
t = np.arange(0, end, step)
w = np.empty_like(t)
theta = np.empty_like(t)


# ------- Fonctions utiles ------- #

def get_h_im(mtot, l):
    return mtot  / (1000*(l**2))

def get_inertie() :
    return 7


def get_max_angle():
    """Détermine l'angle maximal d'inclinaison
    
    Returns:
        int : l'angle [degré]
    """
    a1 = np.arctan(2*(h1 - h_im)/l)
    a2 = np.arctan(2*h_im/l)

    if a1 < a2:
        return a1
    
    return a2

def get_couple_red(theta,m_charge,dist):
    """Calcule le couple de redressement de la grue

    Args:
        theta (float): l'angle d'inclinaison en radians

    Returns:
        couple (float): le couple de redressement in N.m.
    """

    # Diviser l'équation en plusieurs load (pour éviter de se tromper)
    load1 = (h1 - h_im) * np.sin(theta) * m1
    load2 = (h1 - h_im + h2) * np.sin(theta) * (m2- m_charge)
    load3 = dist * m_charge
    load4 = m1 + (m2-m_charge)+ m_charge

    # calculer le load total
    xg = (load1 + load2 + load3) / load4

    # Calculer le centre de gravité
    xc = (l ** 2) * np.tan(theta) / (12 * h_im)

    # Calculer la force d'archimède
    f = 1000 * 9.81 * (l ** 2) * h_im

    # Calculer le couple de redressement
    return f * (xc - xg)




def get_couple_a():
    """retourne le couple déstablisiateur Ca

    Args:
        g (float): constante de graviation de la Terre
        m_charge (float): masse de la charge déplacée [kg]
        d (float): distance déplacée [m]

    Returns:
        float: Ca
    """
    return g*m_charge*dist


def get_angle_deplacement(dist,m_charge):
    """Donne l'angle d'inclinaison de l'a grue quand elle déplace une masse de plusieurs cm.
    Elle retourne une solution aproximée mais néanmoins assez précise pour lui faire confiance.

    Returns:
        float: l'angle [rad]
    """
    # Interval initial
    sup = 0.0  # Limite inférieure
    inf = 2.0  # Limite supérieure

    # Critère d'arrêt
    arret = 1e-6

    # Itération
    while (inf - sup) > arret:
        res = (sup + inf) / 2
        if get_couple_red(res,m_charge,dist) == 0:
            break
        elif get_couple_red(res,m_charge,dist) * get_couple_red(sup,m_charge,dist) < 0:
            inf = res
        else:
            sup = res
    return res



# ----- Deuxième initialisation de variables grâce aux fonctions utiles ---- #

h_im = get_h_im(mtot,l) #hauteur immergée [m]
theta_0 = get_angle_deplacement(0.001, m_charge)
inertie = get_inertie()



# ------- Fonctions pour les graphes -------- #

def couple_red_theta(theta) :
    """
    Calcul le couple de redressement en fonction d'un angle theta en radians

    Args:
        theta (int): Angle d'inclinaison [rad]

    Returns:
        int : le couple de redressement
    """
    f = 1000*9.81*(l**2)*h_im   #rho * g * volume immergé (longueur**2 *hauteur_immergée)

    h1 = h_im + l*np.tan(theta) / 2
    h2 = h_im - l*np.tan(theta) / 2
    hc = (h1**2 + h1*h2 + h2**2)  / (3*(h1+h2))
    lc = l*(h1 + 2*h2) / (3*(h1+h2))
    yc = (-l/2) + lc
    zc = (-h_im) + hc
    ycp = yc*np.cos(theta) - zc*np.sin(theta)
    total_load = f * ycp
    return total_load




def couple_a_theta(theta) :
    """
    Calcule le couple destabilisateur en fonction d'un angle theta

    Args:
        theta (int): angle d'inclinaison à l'instant t [rad]

    Returns:
        int: couple destabilisateur
    """
    ca = (m_charge * g * np.cos(theta))
    return ca



def get_angle():
    """
    Remplit les tableaux w et theta des bonnes valeurs en fonction de l'angle de départ
    
    """

    #Conditions initiales
    theta[0] = theta_0 
    dt = step
    for i in range (len(t)-1):
        loadw = (-d*w[i])+ couple_red_theta(theta[i]) + couple_a_theta(theta[i])
        w[i+1] = w[i] + (loadw / inertie) * dt
        theta[i + 1] = theta[i] + w[i]*dt
        
        




# ----- Création des tableaux ----- #

angle_max = get_max_angle()
        
def graphique_masses() :
    angle_max = get_max_angle()
    t2 = np.arange(0, 2, 0.1)
    for i in range(2,7,2) :
        m = np.empty_like(t2)
        for j in range(1,len(t2)):
            dist = t2[j]
            m[j] = get_angle_deplacement(dist,i)
        plt.plot(t2,m, label = f" m = {i}")
        
        plt.legend()
    plt.axhline(y=angle_max, color = "red", label = "angle maximum possible", linestyle = '--')
    plt.grid(True)
    plt.show()

def graphique_incl():
    """
    Fait deux sortes de graphiques :
        -un graphique qui montre l'inclinaison de la grue en fonction du temps pour un dépalcement 
        -un graphique qui montre l'inclinaison maximum de la grue pour un déplacement d pour chaque changement de paramètre
    
    """
    angle_max = get_max_angle()
    get_angle()

    plt.figure(1)
    plt.plot(t,theta, label = "angle d'inclinaison")
    plt.axhline(y=angle_max, color = "red", label = "angle maximum possible", linestyle = '--')
    plt.xlabel("temps [s]")
    plt.ylabel("angle d'inclinaison [rad]")
    plt.title("Angle d'inclinaison en fonction du temps")
    plt.legend()
    plt.grid(True)
    plt.show()

def graphiques_angle_vitesse_accélération():
    """Produit 3 graphiques avec l'angle en fonction du temps, la vitesse angulaire et l'accélération angulaire
       Unités utilisées dans le graphiques : radians, secondes, radians/seconde, radians/secondes^2"""
 
    
    # Sous-graphique : 1. Angle [rad] en fonction du temps [s]
    
 
    plt.subplot(3,1,1)
    plt.plot(t,theta, label="Angle")
    plt.legend()
    plt.xlim(0,end)
    
    plt.plot([0,end], [angle_max,angle_max], '--r', label='submersion') 
    plt.plot([0,end], [-angle_max,-angle_max], '--r')                  
    plt.ylabel('Angle \n (rad)')
    
 
    # Sous-graphique : 2. Vitesse angulaire [rad/s] en fonction du temps [s]
    
    plt.subplot(3,1,2)
    plt.plot(t,theta,'g', label="Vitesse angulaire")
    plt.legend()
    plt.xlim(0,end)
    plt.grid(True)
    plt.ylabel('Vitesse angulaire \n (rad/s)')
    
    # Sous-graphique : 3. Accélération angulaire[rad/s**2] en fonction du temps [s]
    
    plt.subplot(3,1,3)
    plt.plot(t,w, label="Accélération angulaire")
    plt.legend()
    plt.grid(True)
    plt.xlim(0,end)
    plt.ylabel('Accélération angulaire \n (rad/s^2)')
    
    # Paramètres généraux

    plt.xlabel('Temps (s)')
    plt.grid(True)
    plt.show()

def diagramme_phase():

    plt.plot(theta, w, label="Diagramme de phase")
    plt.xlabel("Angle d'inclinaison [rad]")
    plt.ylabel("Vitesse angulaire [rad/s]")
    plt.title("Diagramme de phase")
    plt.legend()
    plt.grid(True)
    plt.show()  

if __name__ == "__main__":
    get_angle()
    graphique_incl()
    graphique_masses()
    graphiques_angle_vitesse_accélération()
    diagramme_phase()


