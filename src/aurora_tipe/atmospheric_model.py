# -*- coding: utf-8 -*-
"""
Created on Thu May 15 11:21:21 2025

@author: matei
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
import random as rd


R = 8.314       #constantes des gaz parfaits
g=9.817         #accéleration de la pesanteur
k= 1.38e-23         #constante de Boltzmann

#rayon atomique
r_O = 66e-12
r_O2 = 146e-12
r_N2 = 158e-12

#composition de l'atmosphère
def compo_atm(h):

    if h < 100 : 
        prop_N2 = 78
        prop_O2 = 22
        prop_O = 0
    else :
        prop_N2 = -0.38*h + 116
        prop_O2 = -0.14*h + 36
        prop_O = 0.57*h - 57
    return (prop_O,prop_O2,prop_N2)



#profil de la temperature en K
def temperature(h):
    if 51 <= h < 71:
        return 270.7 - 2.8 * (h - 51)
    elif 71 <= h < 85:
        return 214.7 - 2.0 * (h - 71)
    elif 85 <= h < 91:
        return 186.9
    elif 91 <= h <= 250:
        return 186.9 + 3.67 * (h - 91)

#masse molaire effective de l'air en kg/mol
def massemol_eff (h):
    compo = compo_atm(h)
    return compo[0]*16e-5 + compo[1]*32e-5 + compo[2]*28e-5 
    

#facteur correctif pour g
def l(h):
    return ( 6378 / (6378 + h) )**2


#pression à l'altitude h en Pa
def Pression (h):
    
    def f(z):                                                           #pour pouvoir intégrer grace à quad
        return (l(z) * g * massemol_eff(z)) / (R * temperature(z))
    
    if 51 <= h < 71 :
        res, err = quad(f, h, 51)
        return 70.5 * np.exp(res * 10**3)

    elif 71 <= h < 85:
        res, err = quad(f, h, 71)
        return 4.48 * np.exp(res * 10**3)
        
    
    elif 85 <= h < 91 :
        res, err = quad(f, h, 85)
        return 0.446 * np.exp(res * 10**3)
        
    elif 91 <= h <= 250:
        res, err = quad(f, h, 91)
        return 0.154 * np.exp(res * 10**3)





#expression du libre parcours moyen en km
def libre_parcours_moy(h):
    compo = compo_atm(h)
    d = (r_O * compo[0] + r_O2 * compo[1] + r_N2 * compo[2]) * 10**-2
    return (k*temperature(h) * 10**-4) / (np.pi * Pression(h) * d**2)




#simule le parcours d'un seul électron
def parcours_electron (h, E):
    compteur_coll =[]
    ha = h
    for i in range (10000) :
        lambd = libre_parcours_moy(ha)
        dist_parc =  min (lambd, 5)         #pour qu'il y ait moins de calculs, on contrôle la distance parcouru à chaque itération
        print (lambd, dist_parc)
        P_collision = 1 - np.exp(-1)                  #la probabilité est toujours la même
        r = rd.random ()
        if r < P_collision :
            compteur_coll.append((ha,1))
        ha = ha - dist_parc
        if ha < 80 :
            return compteur_coll, ha
            
    return compteur_coll, ha



def main():
    parcours_electron(250, 50)


    H = np.linspace(250, 51,400)

    L = []
    for i in H:
        L.append(libre_parcours_moy(i))
    

    
    plt.plot(L,H)
    plt.semilogx(L,H)
    plt.ylabel('hauteur en km')
    plt.xlabel('libre parcours moyen (en m)')
    plt.title("Libre parcours moyen")
    plt.legend()
    plt.grid()
    plt.show()

    libre = [L,H]


if __name__ == "__main__":
    main()

"""
#tracé de la pression atmosphérique et du libre parcours moyen
H = np.linspace(250, 51, 400)
P = []
for i in H:
    P.append(Pression(i))
       
L = []
for i in H:
    L.append(libre_parcours_moy(i))
    
plt.plot(P,H)
plt.semilogx(P,H)
plt.ylabel('hauteur en km')
plt.xlabel('pression en Pa (échelle log)')
plt.title('représentation pression atmosphérique')
plt.legend()
plt.grid()
plt.show()

plt.plot(L,H)
plt.semilogx(L,H)
plt.ylabel('hauteur en km')
plt.xlabel('libre parcours moyen(en km)')
plt.title("représentation du libre parcours moyen d'un électron")
plt.legend()
plt.grid()
plt.show()
"""

"""

#simule le parcours d'un seul électron
def parcours_electron (h, E):
    c_coll =[[],[0]]                          #tableau compteur de collisions
    while h >= 100 and E >=  9.43 :
        lambd = libre_parcours_moy(h)
        dist_parc = 1               #pour qu'il y ait moins de calculs, on contrôle la distance parcouru à chaque itération
        P_collision = 1 - np.exp(-dist_parc/lambd)                  #la probabilité est toujours la même
        r = rd.random ()
        if r < P_collision :
            c_coll[0].append (h)
            c_coll[1].append(c_coll[1][-1] + 1)
    return c_coll

"""



"""
H =np.linspace(250, 51, 250)
compo = [[],[],[]]

for i in range(len(H)) :
    compo[0].append(compo_atm(H[i])[0])
    compo[1].append(compo_atm(H[i])[1])
    compo[2].append(compo_atm(H[i])[2])




plt.plot(H, compo[1], label="O2")
plt.plot(H, compo[0], label = 'O')
plt.plot(H, compo[2], label = 'N2')
plt.xlabel('hauteur en km')
plt.ylabel('composition atmosphérique')
plt.title('représentation des espèces présentes')
plt.legend()
plt.grid()
plt.show()

"""


