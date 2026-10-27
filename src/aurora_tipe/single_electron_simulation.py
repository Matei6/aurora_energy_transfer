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

#O
O_E_ioni = 13.6   #en eV
O_E_1exc = 9.43
O_E_2exc = 11.63

#O2
O2_E_ioni = 12.1   #en eV
O2_E_1exc = 9.7

#N2
N2_E_ioni =  15.6
N2_E_1exc = 13





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





#expression du libre parcours moyen
def libre_parcours_moy(h):
    compo = compo_atm(h)
    d = (r_O * compo[0] + r_O2 * compo[1] + r_N2 * compo[2]) * 10**-2
    return (k*temperature(h) * 10**-3) / (np.pi * Pression(h) * d**2)





c_ioni = []           #compteur d'atome ionisé
c_exc = 0           #compteur d'atome excité




#input : hauteur de collision h en km, énergie Ei de l'électron en eV

def atome_collision (h):
    #renvoie avec quelle espèce l'électron est rentré en contact
    # 1 : O      2 : O2       3 : N2
    P = rd.randint(0, 100)
    if h < 100 :
        if P > 22 :
            return (3)
        else : 
            return (2)
    else :
        compo = compo_atm(h)
        if P < compo[0]:
            return(1)
        elif P < compo[0] + compo[1] :
            return(2)
        elif P < compo[0] + compo[1] +compo[2] :
            return(3)



def ener_coll (Ei, h):
    #renvoie la nature de la collision, la nouvelle énergie de l'électron
    atome = atome_collision (h)
    global c_ioni
    global c_exc
    if atome == 1 :
        if Ei > O_E_ioni : 
            c_ioni.append(h)
            Ef = Ei - O_E_ioni
        elif O_E_1exc < Ei :
            c_exc += 1
            Ef = Ei - O_E_1exc
        else : Ef = Ei
            
            
            
            
    elif atome == 2 :

        if Ei > O2_E_ioni : 
            c_ioni.append(h)
            Ef = Ei - O2_E_ioni
        elif O2_E_1exc < Ei :
            c_exc += 1
            Ef = Ei - O2_E_1exc
        else : Ef = Ei            
            
    elif atome == 3 :

        if Ei > N2_E_ioni : 
            c_ioni.append(h)
            Ef = Ei - N2_E_ioni
        elif N2_E_1exc < Ei :
            c_exc += 1
            Ef = Ei - N2_E_1exc
        else : Ef = Ei
        
    return (Ef)





#simule le parcours d'un seul électron
def parcours_electron (h, E):
    compteur_coll =[]
    ha = h
    E_elec = E
    Tabl = [[ha], [E]]
    while ha > 80 and E_elec > 9.43 :
    
        lambd = libre_parcours_moy(ha)
        dist_parc =  min (lambd, 5)         #pour qu'il y ait moins de calculs, on contrôle la distance parcouru à chaque itération
        print (lambd, dist_parc)
        P_collision = 1 - np.exp(-1)                  #la probabilité est toujours la même
        r = rd.random ()
        
        if r < P_collision :
            compteur_coll.append((ha,1))
            E_elec = ener_coll(E_elec, ha)
            
        ha = ha - dist_parc
        Tabl[0].append(ha)
        Tabl[1].append(E_elec)
            
    return compteur_coll, Tabl




#nombre de molécules ionisés

H= [ h for h in range (250, 80, -1)]
ioni_rate = [0 for i in H]

    
    
def taux_ioni(h, E):
    parcours_electron(h, E)
    c_io = c_ioni.copy()
    for i in range(len( H)):
        while c_io[0] > H[i] and len(c_io)!=1:
            ioni_rate[i] +=1
            c_io.pop(0)
    return ioni_rate



def main():
    a, b = parcours_electron(250, 5000)
    plt.plot(b[1], b[0])
    plt.ylabel('hauteur en km')
    plt.xlabel("énergie de l'électron (en eV)")
    plt.title("énergie d'un électron")
    plt.legend()
    plt.grid()
    plt.show()


if __name__ == "__main__":
    main()