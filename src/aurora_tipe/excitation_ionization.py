# -*- coding: utf-8 -*-
"""
Created on Wed May 14 10:02:41 2025

@author: matei
"""
import random as rd


#dans ce cas, il y a collision entre l'électron et un atome de l'atmosphère. Il s'agit de déterminer quel est cet atome, et s'il y a excitation ou ionisation. 



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




c_ioni = []         #compteur d'atome ionisé
c_exc = 0           #compteur d'atome excité



def compo_atm(h):
    #à partir de 100 km jusqu'à 250 km
    prop_N2 = -0.38*h + 116
    prop_O2 = -0.14*h + 36
    prop_O = 0.57*h - 57
    return (prop_O,prop_O2,prop_N2)



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
        print(compo)
        print(P)
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
    
    
