import argparse
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
import random as rd
from pathlib import Path


R = 8.314       #constantes des gaz parfaits
g=9.817         #accéleration de la pesanteur
k= 1.38e-23         #constante de Boltzmann

#rayon atomique
r_O = 66e-12
r_O2 = 146e-12
r_N2 = 158e-12

#O
O_E_ioni = 13.6   #en eV
O_E_exc = 9.43


#O2
O2_E_ioni = 12.1   #en eV
O2_E_exc = 9.7

#N2
N2_E_ioni =  15.6
N2_E_exc = 13





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
    return (k*temperature(h) * 10**-3) / (np.pi * Pression(h) * d**2)               





c_exc = 0           #compteur d'atome excité


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



def ener_coll (Ei, h, c_ioni):
    #renvoie la nature de la collision, la nouvelle énergie de l'électron
    atome = atome_collision (h)
    global c_exc
    if atome == 1 :
        if Ei > O_E_ioni : 
            c_ioni.append(h)
            Ef = Ei - O_E_ioni
        elif O_E_exc < Ei :
            c_exc += 1
            Ef = Ei - O_E_exc
        else : Ef = Ei
            
            
            
            
    elif atome == 2 :

        if Ei > O2_E_ioni : 
            c_ioni.append(h)
            Ef = Ei - O2_E_ioni
        elif O2_E_exc < Ei :
            c_exc += 1
            Ef = Ei - O2_E_exc
        else : Ef = Ei            
            
    elif atome == 3 :

        if Ei > N2_E_ioni : 
            c_ioni.append(h)
            Ef = Ei - N2_E_ioni
        elif N2_E_exc < Ei :
            c_exc += 1
            Ef = Ei - N2_E_exc
        else : Ef = Ei
        
    return (Ef, c_ioni)




#def de la fonction factorielle
def fact(n):
    if n == 0:
        return 1
    else:
        return n  * fact(n-1)
    
    
    
#def de la probabilité de Poisson
def P(k, mu):
    return np.exp(-mu) * mu**k/ fact(k)





#donne les hauteurs auquels il ya collision
def nb_coll(h, dist_parc):
    lambd = libre_parcours_moy(h)
    mu = dist_parc / lambd
    r=rd.random()
    hc = []
    if r < P(0,mu) :
        return hc
    elif r < P(0,mu) + P(1,mu) :                # une collision
        f = rd.random()              #pour déterminer hauteur de la collision    
        hc.append(h - f*dist_parc)
        return hc
    
    elif r < P(0,mu) + P(1,mu) +P(2,mu) :        # deux collisions
        for i in range (2):
            f = rd.random()                 
            hc.append(h - f*dist_parc)
        return hc
            
    elif r < P(0,mu) + P(1,mu) +P(2,mu)+P(3,mu) :        # trois collisions
        for i in range (3):
            f = rd.random()                 
            hc.append(h - f*dist_parc)
        return hc  
    
    elif r < P(0,mu) + P(1,mu) +P(2,mu)+P(3,mu)+ P(4,mu) :  # quatre collisions
        for i in range (4):
            f = rd.random()                 
            hc.append(h - f*dist_parc)
        return hc
    return hc







#simule le parcours d'un seul électron
def parcours_electron (h, E):
    compteur_coll =[]                       #nombre de collision
    ha = h
    E_elec = E  
    c_ioni = []                            
    Tabl = [[ha], [E]]           #Tabl[0] = hauteur, Tabl[1] = énergie à h
    
    while ha > 80 and E_elec > 9.43 :
        lambd = libre_parcours_moy(ha)
        dist_parc =  min (lambd, 5)         
        H_coll = nb_coll(ha, dist_parc)
        if len(H_coll) != 0 :
            H_coll.sort(reverse = True )
            while len(H_coll) != 0:
                compteur_coll.append(H_coll[0])
                E_elec, c_ioni = ener_coll(E_elec, H_coll[0], c_ioni)
                Tabl[0].append(H_coll[0])
                Tabl[1].append(E_elec)
                H_coll.pop(0)

            
        ha = ha - dist_parc

            
    return compteur_coll, Tabl, c_ioni



H= [ h for h in range (250, 80, -1)]


#calcul du dépot d'énergie en eV/km
def depot_ener(h, E) :
    a, Tabl, c_ioni = parcours_electron(h, E)
    n = len(Tabl[1])
    depot = [0 for i in H]

    for i in range(len(H)) :
            while Tabl[0][0] > H[i] and len(Tabl[0])!=1 :
                depot[i]+= Tabl[1][0] -Tabl[1][1]
                Tabl[0].pop(0)
                Tabl[1].pop(0)
    return depot
        


    

def moy_depot (h, E):
    t_tot = []
    for i in range (1000):
        t_tot.append(depot_ener(h, E))
    n = len(t_tot[0])
    t=[0 for i in range (n)]


    for i in range (n):
        for j in range(1000):
            t[i] += t_tot[j][i]
        t[i] = t[i]* (1/1000)
        
    return (t, t_tot)



#taux d'ionisation pour un seul électron
def taux_ioni(h, E):
    a, b, cio =  parcours_electron(h, E)
    c_io = cio.copy()
    ioni_rate = [0 for i in H]
    for i in range(len( H)):
        while c_io and c_io[0] >= H[i]:
            ioni_rate[i] +=1
            c_io.pop(0)
    return ioni_rate


#moyenne pour 1000 électrons
def moy1000 (h, E):
    t_tot = []
    for i in range (1000):
        t_tot.append(taux_ioni(h, E))
    n = len(t_tot[0])
    t=[0 for i in range (n)]


    for i in range (n):
        for j in range(1000):
            t[i] += t_tot[j][i]
        t[i] = t[i]* (1/1000)
        
    return (t, t_tot)
'''

a, b    = moy_depot(250, 5000)


plt.plot(a,H)
plt.semilogx(a,H)
plt.ylabel('hauteur en km')
plt.xlabel("dépot d'énergie ($eV.km^{-1}$)")
plt.title("dépot d'énergie")
plt.ylim(100, 245)
plt.legend()
plt.grid()

'''






def main(h=250, E=5000):
    t, t2 = moy1000(h, E)
    if any(value > 0 for value in t):
        plt.semilogx(t, H)
        plt.xlim(10**-4)
    else:
        plt.plot(t, H)
    plt.ylabel('hauteur en km')
    plt.xlabel("nombre de molécules ionisées ($km^{-1}$)")
    plt.title("Taux d'ionisation")
    plt.ylim(100, 200)
    plt.grid()
    output_path = Path(__file__).resolve().parents[2] / "figures" / "ionization_profile.png"
    plt.savefig(output_path)
    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--altitude", type=float, default=250)
    parser.add_argument("--energy", type=float, default=5000)
    arguments = parser.parse_args()
    main(arguments.altitude, arguments.energy)

'''
#moyenne pour 100 électrons
def moy2 (h,E):
    t1 = taux_ioni(250, 5000)
    t2  = taux_ioni(250, 5000) 
    t=[]
    n = len(t1)
    for i in range (n):
        t.append(0.5*(t1[i] + t2[i]))
    return (t,t1,t2)



t,t1, t2=moy(250, 5000)

plt.plot(t,H)
plt.ylabel('hauteur en km')
plt.xlabel("énergie de l'électron (en eV)")
plt.title("énergie d'un électron")
plt.legend()
plt.grid()
plt.show()



def moyT (h,E):

    t_tot = []
    t_tot.append(taux_ioni(250, 5000))
    t_tot.append(taux_ioni(250, 5000))
    t=[]
    n = len(t_tot[0])
    for i in range (n):
        t.append(0.5*(t_tot[0][i] + t_tot[1][i]))
    return (t, t_tot)

'''






'''
plt.plot(t,H)
plt.semilogx(t,H)
plt.ylabel('hauteur en km')
plt.xlabel("nombre d'atome ionisés (km-1)")
plt.title("Taux d'ionisation")
plt.ylim(100, 200)
plt.legend()
plt.grid()
plt.show()

'''

'''



a, b = parcours_electron(250, 5000)


plt.plot(b[1],b[0])
plt.ylabel('hauteur en km')
plt.xlabel("énergie de l'électron (en eV)")
plt.title("énergie d'un électron")
plt.legend()
plt.grid()
plt.show()'''