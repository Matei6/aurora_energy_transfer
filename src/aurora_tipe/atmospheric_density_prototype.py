# -*- coding: utf-8 -*-
"""
Created on Tue Jan 14 10:22:53 2025

@author: matei
"""
import numpy as np
import matplotlib.pyplot as plt

#valeur des constantes
Ei = 5 #en keV, énergie initiale des particules   8e-16 J
vi = 4.2e7 #en m/s, vitesse initiale correspondante à Ei(energie considérée uniquement cinétique)

hi = 250 #en km, hauteur initiale des électrons injectés 
thêta = 65 #en degrés, l'experience se fait à angle constant
me = 9.1e-31 #en kg, masse d'un électron

#Note : on pourra extraire le densité p des différentes molécules (en fonction exponentielle par exemple )
#densité atmosphérique : d'abord modèle linéaire puis pour plus de complexité prendre un modèle exponentiel
#aussi séparer après, dans un modèle plus complexe les densités de O2, O et N2
#atmosphère = plein de tranches et à chaque tranche il y a une nouvelle chance de collision

def d_atm(h) :              #densité de l'atmosphère en km^-3
#h = hauteur en km
    d = h* (-6*10**12) + 1.5 * 10**15
    return d

#Première partie : pour une surface de 10 km² en découpant le tout en tranches de 1 km de hauteur (très grossier)
S = 100 #km²
def nb_molecules (h, S = 100, t = 1):
    #h est la hauteur, S est la surface de la tranche considérée et t est la hauteur de la tranche
    #tout est en km
    n = d_atm (h) * S * t
    return n

#pour l'instant, on considère qu'il n'y a que de l'oxygène et de l'azote dans l'atmosphère, en proportions 1/2
def rapport_surfaceocc_surfvide(h, S = 100, t = 1):
    s1 = 7.8e-23    #surface occupé par N2 en km
    s2 = 1.7e-23    #surface occupé par O en km
    r1 = (nb_molecules(h)/2) * s1 + (nb_molecules(h)/2) * s2
    R = r1/100
    return (R)

def main():
    heights = np.linspace(100, 250, 150)

    plt.plot(d_atm(heights), heights)
    plt.ylabel('hauteur en km')
    plt.xlabel('densité atmosphérique')
    plt.title('représentation linéaire de la densité atmosphérique')
    plt.grid()
    plt.show()

    plt.plot(nb_molecules(heights), heights)
    plt.ylabel('hauteur en km')
    plt.xlabel('nombre de molecules')
    plt.title("représentation du nb de molécules dans l'atmosphère")
    plt.grid()
    plt.show()

    plt.plot(rapport_surfaceocc_surfvide(heights), heights)
    plt.ylabel('hauteur en km')
    plt.xlabel('rapport surface occupée sur la surface totale')
    plt.title("Rapport en fonction de la hauteur")
    plt.grid()
    plt.show()


if __name__ == "__main__":
    main()
