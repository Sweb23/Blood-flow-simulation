import numpy as np
import matplotlib.pyplot as mat

#initialisation
T = 0.0125 #durée d'un battement en minutes
TS = 0.0050 #durée de la systole
Tmax = 0.0020 #point auquel l'écoulement est maximal dans l'artère aorte
Qmax = 28 #Q(Tmax)
Rs = 17.86 #résistance systémique (base 17,86)
Csa = 0.001575 #compliance systémique (valeur à ajuster, valeur de base 0.00175)
dt = 0.01 * T
n = 2000

def ecoulement(t):
    tc = t % T
    if tc < TS :
        if tc < Tmax:
            Q = Qmax *tc/Tmax

        else:
            Q = Qmax * (TS - tc)/(TS - Tmax)
    else:
        Q = 0
    return Q

def Psa_new(Psa_old,QAo):
    Psa  =(Psa_old+dt*QAo/Csa)/(1+dt/(Rs*Csa))
    return Psa

def arteres():
    t_plot = [0]*n
    P_plot = [0]*n
    Q_plot = [0]*n
    P = 80
    for k in range(1,n):
        t = k*dt
        Q = ecoulement(t)
        P = Psa_new(P,Q)
        t_plot[k] = t
        Q_plot[k] = Q
        P_plot[k] = P
    return P_plot,Q_plot,t_plot

def graphe():
    P_plot,Q_plot,t_plot = arteres()
    j = int(input("0 pour Pression, 1 pour ecoulement : "))
    if j == 0 :
        mat.plot(t_plot,P_plot)
        mat.xlabel("Temps")
        mat.ylabel("Psa")
        mat.title("Pression systèmique artiérielle en fonction du temps\nModèle à 1 chambre")
        mat.show()
    else:
        mat.plot(t_plot,Q_plot)
        mat.xlabel("Temps")
        mat.ylabel("Q")
        mat.title("Ecoulement artiériel en fonction du temps\nModèle à 1 chambre")
        mat.show()
