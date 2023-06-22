import numpy as np
import matplotlib.pyplot as mat

#initialisation
T = 0.0125 #durée d'un battement en minutes
TS = 0.0050 #durée de la systole
Rs = 17.86 #résistance systémique (base 17,86)
Csa = 0.001575 #compliance systémique (valeur à ajuster, valeur de base 0.00175)
Rmi = 0.01 #résistance de la valve mitrale (entre l'oreillette gauche et le ventricule gauche)
RAo = 0.6 #résistance de la valve aortique (entre le ventricule gauche et les artères)
dt = 0.01 * T
PLA = 5 #pression de l'oreillette gauche (juste avant le ventricule gauche)
n = 2000
CLVS = 0.00003 #compliance du ventricule gauche lors de la systole (L/mmHg)
CLVD = 0.0146 #compliance du ventricule gauche lors de la diastole (L/mmHg)
PLV0 = 5
Psa0 = 80
Smi0 = int(PLA > PLV0)
SAo0 = int(PLV0 > Psa0)


def C_new(t,CVS,CVD): #calcule les nouvelles valeurs de C pour les ventricules
    tauS = 0.0025
    tauD = 0.0075
    if t % T < TS:
        e = (1-np.exp(-(t%T)/tauS))/(1-np.exp(-TS/tauS))
        CV = CVD*(CVS/CVD)**e
    else:
        e = (1-np.exp(-((t%T)-TS)/tauD))/(1-np.exp(-(T-TS)/tauD))
        CV = CVS*(CVD/CVS)**e
    return CV


def P_new(PLV_old,Psa_old,CLV_old,CLV,Smi,SAo):
    C11 = CLV + dt*((Smi/Rmi)+(SAo/RAo))
    C12 = -dt*(SAo/RAo)
    C22 = Csa + dt*((SAo/RAo)+(1/Rs)) #cf p48
    B1 = CLV_old*PLV_old + dt*(Smi/Rmi)*PLA
    B2 = Csa*Psa_old
    D = C11*C22-(C12**2) #déterminant du système
    PLV = (B1*C22 - B2*C12)/D
    Psa = (B2*C11 - B1*C12)/D
    return PLV, Psa

def S_new(Smi,SAo,PLV_old,Psa_old,CLV_old,CLV):
    done = False
    while not done:
        Smi_prime = Smi
        SAo_prime = SAo
        PLV,Psa = P_new(PLV_old,Psa_old,CLV_old,CLV,Smi,SAo)
        Smi = int(PLA > PLV)
        SAo = int(PLV>Psa)
        done = (Smi == Smi_prime) and (SAo == SAo_prime)
    return PLV,Psa

def systeme():
    PLV = PLV0
    Psa = Psa0
    Smi = Smi0
    SAo = SAo0
    CLV = C_new(0,CLVS,CLVD)
    t_plot = [0]*n
    PLV_plot = [0]*n
    Psa_plot = [0]*n
    CLV_plot = [0]*n
    PLV_plot[0] = PLV
    Psa_plot[0] = Psa
    CLV_plot[0] = CLV
    for k in range(1,n):
        t = k*dt
        PLV_old = PLV
        Psa_old = Psa
        CLV_old = CLV
        CLV = C_new(t,CLVS,CLVD)
        PLV,Psa = S_new(Smi,SAo,PLV_old,Psa_old,CLV_old,CLV)
        PLV_plot[k] = PLV
        Psa_plot[k] = Psa
        CLV_plot[k] = CLV
        t_plot[k] = t
    return PLV_plot,Psa_plot,CLV_plot,t_plot

def graphe():
    PLV_plot,Psa_plot,CLV_plot, t_plot = systeme()
    j = int(input("0 pour PLV, 1 pour Psa, 2 pour CLV : "))
    if j == 0 :
        mat.plot(t_plot,PLV_plot)
        mat.xlabel("Temps")
        mat.ylabel("PLV")
        mat.title("Pression du ventricule gauche en fonction du temps\nModèle à 2 chambres")
        mat.show()
    elif j ==1 :
        mat.plot(t_plot,Psa_plot)
        mat.xlabel("Temps")
        mat.ylabel("Psa")
        mat.title("Pression artiérielle sytémique en fonction du temps\nModèle à 2 chambres")
        mat.show()
    else:
        mat.plot(t_plot,CLV_plot)
        mat.xlabel("Temps")
        mat.ylabel("CLV")
        mat.title("Compliance du ventricule gauche en fonction du temps\nModèle à 2 chambres")
        mat.show()