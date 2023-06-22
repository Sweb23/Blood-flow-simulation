"""version 1.3 : ajout du préconditionneur diag(A) selon la méthode
du gradient conjugué préconditionné"""

import numpy as np
import matplotlib.pyplot as mat

def diago(M):
    """Prend en entrée une matrice colonne
Renvoie la matrice diagonale ayant pour coefficients ceux de la matrice entrée"""
    n = len(M)
    D = np.zeros((n,n))
    for i in range(n):
        D[i][i] = M[i][0]
    return D

def mat_to_float(x):
"""effectue l'assimilation entre une matrice de taille 1,1 et un réel"""
    if x.shape != (1,1):
        return "pas la bonne taille"
    else:
        return x[0][0]

def conjgradPrec(A, b, x, C = np.array([[1,0,0],[0,1,0],[0,0,1]])):
    r = b - A @ x
    inv_C = np.linalg.inv(C)
    z = inv_C @ r
    d = z
    for k in range(len(b)):
        alpha = mat_to_float(np.transpose(z)@z)/mat_to_float(np.transpose(d)@A@d)
        x = x + alpha*d
        r_new = r - alpha*A@d #r(k+1)
        z_new = inv_C@r #z(k+1)
        beta = mat_to_float(np.transpose(z_new)@r_new)/mat_to_float(np.transpose(z)@r)
        d = z_new + beta*d
        z = z_new
        r = r_new
    return x

def P_new(P_old,C_old,C,S):
    """Calcule P(t) en considérant S comme connu. P_old, C_old, C sont des matrices colonnes"""
    A = np.zeros((6,6))
    A = -dt*(((S*G) + np.transpose(S*G)))
    A = A + diago(C - np.transpose(np.array([np.sum(A, axis=0)])))
    inv_A = np.linalg.inv(A)
    P = conjgradPrec(A,C_old*P_old,P_old, np.diag(np.diag(A)))
    #P = inv_A @ (C_old*P_old) version sans gradient conjugué
    return P

def verification_S(S,P,C,k):
    """P_new trouve P en fonction de S, et cette fonction détermine des nouvelles valeurs de S consistantes en calculant P en
fonction de S puis S en fonction de P etc."""
    done = False #si done == True, on a fini.
    while not done:
        S_old = S
        P[:,[k]] = P_new(P[:,[k-1]],C[:,[k-1]],C[:,[k]],S)
        for i in range(6): #on détermine les nouvelles valeurs des valves
            for j in range(6):
                S[i][j] = int(P[i][k] > P[j][k]) #définition de S
        done = np.array_equal(S,S_old)
    return P[:,[k]], S

def C_new(t,CVS,CVD):
    """calcule les nouvelles valeurs de C pour les ventricules"""
    tauS = 0.0025
    tauD = 0.0075
    if t % T < TS:
        e = (1-np.exp(-(t%T)/tauS))/(1-np.exp(-TS/tauS))
        CV = CVD*(CVS/CVD)**e
    else:
        e = (1-np.exp(-((t%T)-TS)/tauD))/(1-np.exp(-(T-TS)/tauD))
        CV = CVS*(CVD/CVS)**e
    return CV

T = 0.75 #période d'un battement de coeur en s
dt = 0.01*T
CLVS = 0.00003 #compliance du ventricule gauche lors de la systole (L/mmHg)
CLVD = 0.0146 #compliance du ventricule gauche lors de la diastole (L/mmHg)
CRVS = 0.0002 #compliance du ventricule droit lors de la systole (L/mmHg)
CRVD = 0.0365 #compliance du ventricule droit lors de la diastole (L/mmHg)
TS = 0.33 #durée de la systole (s)
TD = T - TS #durée de la diastole (s)
R = np.array([[np.inf,0.01,np.inf,np.inf,np.inf,np.inf],[np.inf,np.inf,17.5,np.inf,np.inf,np.inf],[np.inf,17.5,np.inf,0.01,np.inf,np.inf],[np.inf,np.inf,np.inf,np.inf,0.01,np.inf],[np.inf,np.inf,np.inf,np.inf,np.inf,1.79],[0.01,np.inf,np.inf,np.inf,1.79,np.inf]])
P0 = np.array([[5],[80],[2],[2],[8],[5]])
C0 = np.array([[C_new(0,CLVS,CLVD)],[0.00175],[1.75],[C_new(0,CRVS,CRVD)],[0.00412],[0.08]])
G = np.array([[1/R[j][i] for i in range(6)] for j in range(6)])

def simulation_chambres(C_ini,P_ini,R,n):
    S = [[int(P_ini[i]>P_ini[j]) for j in range(6)] for i in range(6)]
    t_plot = [0]*n
    P_plot = np.zeros((6,n))
    C_plot = np.zeros((6,n))
    P_plot[:,[0]] = P_ini[:,[0]]
    C_plot[:,[0]] = C_ini[:,[0]]
    for k in range(1,n):
        t = k*dt
        t_plot[k] = t
        C_plot[1][k],C_plot[2][k],C_plot[4][k],C_plot[5][k] = C_plot[1][k-1],C_plot[2][k-1],C_plot[4][k-1],C_plot[5][k-1]
        C_plot[0][k], C_plot[3][k] = C_new(t,CLVS,CLVD), C_new(t,CRVS,CRVD)
        P_plot[:,[k]], S = verification_S(S,P_plot,C_plot,k)
    return P_plot,C_plot,t_plot """P_plot est la matrice de terme général P_plot[i][j] = pression dans la chambre i à l'étape j
(idem pour C_plot)"""

def graphe(A,B):
    """Renvoie le graphe de la ligne j de A en fonction de B, où A est une matrice carrée et B une liste (utile pour la
simulation des 6 chambres)"""
    longueur = len(B)
    Out = [0]*longueur
    j = int(input("Quelle chambre ? "))
    for i in range(longueur):
        Out[i] = A[j][i]
    mat.plot(B,Out)
    mat.show()

"""indices :
0 : ventricule gauche
1 : artères systémiques
2 : veines systémiques
3 : ventricule droit
4 : artères pulmonaires
5 : veines pulmonaires
"""