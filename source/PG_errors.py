import math

import numpy as np

import PG_maths as pgm
import PG_globalVars as gv
import matplotlib.pyplot as plt

#Here the objective is to make error propagation
# WORK IN PROGRESS

order = 30

file = open("gnv.txt", "rt")

T = []
X = []
Y = []
Z = []
cpt = 0
for line in file:
    if cpt < 2000:
        tokens = line.split(" ")
        assert(tokens[2] == "E")
        assert(tokens[16] == "00000000\n")
        X.append(float(tokens[3]))
        Y.append(float(tokens[4]))
        Z.append(float(tokens[5]))
        T.append(float(tokens[0]))
        cpt += 1

print(len(X))

Ax, Ay, Az = pgm.SG_filt_cartesian([list(X), list(Y), list(Z), list(T)],2)
X, Y, Z = pgm.SG_filt_cartesian([list(X), list(Y), list(Z), list(T)],0)

R = []
THETA = []
PHI = []
for i in range(len(X)):
    R.append(np.sqrt(X[i]**2 + Y[i]**2 + Z[i]**2))
    THETA.append(np.arccos(Z[i]/R[i]))
    PHI.append(math.atan2(Y[i], X[i]))

plt.plot(T,X)
plt.plot(T,Y)
plt.plot(T,Z)
plt.figure()
plt.plot(T,R)
plt.figure()
plt.plot(T,THETA)
plt.plot(T,PHI)

# print(Ax[17])
# Ar = []
# for i in range(len(Ax)):
#     Ar.append((Ax[i]*X[i]+Ay[i]*Y[i]+Az[i]*Z[i])/R[i])
# print(Ar)
# Arth = []
# for i in range(len(Ax)):
#     Arth.append(-gv.GM_e/(R[i]**2))
# print(Arth)
#
# Aspherique = [
#                 + np.sin(THETA[0]) * np.cos(PHI[0]) * Ax[0]
#                 + np.sin(THETA[0]) * np.sin(PHI[0]) * Ay[0]
#                 + np.cos(THETA[0]) * Az[0],
#                 + np.cos(THETA[0]) * np.cos(PHI[0]) * Ax[0]
#                 + np.cos(THETA[0]) * np.sin(PHI[0]) * Ay[0]
#                 - np.sin(THETA[0]) * Az[0],
#                 - np.sin(PHI[0]) * Ax[0]
#                 + np.cos(PHI[0]) * Ay[0]
# ]
#
#
# print("vecteur sphérique 0")
# print(Aspherique)
# print('norme')
# print(np.sqrt(Aspherique[0]**2 + Aspherique[1] **2 + Aspherique[2]**2))
# print("valeur de Ax avec seulement Ar")
# print(Ar[0]*np.sin(THETA[0]) * np.cos(THETA[0]))
# print(pgm.vector_from_spherical_to_cartesian(Aspherique,THETA[0], PHI[0]))
# print(Ax[0], Ay[0], Az[0])
# print(pgm.vector_from_spherical_to_cartesian(Aspherique,THETA[0], PHI[0])[0]-Ax[17], pgm.vector_from_spherical_to_cartesian(Aspherique,THETA[0], PHI[0])[1]-Ay[17])
# print(X[0],Y[0],Z[0])
# print(THETA[0]*180/np.pi,PHI[0]*180/np.pi,R[0])
#
# coeffs = pgm.linear_regression_solution([list(Ax), list(Ay), list(Az)], [list(R), list(THETA), list(PHI), list(T)], 2)
# print([x/coeffs[0][0] for x in coeffs[0]])
# print(coeffs)


