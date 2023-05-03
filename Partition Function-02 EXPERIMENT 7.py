# -*- coding: utf-8 -*-
"""
Created on Fri Mar 24 08:53:04 2023

@author: LENOVO
"""

# Experiment-7
#Simranjeet Chhabra
#2020PHY1028
#Partion function-02

#importing libraries
import numpy as np
from scipy.integrate import quad
from scipy.stats import linregress
import matplotlib.pyplot as plt
import pandas as pd

#constants
h=6.626e-34     #Js
k=1.38e-23      #J/K
m=9.1e-31      # mass of electron (in kg)
N=6.022e23      #N_A

V = np.linspace(20e-3, 50e-3, 50)
T = np.linspace(150, 450, 50)

#solving integeration
def partition_fnc(V, T):
    z_v = []
    
    for v in V:
        z_t = []
        
        for t in T:
            
            integrand =lambda n: (n**2)*np.exp(-((h*n)**2)/(8*m*k*t*v**(2/3)))
            integral = (np.pi/2)*quad(integrand, 0, 10**11)[0]
            
            z_t.append(integral)
        z_v.append(z_t)
        
    return np.array(z_v)
   
#forward method
def derivative(y, x):
    
    dy_dx = [(y[i+1]-y[i])/(x[i+1]-x[i]) for i in range(len(x)-1)]
        
    return np.array(dy_dx)
    

Z = partition_fnc(V, T)
lnZ = np.log(Z)

#pressure
# taking transpose of the matrix for different volumes
P = np.array([N*k*T[i]*derivative(lnZ.T[i], V) for i in range(len(T)-1)]) 

# Average energy(U/N)
E = k*(T[:len(T)-1]**2)*np.array([derivative(lnZ[i], T) 
                                  for i in range(len(T)-1)])[0] 

#Entropy
S = [N*E[i]/T[i] + N*k*(lnZ[i]-np.log(N)+1) for i in range(len(T)-1)]

plt.figure(figsize = (20, 15))

plt.subplot(2, 2, 1)
for i in [1, 10, 45, -6]: # taking any 4 values from matrix
    plt.plot(T, lnZ[i], label = f'V = {np.round(V[i], 3)} $m^3$')
plt.grid() 
plt.legend()
plt.xlabel('T')
plt.ylabel('ln Z')
plt.title('ln Z vs T')
    
plt.subplot(2, 2, 2)
for i in [1, 10, 45, -6]:
    plt.plot(np.log(T), lnZ[i], label = f'V = {np.round(V[i], 3)} $m^3$')
plt.grid() 
plt.legend()
plt.xlabel('ln T')
plt.ylabel('ln Z')
plt.title('ln Z vs ln T')

plt.subplot(2, 2, 3)
for i in [4, 20, 44, -9]:
    plt.plot(V, lnZ.T[i], label = f'T = {np.round(T[i], 3)} $k$')

plt.legend()
plt.xlabel('V')
plt.ylabel('ln Z')
plt.title('ln Z vs V')
plt.grid() 

plt.subplot(2, 2, 4)
for i in [4, 20, 44, -9]:
    plt.plot(np.log(V), lnZ.T[i], label = f'T = {np.round(T[i], 3)} $K$')

plt.legend()
plt.xlabel('ln V')
plt.ylabel('ln Z')
plt.title('ln Z vs ln V')
plt.grid() 

plt.show()

plt.figure(figsize = (20, 7))

plt.subplot(1, 3, 1)
for i in [1, 12, 35, -2]:
    plt.plot(T[:len(T)-1], P.T[i], label = f'V = {np.round(V[i], 3)} $m^3$')
 
plt.legend()
plt.xlabel('T')
plt.ylabel('P')
plt.title('Pressure vs Temperature')
plt.grid() 

plt.subplot(1, 3, 2)
for i in [1, 12, 35, -2]:
    plt.plot(V[:len(T)-1], P[i], label = f'T = {np.round(T[i], 3)} K')
 
plt.legend()
plt.xlabel('V')
plt.ylabel('P')
plt.title('Pressure vs Volume')
plt.grid() 

plt.subplot(1, 3, 3)
plt.plot(T[:len(T)-1], E, label = f'V = {V[0]} $m^3$')
plt.xlabel('T')
plt.ylabel('<E>')
plt.legend()
plt.title('<E> vs T')
plt.grid() 

plt.show()

plt.figure(figsize = (15, 7))

plt.subplot(1, 2, 1)
for i in [1, 12, 35, -2]:
    plt.plot(T[:len(T)-1],np.array(S).T[i],label=f'V ={np.round(V[i],3)}$m^3$')
plt.grid()  
plt.legend()
plt.xlabel('T')
plt.ylabel('S')
plt.title('Entropy vs Temperature')

plt.subplot(1, 2, 2)
for i in [1, 12, 35, -2]:
    plt.plot(V, np.array(S)[i], label = f'T = {np.round(T[i], 3)} K')
 
plt.legend()
plt.xlabel('V')
plt.ylabel('S')
plt.title('Entropy vs Volume')
plt.grid() 
plt.show()

nv = 10 # no of rows for volume
nt = 6 #no of columns for temp

table = pd.DataFrame(lnZ[:nv,:nt], index = [f'V{i+1}' for i in range(nv)])
table.columns = [f'T{i+1}' for i in range(nt)]
print("Matrix for Z :")
print(table)
print("................................................................")
slope = linregress(T[:len(T)-1], E)[0]
print('\nSlope of <E> vs T plot:', slope,'J/K')

#specific heat
C_v = derivative(E, T[:len(T)-1]) 
print('Calculated value of Cv:', C_v[-1],'J/K')
print('Theoretical value of Cv:', 1.5*8.314/(6.022*10**(23)),'J/K')
print(".................................................................")
print('\nEnergy fluctuations:', k*(T[:len(T)-2]**2)*C_v)