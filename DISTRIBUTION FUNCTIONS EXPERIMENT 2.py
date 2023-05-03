# -*- coding: utf-8 -*-
"""
Created on Fri Jan 13 08:54:11 2023

@author: LENOVO
"""

#Simran jeet Chhabra
#2020PHY1028
#Prac2:Plot the following functions at different temp and energy of the system.
'''
                (1) Maxwell-Boltzmann Distribution
                (2) Bose-Einstein Distribution
                (3) Fermi-Dirac Distribution
'''
# importing libraries
import numpy as np
import matplotlib.pyplot as plt

# Maxwell-Boltzmann Distribution 
def max_bolt(E,T):
    f=[]
    a=1 #constant
    for i in range(len(E)):
        e=(-E[i]/(k*T))  #exp(-x)
        f.append(np.exp(e))  
    return f

k_eV= 8.617333262*10**(-5)        # eV/K (boltzman const)
k=k_eV
k_j=1.38064852*10**(-23)          # J/K
t=11805                          # K                        
temp_in_eV=t*k_eV                 # eV
temp_in_j=t*k_j               # J
temp_conv=temp_in_eV
yl=5
E=np.linspace(-yl,yl,1000)
x=E*temp_conv
y1=max_bolt(E, t)

# Bose -Einstien Distribution
def bos_ein(E,T):
    b=[]    
    for i in range(len(E)):
        e=(E[i]/(k*T))-(U/(k*T))  
        b.append(1/(np.exp(e)-1))  #1/(exp(alpha+x) - 1)
    return b 
#q=0    
#U=q*temp_conv 
U=0
y2=bos_ein(E, t)

#Fermi Dirac
def fer_dir(E,T):
    f=[]
    for i in range(len(E)):
        e=(E[i]/(k*T))-(Ef/(k*T))  #1/(exp(x)*np.exp(alpha)+1)
        f.append(1/(np.exp(e)+1))
    return f
#Ef=q*temp_conv 
Ef=0 
y3=fer_dir(E, t)

#plotting
fig,axes=plt.subplots(1,3,figsize=(10,4))
axes[0].plot(E,y1,"red")
axes[0].set_title("Maxwell-Boltzmann Distribution")
axes[0].set_xlabel("$\epsilon/KT$")
axes[0].set_ylabel("F_mb($\epsilon$)")
axes[0].grid()

axes[1].plot(E,y2,"green")
axes[1].set_title(" Bose -Einstien Distribution")
axes[1].set_xlabel("$\epsilon/KT$")
axes[1].set_ylabel("F_be($\epsilon$)")
#axes[1].set_ylim(0,200)
#axes[1].set_xlim(0,5)
axes[1].grid()

axes[2].plot(E,y3,"orange")
axes[2].set_title(" Fermi -Dirac Distribution")
axes[2].set_xlabel("$\epsilon/KT$")
axes[2].set_ylabel("F_fd($\epsilon$)")
axes[2].grid()

fig.tight_layout()
plt.show()

#all three in one graph
plt.plot(x,y1,".",label="Maxwell-Boltzmann")
plt.plot(x,y2,".",label="Fermi-Dirac")
plt.plot(x,y3,".",label="Bose-Einstein ")
plt.grid()
plt.title("Distribution Graph")
plt.xlabel("$\epsilon/KT$")
plt.ylabel("F($\epsilon$)")
plt.ylim(0,1)
plt.legend()
plt.show() 

Tlist=[5000,1000,100,10]       #list of temp     
Tlist_Energy=k*np.array(Tlist)
Tmax=Tlist_Energy.max()

#Maxwell-Boltzmann as a function of energy
T1 = [500, 1000, 10000]
e = np.linspace(-0.1, 0.1, 200)

for T in T1:
    z0=max_bolt(e, T)
    plt.plot(e,z0,label=f"for Temp={T} K")
plt.xlabel("$\epsilon $")
plt.ylabel("$ F_{MB}(\epsilon) $")
plt.title('Maxwell Boltzmann as function of energy')
plt.legend()
plt.grid()
plt.show()

#fermi-dirac as a function of energy
for T in Tlist:    
    Ef=1
    z1=fer_dir(E, T)
    plt.plot(E,z1,label=f"for Temp={T} K")
plt.grid()
plt.title("Fermi-Dirac for $E_f$ = 1 eV")
plt.xlabel("$\epsilon $")
plt.ylabel("$ F_{FD}(\epsilon) $")
plt.legend(loc="best")
plt.show()    

#Bose-Einstein as a function of energy
for T in Tlist:
    U=1
    E=np.linspace(U,U+2*Tmax,1000)
    w2=bos_ein(E, T)
    plt.plot(x,w2,label=f"for Temp={T} K")
plt.grid()
plt.title("Bose-Einstein for $ \mu $ = 1 eV ")
plt.xlabel("$\epsilon $")
plt.ylabel("$ F_{BE}(\epsilon) $")
plt.ylim(0,5)
plt.legend()
plt.show() 


def f_mb(x, alpha):
    return np.exp(-x)

def f_be(x, alpha):
    return 1/(np.exp(x+alpha) -1)

def f_fd(x, alpha):
    return 1/(np.exp(x+alpha) +1)

def distribution_analysis(fn, alpha, E_system, ptitle):    
   
    x = np.linspace(- alpha-3, -alpha+3, 1000)
    fx = [fn(i, alpha) for i in x]

    fx_arr = [i for i in fx if i>0]
    x_arr = [ x[j] for j in range(len(x)) if fx[j] > 0]
    
   
# Temperature Variation of the probability function
    theta = E_system/k
    T_by_theta = [1/i for i in x_arr]
    T_arr = []
    y_arr = []
    for i in range(len(x_arr)):
        t = (1/x_arr[i])
        if t >= 0:
            T_arr.append(t)
            y_arr.append(fx_arr[i])
    plt.xlabel(r'$T/\theta$')
    plt.ylabel(r'$f(t)$')
    plt.title(r'Temperature Variation')
    plt.plot(T_arr, y_arr)
    plt.grid()
    plt.show()

distribution_analysis(f_mb, 0, 0.1, 'Maxwell')
distribution_analysis(f_be, 0, 0.1, 'Bose-Einstein')
distribution_analysis(f_fd, 0, 0.1, 'Fermi-Dirac')

