# -*- coding: utf-8 -*-
"""
Created on Sat Feb  4 11:30:56 2023

@author: LENOVO
"""

#Simran jeet Chhabra
#2020PHY1028
# Pract-04
#The Laws of Radiations
#Rayleigh-jeans and Planck's(Energy Density)
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          
# importing libraries
import numpy as np
import matplotlib.pyplot as plt

# Density of states(both)
def G(x):
    glist=[]
    for i in range(len(x)):
        g=np.pi*(x[i]**2)
        glist.append(g)
    return glist

#The Rayleigh function
def Rayleigh(x): 
    rlist=[]
    for i in range(len(x)):
        r=(x[i]**2)
        rlist.append(r)
    return rlist
    
lo=1e-10        #m (given)
h=6.62e-34      #Js
k=1.38e-23      #J/K
c=3e8           #m/s2
e = 1.6 *1e-19   # C
nu_o=c/(2*lo)   #hertz
nu_p=np.logspace(10,30,50) #logspace distributes the equally spaced points
visrange=np.linspace(1e14,10e14,20) #for visible range
xv=visrange/nu_o

#plotting

#a(i)
plt.plot(nu_p,G(nu_p),'o-',label="Density of states in complete range",c='purple')
plt.title("Density of states VS frequency")
plt.ylabel("Density of states")
plt.xlabel(" Frequency")
plt.xscale("log")
plt.legend()
plt.grid()
plt.show()


#a(ii)
plt.plot(xv,G(xv),label="Density of states in visible range",c='deeppink')
plt.title("Density of states VS frequency")
plt.ylabel("Density of states")
plt.xlabel("Frequency")
plt.legend()
plt.grid()
plt.show()

#a(iii)
plt.plot(nu_p,G(nu_p),c='red')
plt.title("Density of states VS frequency")
plt.ylabel("Density of states in log scale")
plt.xlabel("Frequency in log scale")
plt.xscale("log")
plt.yscale("log")
plt.grid()
plt.show()

#Rayleigh jeans law,energy of states
e_dim=np.linspace(0,12,100) #It was given that the energy varies from 0 to 12 
#(b)
fr=np.array(Rayleigh(e_dim))
plt.plot(e_dim,fr,c='brown')
plt.title("Frj(x) VS x")
plt.ylabel("Frj(x)")
plt.xlabel("x")
plt.legend()
plt.grid()
plt.show()

#b(1)
#for temp
def constant(T):
    estar=k*T
    lstar=(h*c)/estar
    cons=(8*np.pi*estar)/(lstar**3)
    return cons

#Temp=[1200,1500,1800]  #K
Temp=[4000,5000,6000,8000,10000]
rayleigh=[]

for T in Temp:
    estar=k*T
    freq=np.linspace(0,12,100)  #the dimensional frequency is converted to dimensionless energy
    c1=constant(T) 
    c2=np.array(Rayleigh(e_dim)) 
    listr=c1*c2
    rayleigh.append(c1*c2)
    plt.plot(freq,listr,label=f"For T = {T}K")
    plt.axvspan(1.25,1.55,facecolor='r',alpha=0.1)
plt.title("Rayleigh Jean's Energy Spectrum Density vs freq in complete range")
plt.ylabel("Rayleigh Jean's Spectrum ")
plt.xlabel("frequency")
plt.legend()
plt.grid()
plt.show()

# Plancks radiation law ,energy of states
    
def Planck(x):
    plist=[]
    for i in range(len(x)):
        p=(x[i]**3)/(np.exp(x[i])-1)
        plist.append(p)
    return plist
# (c)
freq_dim=np.linspace(0,12,100)
fp=np.array(Planck(freq_dim))
plt.plot(freq_dim,fp,c='green')
plt.title("Fp(x) vs x")
plt.ylabel("Fp(x)")
plt.xlabel("x")
plt.grid()
plt.show()

planck=[]

for T in Temp:
    estar=k*T
    e_dim=np.linspace(0,12,100) 
    c1=constant(T)
    nu=e_dim*estar/h
    c3=np.array(Planck(e_dim))
    planck.append(c1*c3)   
    plt.plot(nu,c1*c3,label=f"For T = {T}K")   
    plt.axvspan(1,1.55,facecolor='r',alpha=0.1)
plt.title("Planck's Energy of states VS Frequency")
plt.ylabel("Planck's Energy states")
plt.xlabel("Frequency")
plt.legend()
plt.grid()
plt.show()

# Visible Range frequency
nu_vs = np.linspace(4*1e+14, 8*1e+14, 100)
def planck_energy(nu, T):
    return h*nu/(np.exp(h*nu/(k*T)) -1) / e

E_planck = np.vectorize(planck_energy)
for temp in Temp:
    E_arr = E_planck(nu_vs, temp)
    plt.plot(nu_vs, E_arr, label='T = '+str(temp)+' K')
plt.xlabel(r'$\nu$ (Hz)')
plt.xscale('log')
plt.ylabel(r'Avg. energy, $\epsilon_P(\nu)$ (eV)')
plt.title('Planck\'s average energy')
plt.grid()
plt.legend()
plt.show()