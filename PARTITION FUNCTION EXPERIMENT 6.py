# -*- coding: utf-8 -*-
"""
Created on Sun Feb 26 16:16:15 2023

@author: LENOVO
"""

#Simran jeet Chhabra
#2020PHY1028
# Pract-06
#partition function

#importing libraries
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress

# boltzman const
k = 8.61733 * 10**(-5)  # ev K^(-1)

#partition function
def f(e, g,N):
    def Z(g, e, T):
        
        p_fnc = [ g[i] * np.exp(-e[i] / (k*T)) for i in range(len(e)) ]
        
        return np.sum(np.array(p_fnc))
    
    cT = e[-1]/k  #charaterstic temp 
    print("charaterstic temp : ",cT)    
    T = np.linspace(0.001, cT*100, 10000) 
    z = [Z(g, e, i) for i in T]
    b = list(T).index(min(list(T), key = lambda x:abs(x-(cT/4))))
    
    plt.figure(figsize = (14,8))
    plt.subplot(1,2,1)
    plt.plot(T[:b], z[:b],label = 'Low Temprature',c='red')
    plt.legend(loc='best',fontsize='20')
    plt.xlabel('Temperature',fontsize='20')
    plt.ylabel('Partition Function(Z)',fontsize='20')
    plt.grid()
   
    plt.subplot(1,2,2)
    plt.plot(T[b:], z[b:],label = 'High Temprature',c='green')
    plt.legend(loc='best',fontsize='20')
    plt.xlabel('Temperature',fontsize='20')
    plt.ylabel('Partition Function(Z)',fontsize='20')
    plt.suptitle('Partition Function VS Temperature',fontsize='20')
    plt.grid()
    plt.show()
    
# Fractional polpulation
   
    plt.figure(figsize = (14,8))
    plt.subplot(1,2,1)
    
    for i in range(len(e)):
        fp = g[i]*np.exp(-e[i]/(k*T))/z
        
        plt.plot(T[:b], fp[:b], label = f'e={e[i]}, g={g[i]}')
       
    plt.plot(T[:b], np.ones(len(T[:b]))*0.5, ls='dotted',c='deeppink')
    plt.legend(loc='best',fontsize='20')
    plt.xlabel('T',fontsize='20') 
    plt.ylabel(r'$N_i/N$',fontsize='20')
    plt.grid()
    plt.title(f'low temperatures for {len(e)}-level system',fontsize='20')
    
    plt.subplot(1,2,2) 
    
    for i in range(len(e)):
        fp = g[i]*np.exp(-e[i]/(k*T))/z
        
        plt.plot(T[b:], fp[b:], label = f'e={e[i]}, g={g[i]}')
        
    plt.plot(T[b:],np.ones(len(T[b:]))*(1/(len(e))), ls='dotted',c='deeppink')
    plt.xlabel('T',fontsize='20') 
    plt.ylabel(r'$N_i/N$',fontsize='20')
    plt.grid()
    plt.legend(loc='best',fontsize='20')
    plt.title(f'high temperatures for {len(e)}-level system',fontsize='20')    
    plt.show()

#Internal energy        
    def fnc(g, e, T):
        f = [g[i]*e[i]*np.exp(-e[i]/(k*T)) for i in range(len(e))]
        
        return np.sum(np.array(f))
    
    u = [fnc(g, e, i) for i in T]  #u/n
    u = N*np.array(u)/z
    
    plt.figure(figsize = (14,8))
    plt.subplot(1,2,1)
    plt.plot(T[:b], u[:b],label = 'Low Temprature',c='purple')
    plt.xlabel('T',fontsize='20')
    plt.ylabel('U/N',fontsize='20')
    plt.legend(loc='best',fontsize='20')
    plt.grid()
    
    plt.subplot(1,2,2)
    plt.plot(T[b:], u[b:],label = 'High Temprature',c='darkblue')
    plt.xlabel('T',fontsize='20')
    plt.ylabel('U/N',fontsize='20')
    plt.legend(loc='best',fontsize='20')
    plt.suptitle('Average energy vs temperature',fontsize='20')
    plt.grid()   
    plt.show()
 
 #Entropy   
    #s= (N*k*np.log(np.array(z)/N)) + (u/(N*T)) +N*k 
    s = N*k*np.log(np.array(z)) + u/T
    
    plt.figure(figsize =  (14,8))
    plt.subplot(1,2,1)
    plt.plot(T[:b], s[:b],label = 'low Temprature',c='brown')
    plt.xlabel('T',fontsize='20')
    plt.ylabel('S/N',fontsize='20')
    plt.legend(loc='best',fontsize='20')
    plt.grid()
    
    plt.subplot(1,2,2)
    plt.plot(T[b:], s[b:],label = 'High Temprature',c='deeppink')
    plt.xlabel('T',fontsize='20')
    plt.ylabel('S/N',fontsize='14')
    plt.suptitle('Entropy vs temperature',fontsize='20')
    plt.legend(loc='best',fontsize='20')
    plt.grid()
    plt.show()
 
#free energy    
    fe = -N*k*T*np.log(z)
    
    plt.figure(figsize = (14,8))
    plt.subplot(1,2,1)
    plt.plot(T[:b], fe[:b],label='Low Temp',c='olive')
    plt.xlabel('T',fontsize='20')
    plt.ylabel('F/N',fontsize='20')
    plt.legend(loc='best',fontsize='20')
    plt.grid()
    
    plt.subplot(1,2,2)
    plt.plot(T[b:], fe[b:],label='High Temp',c='grey')
    plt.xlabel('T',fontsize='20')
    plt.ylabel('F/N',fontsize='20')
    plt.suptitle('Free energy vs temperature',fontsize='20')
    plt.legend(loc='best',fontsize='20')
    plt.grid()    
    plt.show()
    
    print('\nslope of free energy vs temp graph: ',-linregress(T[b:]
                                                               , fe[b:])[0])
    print('Entropy at Tmax: ', max(s))
    print('\nintercept of free energy vs temp graph: ',linregress(T[b:],
                                                                  fe[b:])[1])
    print('Internal energy at Tmax: ', max(u))
    
#2-level system 
print("2-level system : ")  
f(e = [0, 1], g = [1,1],N=10)
print(".............................")
#3-level system   
print("3-level system : ")  
f(e = [0, 1,2], g = [1,1,1],N=10)
