# -*- coding: utf-8 -*-
"""
Created on Tue Feb 14 21:12:18 2023

@author: LENOVO
"""
#Simran jeet Chhabra
#2020PHY1028
# Pract-05
#The Laws of Radiations
#Stefan-Bolztmann Law ( Radiant Flux )

#important libraries
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress

#Wiens Displacement law 
def weins(x):
    f_p=(x**3)/(np.exp(x)-1)
    return f_p

#constants
a1=0.00000000001 #lower limit
b1=100     #upper limit
x=np.linspace(a1,12,100)
xp=np.linspace(a1,b1,100)
y=weins(x).tolist()

#plotting
plt.plot(x,y,c='darkblue')
plt.text(2.78,1.4,"xp=2.78",c='black',fontsize='14')
plt.text(3.5,1.33,"x_mid=3.5",c='black',fontsize='12')
plt.vlines(x=2.78,ymin=0,ymax=1.4,colors='red',linestyle='dotted')
plt.hlines(1.4,xmin=0,xmax=2.78,colors='red',linestyle='dotted')
plt.hlines(1.334,xmin=0,xmax=3.5,colors='red',linestyle='dotted')
plt.vlines(3.5,ymin=0,ymax=1.334,colors='red',linestyle='dotted')
plt.title("f_p(x) vs x")
plt.ylabel("f_p(x)")
plt.xlabel("x")
plt.grid()
plt.show() 

#peak value
h=6.626e-34      #Js
k=1.38e-23      #J/K
c=2.99e8           #m/s2

def Weins_law(x,y):
   max_value=max(y)
   max_index=y.index(max_value)
   max_peak=x[max_index]
   b=(h*c)/(k*max_peak)
   return max_peak,b,max_value
x_p,b,e =Weins_law(x, y) #peak value and weins constant
print("----------------------------------")
print("\nThe value of peak x (i.e x_p) : ",x_p)
print("\nThe value of fp : ",e)
print("The value of Wein's constant b in meter Kelvin :",b)

#Stefan-Bolztmann Law 

def simp(quad,a,b):
    n=1000
    h = float((b-a)/n)
    result = (1/3)*(quad(a)+quad(b))
    for i in range(1,n,2):
        result+= 4*(quad(a+i*h))
    for j in range(2,n-1,2):
        result+=2*(quad(a+j*h))
    result*=h/3
    return result

def constant(T):
    estar=k*T
    lstar=(h*c)/estar
    cons=(8*np.pi*estar)/(lstar**3)
    return cons

Itrue=(np.pi**4)/15  #actual value of integeration
print("----------------------------------")
print('\n\N{greek small letter pi}\N{superscript four}/15 :',Itrue)

def Integration(weins,a1,b1):
    I=simp(weins, a1, b1)
    error=(abs(I-Itrue)/I)
    return I,error
I,err=Integration(weins,a1,b1)

print("\nIntegration Value : ",I)
print("\nError in Integrated value :",err,"\n\n\n")

def Energy_density(T):
    I,err=Integration(weins, a1, b1)
    u=constant(T)*I
    return u

Temp=np.arange(1,10000,500)  #K

def Radiant_flux(T):
    f=(c/4)*Energy_density(T)
    return f
rlist=[]
for T in Temp:
    r=Radiant_flux(T)
    rlist.append(r)

plt.plot(Temp,rlist,"o-",c='purple')
plt.title("Radiant Flux VS Temprature")
plt.ylabel("Radiant flux")
plt.xlabel("Temprature")
plt.grid()
plt.show()

plt.plot(np.log(Temp),np.log(rlist),'o',c="deeppink")
slope,intercept=np.polyfit(np.log(Temp),np.log(rlist),1)
plt.plot(np.log(Temp),slope*np.log(Temp)+intercept,c='green')
plt.title("lnF vs lnT")
plt.ylim(-17,20)
plt.xlim(0,10)
plt.ylabel("lnF")
plt.xlabel("lnT")
plt.grid()
plt.show()

slope=linregress(np.log(Temp),np.log(rlist))[0]
intercept=linregress(np.log(Temp),np.log(rlist))[1]
print("----------------------------------")
print("from plot of lnF vs lnT :")
print("slope :",slope)
print("intercept",intercept)

print("----------------------------------")
print("Stenfan-Bolztmann constant : ")
print("Numerical value from graph : ",np.exp(intercept))
sigma=2*(np.pi**5)*(k**4)/(15*(c**2)*(h**3))
print("Standard Value : ",sigma)

#divinding the graph in equal parts

def xm(x,weins):
    Itr=Itrue/2
    toll=0.001
    It=round(Itr,3)
    Ilist=[]
    for i in range(len(x)):
        I=simp(weins, a1, x[i])
        I=round(I,3)
        Ilist.append(abs(It-I))
    min_I=min(Ilist)
    max_index=Ilist.index(min_I)
    max_peak=x[max_index]
    b=(h*c)/(k*max_peak)
    return max_peak,b

t1,t2=xm(x,weins)
print("--------------------------")
print("xm:",t1)
print("The value of Wein's constant b in meter Kelvin :",t2)

#wavelength
print("\nwavelength (T=6000K) from the median value =", t2/5778)
print("wavelength (T=6000K) from the peak value =", b/5778)