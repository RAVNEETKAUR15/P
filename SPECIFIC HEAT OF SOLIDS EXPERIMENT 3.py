# importing libraries
import numpy as np
import matplotlib.pyplot as plt

def delta(x):
    if x==0:
        return f*Na
    else:
        return 0

# for Debye integral
def simp(quad,T,a,b): #simpson
    n=1000
    h = float((b-a)/n)
    result = (1/3)*(quad(a,T)+quad(b,T))
    for i in range(1,n,2):
        result+= 4*(quad(a+i*h,T))
    for j in range(2,n-1,2):
        result+=2*(quad(a+j*h,T))
    result*=h/3
    return result

def trap(quad,a,b): #trapezoidal
    n=1000
    h=float((b-a)/n)
    result = 0.5*(quad(a))+0.5*(quad(b))
    for i in range(1,n):
        result+=(quad(a+i*h))
    result *= h
    return result

def Int_Cv(x):
    I=(x**3)/(np.exp(x)-1)
    return I

def Ein_G(v): #Einstein dos
    Glist=[]
    for i in range(len(v)):
        G=f*Na*delta(v[i]-1)
        Glist.append(G/(f*Na))
    return Glist

def Debye_G(v):   #Debye dos
    Glist=[]
    for i in range(len(v)):
        if v[i]<=1:
            G=(9*Na*v[i]**2)
            Glist.append(G)
        else:
            Glist.append(0)
    return Glist

# func for specific heat
def Duo_Cv(y): 
    Cvlist=[]
    for k in range(len(y)):
        Cv=3*R
        Cvlist.append(Cv/(3*R))
    return Cvlist

def Ein_Cv(y):
    Cvlist=[]
    for k in range(len(y)):
        x=1/y[k]
        Cv=((x**2)*np.exp(x)/(np.exp(x)-1)**2)*3*R
        Cvlist.append(Cv/(3*R))
    return Cvlist

def Debye_Cv(y):
    Cvlist=[]
    for k in range(len(y)):
        x=1/y[k]
        c1=(-3*x)/(np.exp(x)-1)
        c2=(12/(x**3))
        c3=trap(Int_Cv, 0.0000000001, x)
        Cv=(c1+c2*c3)*3*R
        Cvlist.append(Cv/(3*R))
    return Cvlist

def Ein_U(T):   #Einstein internal energy
    a=3*Na*h*V_e
    b=(h*V_e)/k
    U=(a)/(np.exp(b/T)-1)

    return U/(3*R)

def Int_U(v,T): #average energy
    a=h*v**3
    b=(h*v)/k
    I=(a)/(np.exp(b/T)-1)
    return I

def Debye_U(T):  #Debye internal energy
    
    c=9*Na/V_d**3

    U=(c)*(simp(Int_U, T,0.01, V_d))
        
    return U

#constants

R=8.31       #J/Kmol
k=1.38e-23   #J/K
h=6.626e-34  #Js
f=3
Na=6.022e23
V_e=6.4e12   # hertz
V_d=1.54e13   # hertz


ve=np.arange(0,2*V_e,V_e/100)
vd=np.arange(0,2*V_d,V_d/100)

xe=ve/V_e
xd=vd/V_d


theta1=1
theta2=h*V_e/k
theta3=h*V_d/k

n=1000

T1=np.linspace(0,2*theta1,n+1)
T2=np.linspace(0,2*theta2,n+1)
T3=np.linspace(0,2*theta3,n+1)

x1=T1/theta1
x2=T2/theta2
x3=T3/theta3

delT=T3[3]-T3[2]
dU_Einc=[]
dU_Debye=[]
 
def central_d(f,x):
    hstep=theta2/n                   
    central = ((f(x + hstep) - f(x-hstep))/(2*hstep))                 
    return central

def forward_d(f,x):
    hstep=theta2/n                  
    forward = ((f(x+hstep) - f(x))/(hstep))                 
    return forward

def backward_d(f,x):
    hstep=theta2/n                  
    backward = ((f(x) - f(x-hstep))/(hstep))                 
    return backward


for i in range(len(T2)):
    Ud=Debye_U(T3[i])-Debye_U(T3[i-1])
    dU_Debye.append(Ud/(delT*3*R))   
    dUe=central_d(Ein_U, T2[i])
    dU_Einc.append(dUe)

z1=dU_Einc
z2=dU_Debye

y1=Ein_Cv(x1)
y2=Debye_Cv(x2)
y3=Duo_Cv(x3)

w1=Ein_G(xd)
w2=Debye_G(xe)

xe=xe*V_e/V_d

#(a) : Density of states
plt.plot(xe,w1,label="Einstein",c='purple')
plt.plot(xd,w2,label="Debye",c='red')
plt.title(" Density of states")
plt.xlabel("V/Vx")
plt.ylabel("G(V/Vx)")
plt.legend()
plt.grid()
plt.show()

#(b): Cv/3R vs T/0
plt.plot(x1,y1,'.',label="Einstein")  #b
plt.plot(x2,y2,"*",label="Debye",c='deeppink')  #c
plt.plot(x3,y3,"^",label="Dulong")  #a
plt.title("Distribution functions")
plt.xlabel("T/$ \Theta $")
plt.ylabel("Cv/3R")
plt.grid()
plt.legend()
plt.show()

#d
plt.plot(x2[1:],z1[1:],"*",c='yellow',label="central derivated Einstein")
plt.plot(x2[1:],y1[1:],"--",c='blue',label="Einstein")
plt.plot(x3[1:],z2[1:],"*",c='cyan',label="central derivated Debye")
plt.plot(x3[1:],y2[1:],"--",c='deeppink',label="Debye")
plt.title("Distribution functions")
plt.xlabel("T")
plt.ylabel("Cv/3R")
plt.grid()
plt.legend()
plt.show()


