import numpy as np
import pandas as pd
import random as r
import matplotlib.pyplot as plt
from scipy.stats import binom

def Prob(Nt,Nc):
    no_tails=[]
    no_heads=[]
    possible_outcomes=[]
    frequency_outcomes=[]
    probability=[]
    sample_space=[]
    com_fre_heads=[]
    com_prob_tails=[]
    com_fre_tails=[]
    com_prob_heads=[]
    N_list=[]
    bin_dis=[]
    fluctuation=[]
    
    for j in range(Nt):
        outcomes=[]
        for i in range(Nc):   
            t=r.choice(["H","T"])
            outcomes.append(t)
            t=outcomes.count("T")
            h=outcomes.count("H")
        sample_space.append(outcomes)
        no_tails.append(t)
        no_heads.append(h)
    
    for k in range(Nc+1):
        possible_outcomes.append(k)
        out=no_heads.count(k)
        frequency_outcomes.append(out)
        probability.append(out/Nt)
        bin_dis.append(binom.pmf(k, Nc, 1/2))

    for l in range(Nt):
        N=Nc*(l+1)
        N_list.append(N)
        if l==0:
            com_fre_heads.append(no_heads[l])
            com_fre_tails.append(no_tails[l])

        else:
            com_fre_heads.append((no_heads[l]+(com_fre_heads[-1])))
            com_fre_tails.append((no_tails[l]+(com_fre_tails[-1])))

    for m in range(len(com_fre_tails)): 
        com_fre_heads=np.array(com_fre_heads)
        com_fre_tails=np.array(com_fre_tails)
        N_list=np.array(N_list)
        cph=com_fre_heads[m]/N_list[m]
        cpt=com_fre_tails[m]/N_list[m]
        com_prob_heads.append(cph)
        com_prob_tails.append(cpt)
        
    # mean=sum(com_prob_heads)/len(com_prob_heads)
    com_prob_heads_sqr=np.array(com_prob_heads)**2
    mean2=sum(com_prob_heads_sqr)/len(com_prob_heads_sqr)
    for o in range(len(com_prob_heads)):
        fluctuation.append((abs(mean2-com_prob_heads[o]**2)))
    
    Not=range(Nt)
        
    data={"No. of trials ":range(1,Nt+1),"Outcome":sample_space,"No. of Heads":no_heads,"No. of Tails":no_tails}
    df=pd.DataFrame(data)
    df=df.set_index("No. of trials ")
    
    return df,probability,possible_outcomes,com_prob_tails,com_prob_heads,Not,bin_dis,fluctuation

Nc1=10
Nt2=1000
Nt3=10000
Nc3=10
Nt_list=[Nc1,10*Nc1,100*Nc1,1000*Nc1]
Nc_list=[3,4,5,6]

for Nt1 in Nt_list:    
    df1,probability1,possible_outcomes1,com_prob_tails1,com_prob_heads1,Not1,bin_dis1,fluctuation1=(Prob(Nt1, Nc1))
    #print("\nFor ",Nt1,"No. of trials and ",Nc1,"No. of Coins")
    #print(df1)
    plt.plot(possible_outcomes1,probability1,label=f"{Nt1} no. of Trials")
    plt.scatter(possible_outcomes1,probability1)
plt.plot(possible_outcomes1,bin_dis1,color="black",label="Binomial distribution")
plt.scatter(possible_outcomes1,bin_dis1,color="black")
plt.minorticks_on()
plt.grid(b=True,which="both",axis="both")
plt.title("Probability Vs No. of Heads")
plt.xlabel("No. of Heads")
plt.ylabel("Probability")
plt.legend()
plt.show()

for Nc2 in Nc_list:   
    df2,probability2,possible_outcomes2,com_prob_tails2,com_prob_heads2,Not2,bin_dis2,fluctuation2=(Prob(Nt2, Nc2))
    #print("\nFor ",Nt2,"No. of trials and ",Nc2,"No. of Coins")
    #print(df2)
    plt.plot(possible_outcomes2,probability2,label=f"{Nc2} no. of Coins")
    plt.scatter(possible_outcomes2,probability2)
plt.minorticks_on()
plt.grid(b=True,which="both",axis="both")
plt.title("Probability Vs No. of Heads")
plt.xlabel("No. of Heads")
plt.ylabel("Probability")
plt.legend()
plt.show()

df3,probability3,possible_outcomes3,com_prob_tails3,com_prob_heads3,Not3,bin_dis3,fluctuation3=Prob(Nt3, Nc3)
print("\nFor ",Nt3,"No. of trials and ",Nc3,"No. of Coins")
print(df3)
plt.plot(Not3,com_prob_heads3)
plt.plot(Not3,com_prob_tails3)
plt.minorticks_on()
plt.grid(b=True,which="both",axis="both")
plt.title("Probability Vs No. of Trials")
plt.xlabel("No. of Trials")
plt.ylabel("Probability")
plt.legend(["Probability of Head","Probability of Tail"])
plt.show()

plt.plot(Not3[:500],com_prob_heads3[:500])
plt.plot(Not3[:500],com_prob_tails3[:500])
plt.minorticks_on()
plt.grid(b=True,which="both",axis="both")
plt.title("Probability Vs No. of Trials")
plt.xlabel("No. of Trials")
plt.ylabel("Probability")
plt.legend(["Probability of Head","Probability of Tail"])
plt.show()

Not3=np.array(Not3)
fluctuation3=np.array(fluctuation3)
diff=[]
klist=np.linspace(0.5,1,5)
for k in klist:
    y=1/Not3**k
    d=y-fluctuation3
    diff.append(d)
    
for i in range(len(diff)):    
    plt.plot(Not3,diff[i],label=f"for {klist[i]}")
plt.plot(Not3,fluctuation3,color="black")
plt.minorticks_on()
plt.grid(b=True,which="both",axis="both")
plt.title("Fluctuation Graph")
plt.xlabel("No. of Trials")
plt.ylabel("Fluctuation")
plt.legend()
plt.show()
