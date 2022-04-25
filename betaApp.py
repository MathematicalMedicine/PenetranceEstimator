#!/usr/bin/python
# Developed with python 2.7
# This program calculates Beta_Hat and Beta_HatStar 
# Written by Sangcheol Seok 01/03/2022
#  python BetaApp.py 

"""
% Assumes parents are phenotifically unknown and 11 x 12
%  allele '2' is the disease allele
%  The disease is rare dominant disease
%  There is at least one kid with '12' and affected.
%    That is, the first kid is '12' and affected.
%  Apply the k-model in the Hodge(1998) paper page1218.
%  Prob(family ascertained | r affected)
% =  const * (r^k ) if ModelId=1
% =  const * (r^k + t) if ModelId=2
% =  const * (r^k + t) if ModelId=3 HetModel, prob(aff|12) = alpha + beta - alpha * beta ;
%     where r = # of kids with 'aff' and '12'
%     Use const = 1/1000 and k= -1,0, 1, 2

%Inputs
%   ModelId : 1  probAsct =   const * r ^k ;
%             2  probAsct =   const * (r ^k + t ) prob(aff|12)=beta
%             3  probAsct =   const * (r ^k + t ) prob(aff|12)=alpha + beta - alpha * beta
%   s: # of kids
%   alpha: prob(aff | '11')
%    beta: prob(aff | '12')
%Outputs
%  beta_hat
%  beta_hatStar
"""
import numpy as np
import sys   # for command line arguments
import math
import string
import random # random number generator
import operator # to use key to sort list
from pylab import *
from numpy import *
from pylab import load
import matplotlib.pyplot as plt
#from math import *



def nCr(n,r):
    f = math.factorial
    return f(n) / f(r) / f(n-r)
    
""" Create all possible phenotypes for the s-1 kids
"""
def binarySequence(curSeq,depth, Pool):
    
    #print ("depth is ",depth, "  curSeq is ", curSeq)
    if depth > 0:
        newSeq=list(curSeq)
        newSeq.append(1)
        Pool = binarySequence(newSeq,depth-1,Pool)
        newSeq=list(curSeq)
        newSeq.append(2)
        Pool = binarySequence(newSeq,depth-1,Pool)
    else:
        Pool.append(curSeq)

    return Pool

def exactNucPedigreeAllModels(ModelId, s, alpha, beta, kvalue):

    hetProbaff= alpha + beta - alpha * beta
    const = 0.001
    
    #loop over all possible cases of phenotypes for the s-1 kids
    Pools= binarySequence([],s,[])
    numRows=len(Pools);
    #print Pools
    GTs=[1] * s
    totalProbGT = [0] * s  #zeros(s,1)
    Aff12 = 0.0
    Aff12NoProb = 0.0
    UnAff12 = 0.0

    for i in range(s):  # we don't need the case for no kid has '12', so we start from at least one '12' case

        GTs[i]=2
        #otherwise, there are 2 phenotype cases
        #GTs
        probGT=float(nCr(s,i+1))/pow(2,s) #2**s # probability of GT is 1/2^s in kids
                                  # nchoosek is the number of liketerms
        #print ("s is ",s," i is ",i," probGT is ",probGT, "nCr(s,i) is ", nCr(s,i), "pow(2,s) is ",pow(2,s))
        totalProb = [0] *  numRows
        #%Step2. Assigning Phenotypes
        for j in range(numRows):
            PTs = Pools[j]
        
            probPT_GT =  1 ; #% prob of PT given GT
            for k in range(s):
                if GTs[k]==1: # % '11'
                    if PTs[k]==1 :
                        probPT_GT = probPT_GT * (1-alpha)
                    else:
                        probPT_GT = probPT_GT * alpha
                    
                else: # % '12'
                
                    if ModelId ==3 :
                        if PTs[k]==1 :
                            probPT_GT = probPT_GT * (1 - hetProbaff); #%(1-beta);
                        else:
                            probPT_GT = probPT_GT * hetProbaff; #%* beta;
                        
                    else:
                        if PTs[k]==1 :
                            probPT_GT = probPT_GT * (1-beta);
                        else:
                            probPT_GT = probPT_GT * beta;
                      
            temp =10*np.array(PTs)+np.array(GTs)
            r = sum(p == 22 for p in temp) # sum((PTs==2)&(GTs==2))
            t = sum(p == 2 for p in PTs)
         
            if ModelId ==1 :
                #%probAsct =   const * r ^ks(2)
                if r> 0 :
                    probAsct =   const * pow(r ,kvalue )#% prob(ascertaint) using k-method const * r^k
                else:
                    probAsct = 0
            
            else: #% ModelId ==2 or 3 %probAsct =   const * (r ^ks + t )
                if r> 0:
                    probAsct =   const * (pow(r, kvalue) + t)  ; #% prob(ascertaint) using k-method const * (r^k +t)
                else:
                    probAsct =   0; #%const * t  ; %  0;
                
            totalProb[j] =   probAsct* probGT*probPT_GT
            Aff12 = Aff12 + r*totalProb[j]
            if r>0:
                Aff12NoProb = Aff12NoProb+ (r-1)*totalProb[j]
            
            #UnAff12 = UnAff12 + sum((PTs==1)&(GTs==2))*probAsct* probGT*probPT_GT
            UnAff12 = UnAff12 + sum(p == 12 for p in temp)*probAsct* probGT*probPT_GT
            #print(totalProb[j] ,  probAsct, probGT,probPT_GT)
            #print("i=", i,"GTs", GTs, " PTs", PTs,"probAsct= ", probAsct," probGT=", probGT, "probPT_GT= ",probPT_GT,  "totalProb ", totalProb[j])
            #print("j=", j,"Aff12= ", Aff12, " Aff12NoProb= ", Aff12NoProb," UnAff12=", UnAff12, "r =",r, "t=",t)
        
        totalProbGT[i]= sum(totalProb)
        #print("i=", i,"Aff12= ", Aff12, " Aff12NoProb= ", Aff12NoProb," UnAff12=", UnAff12)
    #print("Aff12= ", Aff12, " Aff12NoProb= ", Aff12NoProb," UnAff12=", UnAff12)
    beta_hat =  Aff12NoProb/(Aff12NoProb +UnAff12);
    beta_hatStar =  Aff12/(Aff12 +UnAff12);
    #%sprintf('Aff12=%f Aff12NoProb=%f UnAff12=%f   beta_hat = %f  beta_hatStar = %f',Aff12 ,Aff12NoProb, UnAff12, beta_hat,beta_hatStar )

    
    
    
    return beta_hat, beta_hatStar

#  Input values
ModelId=3 # should be one of 1, 2, 3
s=3        # greater than 1
#alpha=0.5  # ranges from 0 to 1
beta= 0.5  # ranges from 0 to 1
#kvalue= 1  # should be one of -1, 0,1,2
alphas=[0.001, 0.01, 0.1, 0.2, 0.3, 0.4,0.5,0.6 ,0.75]
kvalues=[-1, 0, 1,2]
beta_hat = [] #[[0]*len(alphas)]*len(kvalues)
beta_hatStar = [] #[[0]*len(alphas)]*len(kvalues)

for i in range(len(kvalues)):
    bh = [0] *len(alphas)
    bhS = [0] *len(alphas)
    for j in range(len(alphas)):
        [bh[j],bhS[j]]= exactNucPedigreeAllModels(ModelId, s, alphas[j], beta, kvalues[i])
        #print ("i ", i,"j ", j, " kvalues[i] ", kvalues[i], " alphas[j] ", alphas[j])
    beta_hat.append(bh)
    beta_hatStar.append(bhS)

for i in range(len(kvalues)):
    print ["{0:0.4f}".format(i) for i in beta_hat[i]]
print()
for i in range(len(kvalues)):
        print ["{0:0.4f}".format(i) for i in beta_hatStar[i]]

figure(1)
plt.subplot(2, 1, 1)
plt.plot(alphas,beta_hat[0],'ko-',alphas,beta_hat[1],'rx-',alphas,beta_hat[2],'g*-',alphas,beta_hat[3],'cs-')
plt.xlabel('alpha')
plt.ylabel('Beta_Hat')
plt.axis([0,1,0,1])

plt.subplot(2, 1, 2)
plt.plot(alphas,beta_hatStar[0],'ko-',alphas,beta_hatStar[1],'rx-',alphas,beta_hatStar[2],'g*-',alphas,beta_hatStar[3],'cs-')
plt.xlabel('alpha')
plt.ylabel('Beta_HatStar')
plt.axis([0,1,0,1])

plt.show()


