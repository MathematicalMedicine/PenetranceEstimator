#!/usr/bin/python
# Developed with python 2.7
# This program calculates Beta_Hat and Beta_HatStar 
# Written by Sangcheol Seok 01/03/2022
#  python BetaApp.py 

# Being cleaned up and updated to python 3 by Jo Valentine-Cooper 20220106

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
from math import factorial
from pylab import figure
from itertools import permutations, repeat
import matplotlib.pyplot as plt



def nCr(n,r):
    return factorial(n) / factorial(r) / factorial(n-r)
    
def exactNucPedigreeAllModels(ModelId, s, alpha, beta, kvalue):
    # FIXME: Look into ways of optimizing this because it seems to be repeating
    # work a lot.
    
    hetProbaff= alpha + beta - alpha * beta
    const = 0.001
    
    #loop over all possible cases of phenotypes for the s-1 kids
    pools = sorted(list(set(
            permutations([1 for x in range(s)] + [2 for x in range(s)], s))))
    numRows=len(pools);
    #print pools
    GTs=[1] * s
    totalProbGT = [0] * s  #zeros(s,1)
    Aff12 = 0.0
    Aff12NoProb = 0.0
    UnAff12 = 0.0
    
    for numkids in range(s):  # we don't need the case for no kid has '12', so we start from at least one '12' case
        
        GTs[numkids]=2
        #otherwise, there are 2 phenotype cases
        #GTs
        probGT=float(nCr(s,numkids+1))/pow(2,s)
                # 2**s # probability of GT is 1/2^s in kids
                # nchoosek is the number of liketerms
        #print ("s is ",s," numkids is ",numkids," probGT is ",probGT, "nCr(s,numkids) is ", nCr(s,numkids), "pow(2,s) is ",pow(2,s))
        totalProb = [0] *  numRows
        #%Step2. Assigning Phenotypes
        for j in range(numRows):
            PTs = pools[j]
        
            probPT_GT = 1 ; #% prob of PT given GT
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
                    probAsct =   const * pow(r, float(kvalue))#% prob(ascertaint) using k-method const * r^k
                else:
                    probAsct = 0
            
            else: #% ModelId ==2 or 3 %probAsct =   const * (r ^ks + t )
                if r> 0:
                    probAsct =   const * (pow(r, float(kvalue)) + t)  ; #% prob(ascertaint) using k-method const * (r^k +t)
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
        
        totalProbGT[numkids]= sum(totalProb)
        #print("numkids=", numkids,"Aff12= ", Aff12, " Aff12NoProb= ", Aff12NoProb," UnAff12=", UnAff12)
    #print("Aff12= ", Aff12, " Aff12NoProb= ", Aff12NoProb," UnAff12=", UnAff12)
    beta_hat =  Aff12NoProb/(Aff12NoProb +UnAff12);
    beta_hatStar =  Aff12/(Aff12 +UnAff12);
    #%sprintf('Aff12=%f Aff12NoProb=%f UnAff12=%f   beta_hat = %f  beta_hatStar = %f',Aff12 ,Aff12NoProb, UnAff12, beta_hat,beta_hatStar )

    return beta_hat, beta_hatStar


def generate_and_display_test_figure(alphas, beta_hat, beta_hatStar):
    """Given our X-axis values (alphas) and lists of our two sets of Y-axis
    values (beta_hat, beta_hatStar), plot these values and display the results.
    
    """
    # FIXME: right now this is just a "pop up a window" operation, and we want
    # a somewhat customized interactive app. we can and should do better going
    # forward.
    
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
    

def isiter(var):
    """Given a variable, tests to see if it's an iterable."""
    # This is largely so I can do a logical xor test in calc_stats rather than
    # chaining together a bunch of try/except/elses together.
    try:
        iter(var)
    except TypeError:
        return False
    else:
        return True

def repeat_noniter(var):
    """Given a variable, returns it if it's an iterable, otherwise returns a
    repeat() iterable of that variable's value.
    
    """
    return var if isiter(var) else repeat(var)

def format_stats(stats):
    """Formats our stats as a nice string for printing."""
    # just an optional quick display of the values, using the formatting scheme
    # from Sang-cheol's original starting code (list of lists, 4 significant
    # digits for each float)
    
    finalformat = []
    finalformat.append(f"alpha={stats['alpha']}, beta={stats['beta']}, "
            f"s={stats['s']}, modelid={stats['modelid']}")
    for stattype, statcat in stats['stats'].items():
        finalformat.append(f"{stattype}:")
        for catname, statsubcat in statcat.items():
            finalformat.append((f"  {catname}: " +
                    str(["{0:0.4f}".format(val) for val in statsubcat])))
    
    return "\n".join(finalformat)

def calc_stats(modelid, kvals, alphavals, betavals, svals, to_stdout=False):
    """Given a single value for modelid (as documented above), one or more
    values for k, and a list of values for ONE of s, alpha, and beta (and
    single values for the other two), calculates betahat and betahatstar for
    all values of (whichever of s/alpha/beta has multiple values) and returns
    them in a dict.
    
    """
    
    # Test and make sure we only have one iterable
    if not bool(isiter(alphavals) ^ (isiter(betavals) ^ isiter(svals))):
        raise ValueError("one and only one of alphavals/betavals/svals must "
                "be an iterable (list of values)")
    
    allvars = (repeat_noniter(alphavals), repeat_noniter(betavals),
            repeat_noniter(svals))
    
    results = {
            'modelid': modelid,
            'alpha': alphavals,
            'beta': betavals,
            's': svals,
            'k': kvals,
            'stats': {
                'betahat': {},
                'betahatstar': {}
                }
            }
    for k in kvals:
        for alpha, beta, s in zip(*allvars):
                # we have to zip each time thru 'cause each loop through kvals
                # exhausts the generator
            #print(f"k={k}, alpha={alpha}, beta={beta}, s={s}")
            bh, bhs = exactNucPedigreeAllModels(modelid, s, alpha, beta, k)
            results['stats']['betahat'].setdefault(f"k={k}", []).append(bh)
            results['stats']['betahatstar'].setdefault(f"k={k}", []).append(bhs)
    
    if to_stdout:
        print(format_stats(results))
    
    return results


if __name__ == "__main__":
    # Default "test" input values
    modelid = 3     # should be one of 1, 2, 3
    s = 3           # greater than 1
    #alpha = 0.5    # ranges from 0 to 1
    alphas = [0.001, 0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.75]
    beta= 0.5       # ranges from 0 to 1
    #kvalue = 1     # should be one of -1, 0, 1, 2
    kvalues = [-1, 0, 1, 2]
    
    stats = calc_stats(modelid, kvalues, alphas, beta, s, True)['stats']
    
    generate_and_display_test_figure(alphas, list(stats['betahat'].values()),
            list(stats['betahatstar'].values()))

