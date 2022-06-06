#!/usr/bin/python
"""Penetrance Estimator plotting application core; calculates ftilde and
ftildestar penetrance estimates.

"""

# Authors: Sang-Cheol Seok <sangcheol.seok@nationwidechildrens.org>
#                          <scs@mathmed.org>
#          Jo Valentine-Cooper <jvc@mathmed.org>

"""
% Assumes parents are phenotypically unknown and 11 x 12
%  allele '2' is the disease allele
%  The disease is rare dominant disease
%  There is at least one kid with '12' and affected.
%    That is, the first kid is '12' and affected.
%  Apply the k-model in the Hodge(1998) paper page1218.
%  Prob(family ascertained | r affected) = const * (r^k + t)
%  HetModel, prob(aff|12) = alpha + f - alpha * f 
%     where r = # of kids with 'aff' and '12'
%     Use const = 1/1000 and k= -1,0, 1, 2

%Inputs
%   s: # of kids
%   alpha: prob(aff | '11')
%       f: prob(aff | '12')
%Outputs
%  ftilde
%  ftildestar
"""

from math import factorial
from itertools import permutations, repeat, product
import random

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd   # construct dataframe
import seaborn as sns

from sinaplot import sinaplot


# Useful constants so I don't have to keep typing the TeX over and over
FT = r"$\tilde{f}$"
FTS = r"$\tilde{f}^*$"
ALPHA = r"$\alpha$"

# Shorthand for our data printing options context.
DATA_PRINT_OPTIONS = pd.option_context(
        'display.max_rows', None,
        'display.max_columns', None,
        'display.precision', 4)

def nCr(n,r):
    return factorial(n) / factorial(r) / factorial(n-r)

def simNucPedigree(N, s, alpha, f, kvalue):
    # simulate Pedigree
    # FIXME: Look into ways of optimizing this because it seems to be repeating
    # work a lot.
    
    if s < 10 :
        const = 1/100
    else :
        const = 1/1000
        
    hetProbaff= alpha + f - alpha * f

    Aff12 = 0.0
    Aff12NoProb = 0.0
    UnAff12 = 0.0
    
    famId =0
    totalfam=0
    while famId < N :
        totalfam = totalfam +1
        GTs = []
        for i in range(0,s):
            n = random.randint(1,2)
            GTs.append(n)
        PTs = [0] * s
        for k in range(s):
            if GTs[k]==1: # % '11'
                if random.random() < alpha :
                    PTs[k]=2
                else:
                    PTs[k]=1
                    
            else: # % '12'
                if random.random() < hetProbaff :
                    PTs[k]=2
                else:
                    PTs[k]=1

        temp =10*np.array(PTs)+np.array(GTs)
        r = sum(p == 22 for p in temp) # sum((PTs==2)&(GTs==2))
        t = sum(p == 2 for p in PTs)
        if r> 0:
            probAsct =   const * (pow(r, float(kvalue)) + t)   #% prob(ascertaint) using k-method const * (r^k +t)
        else:
            probAsct =   0
        
        if random.random() < probAsct:
            Aff12 = Aff12 + r
            temp =10*np.array(PTs)+np.array(GTs)
            UnAff12 = UnAff12 + sum(p == 12 for p in temp)
            if r > 0 :
                Aff12NoProb = Aff12NoProb + r - 1
            famId = famId + 1
    if (Aff12NoProb +UnAff12) <=0:
        ftilde = 2
    else:
        ftilde =  Aff12NoProb/(Aff12NoProb +UnAff12)
    if (Aff12 +UnAff12) <=0:
        ftildestar = 2
    else:
        ftildestar =  Aff12/(Aff12 +UnAff12)
        
    #print (N, s, GTs, PTs, Aff12NoProb, UnAff12, Aff12,r, t, ftilde , ftildestar)
    #%sprintf('Aff12=%f Aff12NoProb=%f UnAff12=%f   ftilde = %f  ftildestar = %f',Aff12 ,Aff12NoProb, UnAff12, ftilde,ftildestar )

    return ftilde, ftildestar

def exactNucPedigree(s, alpha, f, kvalue):
    # FIXME: Look into ways of optimizing this as well, for the same reasons as
    # above.
    
    hetProbaff= alpha + f - alpha * f
    const = 0.001
    
    #loop over all possible cases of phenotypes for the s-1 kids
    pools = sorted(list(set(
            permutations([1 for x in range(s)] + [2 for x in range(s)], s))))
    numRows=len(pools)
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
        
            probPT_GT = 1  #% prob of PT given GT
            for k in range(s):
                if GTs[k]==1: # % '11'
                    if PTs[k]==1 :
                        probPT_GT = probPT_GT * (1-alpha)
                    else:
                        probPT_GT = probPT_GT * alpha
                    
                else: # % '12'
                    if PTs[k]==1 :
                        probPT_GT = probPT_GT * (1 - hetProbaff) #%(1-f)
                    else:
                        probPT_GT = probPT_GT * hetProbaff #%* f
            temp =10*np.array(PTs)+np.array(GTs)
            r = sum(p == 22 for p in temp) # sum((PTs==2)&(GTs==2))
            t = sum(p == 2 for p in PTs)
            
            if r> 0:
                probAsct =   const * (pow(r, float(kvalue)) + t)   #% prob(ascertaint) using k-method const * (r^k +t)
            else:
                probAsct =   0 #%const * t   %  0
                
            totalProb[j] =   probAsct* probGT*probPT_GT
            Aff12 = Aff12 + r*totalProb[j]
            if r>0:
                Aff12NoProb = Aff12NoProb+ (r-1)*totalProb[j]
            
            UnAff12 = UnAff12 + sum(p == 12 for p in temp)*probAsct* probGT*probPT_GT
        
        totalProbGT[numkids]= sum(totalProb)
    ftilde =  Aff12NoProb/(Aff12NoProb +UnAff12)
    ftildestar =  Aff12/(Aff12 +UnAff12)
    #%sprintf('Aff12=%f Aff12NoProb=%f UnAff12=%f   ftilde = %f  ftildestar = %f',Aff12 ,Aff12NoProb, UnAff12, ftilde,ftildestar )

    return ftilde, ftildestar

def listify_noniter(var):
    """Given a variable, returns it if it's an iterable, otherwise returns that
    variable as a one-item tuple.
    
    """
    try:
        iter(var)
    except TypeError:
        return (var, )
    else:
        return var

class FakeProgress(object):
    """Fakes progress bar interface."""
    def tick(self):
        pass

def calc_stats(alphavals, fvals, svals, kvals, Nvals=None, NumReps=None,
        to_stdout=False, progress=FakeProgress()):
    """Given one or more values each for two of alpha, f, s, and k (and one
    value each for the others) OR given one value each for all four plus number
    of families and number of simulations of families to run, calculates ftilde
    and ftildestar for all values of those variables and returns them in a
    dataframe.
    
    """
    # We have two possible operating modes:
    # 1) Penetrance estimate distributions as a function of number of families.
    # 2) Calculation of penetrance estimation as a function of one of the
    #   component variables (typically alpha, but hypothetically we could do
    #   others) for multiple values of another component variable (typically k)
    # I (Jo) don't actually know the distinctions all that well and was
    # referred to the (at the time unpublished) paper when I asked, so please
    # bear with me. :)
    # From a data structure perspective, though, this can be summed up as
    # "we have two iterables and some constants and they're used to generate a
    # pair of lists of values".
    # For mode 1, the iterables are Nvals and range(NumReps).
    # For mode 2, the iterables are kvals and alphavals.
    # So we make the following assumptions:
    # * If Nvals and/or NumReps are None, we're using operating mode 2 above.
    # * Otherwise, we're using mode 1.
    alphavals = listify_noniter(alphavals)
    for alphaval in alphavals:
        if alphaval < 0 or alphaval > 1:
            raise ValueError("alpha must be between 0 and 1")
    fvals = listify_noniter(fvals)
    for fval in fvals:
        if fval < 0 or fval > 1:
            raise ValueError("f must be between 0 and 1")
    svals = listify_noniter(svals)
    for sval in svals:
        if sval < 2:
            raise ValueError("s must be greater than or equal to 2")
    kvals = listify_noniter(kvals)
    for kval in kvals:
        if kval < -10 or kval > 10:
            raise ValueError("k must be between -10 and 10")
    Nvals = listify_noniter(Nvals)
    for Nval in Nvals:
        if Nval is not None and Nval < 1:
            raise ValueError("N must be a positive integer")
    if NumReps is not None:
        if NumReps < 10:
            raise ValueError("NumReps must be greater than or equal to 10")
    
    itercount = [1 for vals in (alphavals, fvals, svals, kvals, Nvals)
            if len(vals) > 1]
    if NumReps is not None:
        itercount.append(1) # NumReps is always an iterator
    if len(itercount) > 2:
        raise ValueError("only 2 vars may have multiple values")
    
    # k gets reversed so the legend appears in the right order
    allvars = product(alphavals, fvals, svals, reversed(kvals), Nvals,
            (None,) if NumReps is None else range(NumReps))
    progress.total = np.prod([len(val) for val in (alphavals, fvals, svals,
            kvals, Nvals)] + [1 if NumReps is None else NumReps, ])
    
    stats = {}
    for alpha, f, s, k, N, RepNum in allvars:
        #print(f"alpha={alpha}, f={f}, s={s}, k={k}, N={N}, "
        #        f"RepNum={RepNum}")
        stats.setdefault(ALPHA, []).append(alpha)
        stats.setdefault('f', []).append(f)
        stats.setdefault('s', []).append(s)
        stats.setdefault('k', []).append(k)
        ft, fts = (2, 2)
        if RepNum is not None and N is not None:
            progress.message = (f"Calculating: α={alpha}, f={f}, "
                    f"s={s}, k={k}, N={N}, RepNum={RepNum}")
            stats.setdefault("N", []).append(N)
            stats.setdefault("RepNum", []).append(RepNum)
            while (ft >= 2 and fts >= 2):
                ft, fts = simNucPedigree(N, s, alpha, f, k)
        else:
            progress.message = (f"Calculating: alpha={alpha}, f={f}, "
                    f"s={s}, k={k}")
            ft, fts = exactNucPedigree(s, alpha, f, k)
        stats.setdefault(FT, []).append(ft)
        stats.setdefault(FTS, []).append(fts)
        progress.tick()
    
    results = pd.DataFrame(stats)
    if to_stdout:
        with DATA_PRINT_OPTIONS:
            print(results)
    
    return results

def output_summarystats(df, filename):
    """Given one of our figure dataframes and a filename, outputs our summary
    stats info to that file.
    
    """
    
    # The "summary stats" for "as a function of alpha" (fig2) pretty much are a
    # specialized form of the raw output. For the estimate distributions based
    # on number of families, though, we have actual summarization. So we need
    # to figure out which is which.
    
    try:
        df['N']
    except KeyError:
        # "as a function of alpha" - no summarization, just formatting
        ftdf, ftsdf = df.pivot('k', ALPHA, FT), df.pivot('k', ALPHA, FTS)
        with open(filename, "w") as outfile:
            with DATA_PRINT_OPTIONS:
                outfile.write("\n\n".join((
                        # lack of commas for first 3 lines is deliberate
                        "==================================\n"
                        "Figure2. Statistics\n"
                        "==================================\n"
                        f"f={df['f'][0]}, s={df['s'][0]}",
                        f"f_tilde\n{df.pivot('k', ALPHA, FT)}",
                        f"f_tilde_star\n{df.pivot('k', ALPHA, FTS)}"
                        )))
    else:
        # distribution plot - we need to actually summarize :)
        base_ft = df.pivot('RepNum', 'N', FT)
        base_fts = df.pivot('RepNum', 'N', FTS)
        ft_stats = pd.DataFrame({
                "mean": base_ft.mean(0),
                "std": base_ft.std(0),
                "med": base_ft.median(0)
                }).transpose()
        fts_stats = pd.DataFrame({
                "mean": base_fts.mean(0),
                "std": base_fts.std(0),
                "med": base_fts.median(0)
                }).transpose()
        with open(filename, "w") as outfile:
            with DATA_PRINT_OPTIONS:
                outfile.write("\n\n".join((
                        # lack of commas for first 5 lines is deliberate
                        "==================================\n"
                        "Figure1. Statistics\n"
                        "==================================\n"
                        f"alpha={df[ALPHA][0]}, f={df['f'][0]}, "
                        f"s={df['s'][0]}, k={df['k'][0]}, "
                        f"NumReps={df['RepNum'].iat[-1] + 1}",
                        f"f_tilde\n{ft_stats}",
                        f"f_tilde_star\n{fts_stats}"
                        )))

def generate_and_display_test_figure(fig1df, fig2df):
    """Given two dataframes - one with our N values and our two sets of Y-axis
    values (ftilde, ftildestar), and the other with our alpha and k
    values and our two sets of Y-axis values, plots the first dataset's values
    by N in sinaplots and the second dataset's values by alpha in lineplots,
    and display the results.
    
    """
    
    plt.figure(1)
    
    plt.subplot(2, 2, 1)
    sinaplot(x="N", y=FT, data=fig1df, edgecolor="black", alpha=.5,
            violin=False)
    plt.ylabel(FT, rotation=0)
    plt.ylim([0,1])
    plt.axhline(0.5, alpha=0.5, dashes=(5,2))
    
    plt.subplot(2, 2, 2)
    sinaplot(x="N", y=FTS, data=fig1df, edgecolor="black", alpha=.5,
            violin=False)
    plt.ylabel(FTS, rotation=0)
    plt.ylim([0,1])
    plt.axhline(0.5, alpha=0.5, dashes=(5,2))
    
    plt.subplot(2, 2, 3)
    sns.lineplot(data=fig2df.pivot(ALPHA, "k", FT), markers=True)
    plt.ylabel(FT, rotation=0)
    plt.ylim([0,1])
    plt.axhline(0.5, alpha=0.5, dashes=(5,2))
    
    plt.subplot(2, 2, 4)
    sns.lineplot(data=fig2df.pivot(ALPHA, "k", FTS), markers=True)
    plt.ylabel(FTS, rotation=0)
    plt.ylim([0,1])
    plt.axhline(0.5, alpha=0.5, dashes=(5,2))
    
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    # Default "test" input values
    # Figure 1 which goes to the 1st row
    alpha = 0.0
    f = 0.5
    s = 2
    k = 1
    N = [1, 2, 4, 6, 8, 10, 30, 50]
    NumReps = 100
    fig1df = calc_stats(alpha, f, s, k, N, NumReps)
    
    # Figure 2 which goes to the 2nd row.
    alpha = [0.001, 0.01, 0.1, 0.2, 0.3, 0.4, 0.5]
            # now up to 0.5 and these values are fixed.  5/30/2022
    f = 0.5
    s = 3
    k = [-1, 0, 1, 2]
    fig2df = calc_stats(alpha, f, s, k)
    
    generate_and_display_test_figure(fig1df, fig2df)
