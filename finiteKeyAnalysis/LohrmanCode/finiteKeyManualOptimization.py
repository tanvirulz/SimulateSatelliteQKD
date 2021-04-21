# -*- coding: utf-8 -*-
"""
Created on Fri Nov 13 09:09:46 2020

@author: alloh
"""

import numpy
import matplotlib.pyplot as plt

# from "Security analysis of quantum key distribution with small block length and its application to quantum space communications

def shannon(Q):
    return -Q*numpy.log2(Q)-(1-Q)*numpy.log2(1-Q); 


def optKeyFrac(p, m, QBER, s = 6, f = 1.2):
    beta,xi, nu = p
    return 1-keyFraction(m, QBER, beta, xi, nu, s, f)

def keyFraction(m,QBER,beta,xi,nu, s = 6,f = 1.2): # m is the initial raw raw key length
    
    
    delta = QBER #QBER
    t = -numpy.log2(10**-(s+2)) #correctness
    r = f*shannon(delta) #leaked information
    epsQKD = 10**(-s)#2**-16
    #m is the raw raw sifted key
    k = numpy.floor(beta*m)
    n = m-k # that's what's left after paramaeter estimation

    Gamma= 1/(m*(delta+xi)+ 1) +1/(m-m*(delta+xi)+ 1)
    PPElim = numpy.exp(-2*m*k*xi**2/(n+1)) + numpy.exp(-2*Gamma*((n*(nu-xi))**2-1))
    epe = numpy.sqrt(PPElim)
       
    tot = 4*(epsQKD - 2**-t - 2*epe)**2   
    
    l = numpy.log2(tot) - t +n*(1-shannon(delta+nu) - r)    
    epa = 0.5*numpy.sqrt(2**(-n*(1-shannon(delta+nu) - r )  + t + l))
   
    if epsQKD < 2**-t + 2*epe + epa:
        l =  -numpy.inf
        
#    if xi > nu:
#        l =  -numpy.inf
#        
#    if nu>0.5-delta:
#        l =  -numpy.inf
#        
#    if beta > 0.5:
#        l =  -numpy.inf
    return l


def minimize(n, qber, s = 10, f = 1.2, loops = 1, steps =101):
      
    par = [0,0,0]


    bounds= [[0,0.5],[0,0.5-qber],[0,0.5-qber]]


    lMax = -numpy.inf
    for loop in range(loops):
#        if bounds[0][0]<0:bounds[0][0]=0;
#        if bounds[0][0]<0:bounds[1][0]=0;
        arr0 = numpy.linspace(bounds[0][0], bounds[0][1] , steps)
        arr1 = numpy.linspace(bounds[1][0], bounds[1][1] , steps)
        
        for i in range(steps): 
            beta = arr0[i]

            for j in range(steps):
                nu = arr1[j]
#
#                if bounds[2][0] < 0:bounds[2][0] = 0
#                
                if bounds[2][1] > nu:
                    arr2 =  numpy.linspace(bounds[2][0], nu, steps)
                else:
                    arr2 =  numpy.linspace(bounds[2][0], bounds[2][1], steps)
                
                for q in range(steps):
                    xi = arr2[q]
                    
                    l = keyFraction(n, qber, beta,xi,nu,  s= s, f=f)
                    if l > lMax:
                        if beta > 0.5 or beta <= 0:
                            print('beta range error')
                            continue
                        if xi > nu or xi<=0:
                            print('xi<>nu comparison error')
                            continue
                        if nu > 0.5-qber or nu<=0:
                            print('nu error')
                            continue
                        lMax = l
                        par = [1*beta,1*xi,1*nu]
                        

        rng = 2
        bounds[0] = [par[0]-rng*(numpy.diff(arr0)[0]), par[0]+rng*(numpy.diff(arr0)[0])]
        bounds[1] = [par[2]-rng*(numpy.diff(arr1)[0]), par[2]+rng*(numpy.diff(arr1)[0])]
        bounds[2] = [par[1]-rng*(numpy.diff(arr2)[0]), par[1]+rng*(numpy.diff(arr2)[0])]

        steps = 11
        
    return lMax, par
                        
            
            
#        
#L, PAR = minimize(5000, 0.0455)  
#print(L, PAR)   
#

S = [6]
#legend = ['s = ' + str(S[0])]#,'s = ' + str(S[1])]

for s in S:
    #N = 10**numpy.linspace(3,4,31)
    N = numpy.linspace(3100,10000,11)
    QBER = 0.0455
    f = 1.7
 
    L = numpy.zeros(len(N))
    LAss =  numpy.zeros(len(N))
    
    beta = numpy.zeros(len(N))
    xi = numpy.zeros(len(N))
    nu = numpy.zeros(len(N))
    
    
    for i in range(len(N)):
        
        l, par = minimize(N[i], QBER, s=s, f =f)
        print(N[i], l, par)
        L[i] = l
        LAss[i] = N[i]*(1-shannon(QBER)-f*shannon(QBER))
        
        beta[i], xi[i], nu[i] =  par[0], par[1], par[2]

    plt.loglog(N, L/N)

plt.ylim([1e-4, 0.1])
plt.xlim([2000, 8000])
plt.legend(legend)
plt.xlabel('Raw key')
plt.ylabel('$l$')
plt.show() 
















