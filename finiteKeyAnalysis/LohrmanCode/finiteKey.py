# -*- coding: utf-8 -*-
"""
Created on Fri Nov 13 09:09:46 2020

@author: alloh
"""




import numpy
import matplotlib.pyplot as plt
from scipy.optimize import minimize

# from "Security analysis of quantum key distribution with small block length and its application to quantum space communications

def shannon(Q):
    return -Q*numpy.log2(Q)-(1-Q)*numpy.log2(1-Q); 


def optKeyFrac(p, m, QBER, s = 6, f = 1.2):
    beta,xi, nu = p
#   print("ADS")
#    print(beta, xi, nu)
    return 1-keyFraction(m, QBER, beta, xi, nu, s, f)

def keyFraction(m,QBER,beta,xi,nu, s = 6,f = 1.2): # m is the initial raw raw key length
    
    
    delta = QBER #QBER
    t = -numpy.log2(10**-(s+2)) #correctness
    r = 1.18*shannon(delta) #leaked information
    epsQKD = 10**(-s)#2**-16
    #m is the raw raw sifted key
    k = numpy.floor(beta*m)
    n = m-k # that's what's left after paramaeter estimation

    Gamma= 1/(m*(delta+xi)+ 1) +1/(m-m*(delta+xi)+ 1)
    PPElim = numpy.exp(-2*m*k*xi**2/(n+1)) + numpy.exp(-2*Gamma*((n*(nu-xi))**2-1))
    epe = numpy.sqrt(PPElim)
       
    tot = 4*(epsQKD - 2**-t - 2*epe)**2   
    
    
    l = numpy.log2(tot) - t +n*(1-shannon(delta+nu) - r)
    #print(l, xi, nu, beta)
    
    epa = 0.5*numpy.sqrt(2**(-n*(1-shannon(delta) - r )  + t + l))
#    print(epsQKD > 2**-t + 2*epe + epa)
#    print(xi < nu)
#    print(nu<0.5-delta)
    
    if epsQKD < 2**-t + 2*epe + epa:
        l = 0
        
    if xi > nu:
        l = 0
        
    if nu>0.5-delta:
        l = 0
    return l
    




#L = keyFraction(3100, 0.0455,  0.5,  0.069, 0.1141 )


S = [6,10]
legend = ['s = ' + str(S[0]),'s = ' + str(S[1])]

for s in S:
    N = numpy.linspace(2000,10000, 101)
    QBER = 0.0455
    f = 1.2

    
    
    L = numpy.zeros(len(N))
    LAss =  numpy.zeros(len(N))
    
    beta = numpy.zeros(len(N))
    xi = numpy.zeros(len(N))
    nu = numpy.zeros(len(N))
    
    
    for i in range(len(N)):
        x0 = [ 0.5,  0.069, 0.1141]
        bnds = ((1e-8,0.5), (1e-8,0.5-QBER),(0,0.5-QBER))
        par = minimize(optKeyFrac, x0 ,bounds=bnds, args=(N[i], QBER, s , f ), method = 'Nelder-Mead')
        L[i] = keyFraction(N[i], QBER,  par.x[0], par.x[1], par.x[2],  s= s, f=f)
    
    #L = keyFraction(N, QBER,  0.5, 0.0693, 0.114, f = f)
        LAss = N[i]*(1-shannon(QBER)-f*shannon(QBER))
        
        beta[i], xi[i], nu[i] =  par.x[0], par.x[1], par.x[2]
    
    #print(L, L/LAss)
    
    plt.loglog(N, L/N)
    plt.ylim([1e-3, 0.1])
    
plt.legend(legend)
plt.xlabel('Raw key')
plt.ylabel('$l$')


















