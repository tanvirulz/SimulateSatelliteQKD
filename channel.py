#This block is written by Tom Vergoossen
import numpy as np
import math
import mpmath
from scipy import special

def calculateDistance(T = np.linspace(-200,200,201), distMin = 500e3, maxElevation = 60):
    re = 6378
    radius = (re+distMin/1000)
    v= np.sqrt(398600/(6378+distMin/1000))*1000
    omega = v/radius*1e-3
  
    xOGS = 0
    yOGS = distMin/1000/np.tan(maxElevation/180*np.pi)
    zOGS = re

    xSat = radius*np.sin(omega*T)
    ySat = 0
    zSat = radius*np.cos(omega*T)
   
    dist = 1000*np.sqrt((xSat - xOGS)**2+(ySat - yOGS)**2+(zSat - zOGS)**2)
    xdist = 1000*(xSat-xOGS)
    ydist = 1000*(ySat-yOGS)
    zdist = 1000*(zSat-zOGS)
    
    elevationAngle = np.arctan((zSat - zOGS)/np.sqrt((xSat-xOGS)**2+((ySat-yOGS)**2)))
    azimuth = np.arctan((xSat-xOGS)/(ySat-yOGS))
   
    return dist, elevationAngle, v, T, xdist, ydist, zdist, azimuth
    
def BeamWidthReceiver(Distance, BeamWaist, Wavelength, MSquare, Jitter, TurbAngle):
    #Beam divergence including atmospheric and pointing effects
    w_0 = BeamWaist
    l = Wavelength
    M2 = MSquare
    d = Distance

    #Diffraction loss (Gaussian beam theory)
    z_R = np.pi*(w_0**2)/l;                                  #Rayleigh range
    w_d = M2*w_0*np.sqrt(1+np.square(d/z_R));                #Beam width at receiving telescope
    w_d = np.sqrt(w_d**2+(d*np.tan(Jitter))**2+(d*np.tan(TurbAngle))**2)
    return w_d

def ChannelLoss(Distance, WidthAtReceiver, ApertReceive, Elevation, AtmosAtten, PointingOffset = 0):
    
    if Distance.size > 1:
        Off = np.zeros((len(Distance)))
        w_d = np.zeros((len(Distance)))
        Ltr = np.zeros((len(Distance)))
        Losses = np.zeros((len(Distance)))
    
        for i in range(len(Distance)):
            #see https://arxiv.org/ftp/arxiv/papers/1605/1605.04241.pdf
            r = ApertReceive/2
            w_d[i] = WidthAtReceiver[i]
            d = Distance[i]
            e = Elevation[i]    
            Off[i] = np.tan(PointingOffset)*d
        
            Lgeom = np.exp(-2*(Off[i]**2)/(w_d[i]**2))*mpmath.nsum(lambda k:((2**k)*Off[i] **(2*k))/(w_d[i]**(2*k)*math.gamma(k+1)**2)*special.gammainc((float(k)+1), 2*(r**2)/(w_d[i]**2)),[0, np.inf])
            Lgeom = 10*np.log10(float(Lgeom)) 
            Ltr[i] = AtmosAtten*(1/np.cos((90-e)*np.pi/180))
        #    Ltr = -3
            Losses[i] = Lgeom + Ltr[i] 
    else:
        r = ApertReceive/2
        w_d = WidthAtReceiver
        d = Distance
        e = Elevation    
        Off = np.tan(PointingOffset)*d
    
        Lgeom = np.exp(-2*(Off**2)/(w_d**2))*mpmath.nsum(lambda k:((2**k)*Off **(2*k))/(w_d**(2*k)*math.gamma(k+1)**2)*special.gammainc((float(k)+1), 2*(r**2)/(w_d**2)),[0, np.inf])
        Lgeom = 10*np.log10(float(Lgeom)) 
        Ltr = AtmosAtten*(1/np.cos((90-e)*np.pi/180))
        #    Ltr = -3
        Losses = Lgeom + Ltr
    return Losses,w_d,Off,Ltr

######################### Qubesat parameters ##################################
ApertSat = 0.08
MSquare = 1.8
PointingJitter = 5e-6
CoeffWaist = 0.89
PointingOffset = 0
Waist = ApertSat/2*CoeffWaist

#Source parameters
Brightn = 56e6
IntrinsQBER = 0.01
CoincTimeSo = 1e-9
Wavel = 785e-9
DetEffSo = 0.25
DarkCSo = 15000
DeadTSo = 200e-9
DetJitter = 640e-12
Pap = 0.005

#Ground station parameters
ApertGr = 0.6
Backgr = 1.508e12
#Backgr = 0
AtmosAtten = -3
OpticEffGr = 0.5
FocalLenGr = 3.962
TurbAngle = 4e-6

#Quantum receiver parameters
DetEffRe = 0.25
DeadTRe = 27e-9
DarkCRe = 500
SensSizeRe = 500e-6
CoincTimeRe = 1e-9
SensSizeRe = 500e-6

#Protocol parameters
Type = "BBM92"
BasisRecon = 0.5
ErrCoEff = 1.1
CoincEff = 1
SizeTStamp = 64
CompRate = 1/4
beta = 1e-10/10
ecorr = 1e-10

#################################Single Pass###################################
#Scenario
MaxElevation = 88
Altitude = 500e3

T = np.linspace(-300,300,1201)
a = calculateDistance(T,Altitude,MaxElevation)

Distance = a[0]
Elevation = a[1]*180/np.pi
Time = a[3]

b = BeamWidthReceiver(Distance,Waist,Wavel,MSquare,PointingJitter,TurbAngle)
WidthR = b
c = ChannelLoss(Distance,b,ApertGr,Elevation,AtmosAtten,PointingOffset)
#print (len(c))
Losses = c[0]
Ltr = c[3]
#print (len(Distance) )

TotalLossDown = Losses+10*np.log10(OpticEffGr)

#print (len(b))


def ChannelLossScalar(Distance, WidthAtReceiver, ApertReceive, Elevation, AtmosAtten, PointingOffset = 0):
    r = ApertReceive/2
    w_d = WidthAtReceiver
    d = Distance
    e = Elevation    
    Off = np.tan(PointingOffset)*d

    Lgeom = np.exp(-2*(Off**2)/(w_d**2))*mpmath.nsum(lambda k:((2**k)*Off **(2*k))/(w_d**(2*k)*math.gamma(k+1)**2)*special.gammainc((float(k)+1), 2*(r**2)/(w_d**2)),[0, np.inf])
    Lgeom = 10*np.log10(float(Lgeom)) 
    Ltr = AtmosAtten*(1/np.cos((90-e)*np.pi/180))
    #    Ltr = -3
    Losses = Lgeom + Ltr
    return Losses,w_d,Off,Ltr
    
    distance_list = np.array(r_list)
    elev_list = np.array(alt_list)
    
