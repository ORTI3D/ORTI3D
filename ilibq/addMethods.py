# -*- coding: utf-8 -*-
"""
Created on Mon Aug 25 15:49:11 2014

@author: olive
"""

from scipy.special import erf,erfc

def AnLinFini(k,v,D,Lx,Lz,x,z):
    indx = range(len(z));indx = indx[-1::-1]
    c=z*0
    sgn,ilast = 1,0
    for i in range(10):
        indx = indx[-1::-1]
        if ilast == 2: 
            sgn = -sgn;ilast=0
        c+=sgn*Seagren(x,z[indx]+Lz*i,Lx,k,v,D)
        ilast+=1
    return c
    
def Seagren(xm,zm,Lx,k,v,Dz): # only for infinite time 
    Da=k*Lx/v;Pe=v*Lx/Dz
    xm1=xm/Lx;zm1=zm/Lx;sqPD=sqrt(Pe*Da)
    a=.5*exp(-zm1*sqPD)*erfc(zm1/2/sqrt(xm1/Pe)-sqrt(Da*xm1))
    b=.5*exp(zm1*sqPD)*erfc(zm1/2/sqrt(xm1/Pe)+sqrt(Da*xm1))
    return a+b
    
def getDissolRate(Lz,k,v,dxv,aT):
    # returns the dissolution rate v is velocity (one value), 
    #dxv is a vector of cell size
    Lx=sum(dxv);nc=len(dxv)
    nx=100;dx=Lx/nx;x=linspace(dx/2,Lx-dx/2,nx)
    nz=20;dz=float(Lz)/nz;z=linspace(dz/2,Lz-dz/2,nz)
    xm,zm=meshgrid(x,z)
    D=v*aT
    c=AnLinFini(k,v,D,Lx,Lz,xm,zm)
    flux=D*(1-c[0,:])/(dz/2)/Lz
    flmoy=dxv*0
    xv=cumsum(dxv);prev=0.
    for i in range(nc):
        flmoy[i]=mean(flux[(x>prev)&(x<xv[i])]);prev = xv[i]
    return flmoy
