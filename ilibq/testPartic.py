# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 17:36:24 2017

@author: oatteia
"""
from numpy import ones,array,r_,sign,exp
mdir = 'D:\\ipht3d\\libdev'
sys.path.append(mdir)
from ilibq import core
md= core.Core()

from scipy.stats import norm

fdir='D:\\ipht3d\\exv2\\test2d';fname='Test2d'
md.openModel(fdir,fname)

inst=md.addin.fit.fitter
dic_options,obs_value = {'type':'Head','zvalue':0,'dispersion':[10,1]},None
inst.H_bc = (md.getValueLong('Modflow','bas.5',0)*(md.getValueLong('Modflow','bas.3',0)<0))[0]
inst.calculate(dic_options,obs_value)

q = inst.Q
x0,x1,y0,y1 = [float(a) for a in inst.domain]
dx,dy = inst.dx,inst.dy
aL,aT  =  10,1
poro = 0.25
vx,vy = q[0]/poro, q[1]/poro
vx= c_[vx,vx[:,-1:]];vy =r_[vy,vy[-1:,:]]
#from geometry import *
#xv,yv=getXYmeshCenters(md,'Z',0)
yp0=array([0,.08,.2,.3,.4,.5,.58,.65,.7,.75,.8,.85,0.89,.92,.94,.96,0.985,0.999])
xy=[(0.5,30),(0.5,25)];x, y = zip(*xy)
data = (x0,dx,y0,dy,aT,aL,0.);# print data
ypi=r_[-yp0[-1:0:-1],yp0]/2.
xpi=ypi*0.;
xpin = (x[0]+x[1])/2.+ypi*(x[0]-x[1]) # 
ypin = (y[0]+y[1])/2.+ypi*(y[0]-y[1]);
nr,nc = shape(inst.K)
clim=ones((nr,nc));#clim[1:-1,1:-1]=1
alph0=-sign(ypi)*norm.isf(.5+abs(ypi),0.,.25); # distribution of starting point along the source line
alph0=(exp(abs(alph0))-1)*sign(alph0)
import rflowC1b
[xp0,yp0,tp0,cua,cub]=rflowC1b.calcTube1(data,alph0,vx,vy,xpin,ypin,clim);
[xp1,yp1,tp1,cu1,cu2,V,Large,Cinf] = inst.calcConc(xp0,yp0,tp0,cua,cub);
for i in range(34):plot(xp0[:,i],yp0[:,i])