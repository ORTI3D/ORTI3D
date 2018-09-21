# -*- coding: utf-8 -*-
"""
Created on Sat Sep 20 08:55:28 2014

@author: olive
"""
try :
    from fipy import *
    from fipy.terms.cellTerm import CellTerm
except ImportError : 
    pass
from array import array as arr2
from geometry import *
import os

class FpModel:
    """A class that contains the mesh and equations for fipy"""
    def __init__(self, core):
        self.core = core
        
    def buildMesh(self,nlay=1):
        """a method to build a rectangular or triangular mesh, from cell width 
        or sets of points (using Gmesh)"""
        if self.core.dicval['FipyFlow']['domn.1'][0]==0: # rectangular
            #print 'rect'
            g = self.core.addin.getFullGrid()
            dx, dy = g['dx'], g['dy']
            self.nx, self.ny, self.nlay = len(dx),len(dy),nlay
            self.mesh = Grid2D(dx = dx, nx = self.nx,dy=dy,ny=self.ny)
            self.ncells = self.nx*self.ny
        else :
            #print 'tri'
            dicz = self.core.diczone['FipyFlow'].dic['domn.4']
            s = gmeshString(dicz)
            self.mesh=Gmsh2D(s)
            self.ncells = self.mesh.getNumberOfCells()
        
    def run(self,fDir,fName,option):
        self.fDir,self.fName = fDir,fName
        self.initModel(option)
        self.setIC(option)
        self.setBC()
        ttable = self.core.makeTtable();#print 'mfw',self.ttable
        self.modelSolve(ttable['tlist'],10,option)

    def initModel(self,option):
        """ a method to initiate the major elements for fipy"""
        poro = self.core.getValueFromName('FipyFlow','F_PORO')
        Rf = 1.
        self.Ss=1
        self.H = CellVariable(mesh=self.mesh,value=10.,hasOld=True)
        self.K = CellVariable(mesh=self.mesh,value=1.)
        self.Ss = CellVariable(mesh=self.mesh,value=1e-6)
        self.Q = CellVariable(mesh=self.mesh,value=0.,hasOld=True) # hasOld does not seem to be usefull
        self.V = -self.K * self.H.grad/poro/Rf
        self.C = CellVariable(mesh=self.mesh,value=0.,hasOld=True) 

    def setIC(self,option):
        """set the initial conditions"""
        H_init = self.core.getValueLong('FipyFlow','flow.1',0)
        thick = self.core.getValueLong('FipyFlow','domn.5',0) - self.core.getValueLong('FipyFlow','domn.6',0)
        K0 = self.core.getValueLong('FipyFlow','flow.3',0)#*thick
        mtype = self.core.addin.getModelType()
        if mtype == 'confined':
            Ss0 = self.core.getValueLong('FipyFlow','flow.5',0)
        else :
            Ss0 = self.core.getValueLong('FipyFlow','flow.6',0)
        self.Ss.setValue(Ss0)
        self.K.setValue(K0)
        self.H.setValue(H_init)
        C0 = self.core.getValueLong('FipyTrans','trans.6',0)
        self.C.setValue(C0)
    
    def setBC(self):
        """sets the boundary conditions over the domain"""
        self.setBCline(self.H,'flow.2','fixed','Flow')
        self.setBCline(self.C,'trans.7','fixed','Trans')
        self.setBCline(self.C,'trans.9','grad','Trans')
        
    def setBCline(self,variable,line,typ,option):
        if line not in self.core.diczone['Fipy'+option].dic.keys():
            return
        zlist = self.core.diczone['Fipy'+option].dic[line] # H_fixed
        for iz in range(len(zlist['coords'])):
            coo = zlist['coords'][iz]
            val = float(zlist['value'][iz])           
            nn = self.zone2faces(coo,0);#print iz,coo,nn
            if typ=='fixed':
                variable.constrain(val,where = nn)
            elif typ == 'grad':
                grad = self.V.faceValue-(self.V*variable).faceValue
                variable.faceGrad.constrain(grad,where=nn)
                
    
    def modelSolve(self,tlist,nsubsteps,option):
        """solves the current problem for head """        
        f1=open(self.fDir+os.sep+self.fName+'.fhd','wb')
        f2=open(self.fDir+os.sep+self.fName+'.ucn','wb')
        #reach steady state
        (ImplicitDiffusionTerm(coeff = self.K.faceValue) + self.Q).solve(var=self.H)

        eqFlow = TransientTerm(coeff=self.Ss) == \
            ImplicitDiffusionTerm(coeff = self.K.faceValue) + self.Q
        poro, Rf = 0.25, 1.
        
        if option == 'Trans':
            Rf, v0= 1.,1.
            Dff = self.core.getValueFromName('FipyTrans','F_DIFF')
            aL = self.core.getValueFromName('FipyTrans','F_AL')
            aT = self.core.getValueFromName('FipyTrans','F_AT');print aL,aT
            Dsp = CellVariable(mesh=self.mesh,value=[[aL*v0,0.],[0.,aT*v0]])

            eqTrans = TransientTerm() == \
                ImplicitDiffusionTerm(coeff=[Dsp]) - \
                VanLeerConvectionTerm(coeff=self.V) #best for sharp gradient
                #PowerLawConvectionTerm(coeff=self.V)

        tprec,dt2=tlist[0],1e-5
        for iper,t in enumerate(tlist[1:]):
            dt,dt1 = 0,t-tprec
            #dt2 = dt1/nsubsteps;print t,dt2
            #for i in range(nsubsteps):
            while dt<dt1:
                self.H.updateOld()
                q = zone2mesh(self.core,'FipyFlow','flow.7',iper=iper);#print 'fp l68',q
                self.Q.setValue(q/self.mesh.cellVolumes)
                eqFlow.solve(var = self.H,dt = dt2)
                if option == 'Flow': continue
                # to get dispersion tensor
                vx,vy = self.V.value;#print vx[:2]
                v = sqrt(vx**2+vy**2)
                dxx = (aL*vx**2+aT*vy**2)/v+Dff
                dxy = (aL-aT)*vx*vy/v+Dff
                dyy = (aL*vy**2+aT*vx**2)/v+Dff
                Dsp = [[dxx,dxy],[dxy,dyy]] # [[D1,D2]] does not work
                self.C.updateOld()
                #eqTrans.solve(var = self.C,dt = dt2)
                i,res = 1,1.
                while res > 1.0E-4: # below 1e-4 fails
                    res = eqTrans.sweep(dt=dt2,var=self.C);i+=1
                dt += dt2;print t,i,dt
                if i<4: dt2 = dt2 * 2
                elif i<7 : dt2 = dt2 * 1.05
                else : dt2 = dt2*0.95
            tprec=t
            # save the values for each time step
            b = arr2('f');b.extend(self.H.value);b.tofile(f1)
            if option =='Trans':
                b = arr2('f');b.extend(self.C.value);b.tofile(f2)                
        f1.close();f2.close()

    def zone2faces(self,coords,medialist):
        #up to now just for one line with 2 points in x or y dir
        # finds the x,y,z of a zone
        if type(medialist)==type(1): medialist= [medialist] # sometime there is just one layer
        x,y = zip(*coords)
        # returns distance from cell faces to eq of the line
        xf,yf = self.mesh.faceCenters
        if x[1]!=x[0]:
            a = (y[1]-y[0])/(x[1]-x[0])
            b = y[1]-a*x[1]
            dst = abs((yf-a*xf-b)/a)
        else :
            dst = abs(xf - x[0])
        # finds faces from distance
        tol = (max(xf)-min(xf))/200
        xmn,xmx,ymn,ymx = min(x)-tol,max(x)+tol,min(y)-tol,max(y)+tol
        nn0 = (dst<tol)&(xf>xmn)&(xf<xmx)&(yf>ymn)&(yf<ymx)
#        if self.TriDim :
#            laylist = []
#            for m in medialist: laylist.extend(media2layers(self.core,m))
#            nn2d = (self.nx+1)*(self.ny+1)
#            #creates an array contiang nodes for all layers
#            for i,lay in enumerate(laylist):
#                if i==0: nn = list(nn0+nn2d*lay)
#                else : nn.extend(list(nn0+nn2d*lay))
#        else :
#            nn = list(nn0)
        return nn0
