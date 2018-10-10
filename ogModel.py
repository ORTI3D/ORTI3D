# -*- coding: utf-8 -*-
"""
Created on Wed Aug 26 21:52:19 2015

@author: olive
"""
from .geometry import *
import matplotlib.tri as mptri
import os

class OgModel:
    
    def __init__(self, core):
        self.core = core
        
    def buildMesh(self,nlay=1):
        if self.core.dicval['OpgeoFlow']['domn.1'][0]==0: 
            nx,ny,xv,yv = getXYvects(self.core)
            self.nel = nx*ny
            self.nnod = (nx+1)*(ny+1)
            return None # rectangular
        dicz = self.core.diczone['OpgeoFlow'].dic['domn.4']
        s = gmeshString(dicz)
        os.chdir(self.core.fileDir)
        f1 = open('gm_in.txt','w');f1.write(s);f1.close()
        bindir = self.core.baseDir+os.sep+'bin'+os.sep
        os.system(bindir+'gmsh gm_in.txt -2 -o gm_out.msh')
        f1 = open('gm_out.msh','r');a = f1.read();f1.close()
        # nodes
        b = a.split('$Nodes')[1]
        c = b.split('$EndNodes')[0]
        c1 = c.split('\n')
        self.nnod = int(c1[1])        
        c2 = [x.split() for x in c1[2:-1]];#print len(c2)
        arr = array(c2,dtype='float')
        arr[:,0] = arr[:,0]-1 # node number start from 0
        self.nodes = arr
        s1 = self.arr2string(arr)
        self.nodestring = s1.replace('.  ','  ')
        # elements
        b = a.split('$Elements')[1];#print len(b)
        c = b.split('$EndElements')[0];#print len(c)
        c1 = c.split('\n')
        c2 = [x.split() for x in c1[2:] if len(x.split())==8];#print len(c2)
        arr = array(c2,dtype='int')
        arr = arr-1
        arr = arr[:,[0,4,5,6,7]] # 2nd column is useless
        nel,nc = shape(arr)
        arr[:,0]=arange(nel)
        self.nel = nel
        self.elements = arr
        arr[:,1] = -100
        s = self.arr2string(arr)
        s = s.replace('-100','0 -1 tri')
        self.elementstring = s
        self.elx = take(self.nodes[:,1],self.elements[:,-3:])
        self.ely = take(self.nodes[:,2],self.elements[:,-3:])
        self.elcenters = [mean(self.elx,axis=1),mean(self.ely,axis=1)]
        return 
        
    def getCellCenters(self):    
        if self.core.dicval['OpgeoFlow']['domn.1'][0]==0:
            return getXYmeshCenters(self.core,'Z',0)
        else :
            return self.elcenters
            
    def getNumberOfCells(self):
        return self.nel
            
    def arr2string(self,arr):
        s=''
        nr,nc = shape(arr)
        for i in range(nr):
            s += str(int(arr[i,0]))+' '
            for j in range(1,nc):
                s += str(arr[i,j])+' '
            s += '\n'
        return s