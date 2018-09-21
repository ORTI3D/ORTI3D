# -*- coding: utf-8 -*-
"""
Created on Wed Aug 26 21:52:19 2015

@author: olive
"""
from geometry import *
import os

class Opgeo:
    
    def __init__(self, core):
        self.core = core
        
    def buildMesh(self,nlay=1):
        self.core.addin.mesh = self
        self.meshtype,self.threeD = 'mesh',False
        if self.core.dicval['OpgeoFlow']['domn.1'][0]==0: 
            self.meshtype = 'grid'
            nx,ny,xv,yv = getXYvects(self.core)
            self.nel = nx*ny
            self.nnod = (nx+1)*(ny+1)
            return None # rectangular
        dct = self.core.diczone['OpgeoFlow'].dic
        dicD = dct['domn.4']
        if dct.has_key('flow.5'): dicK = dct['flow.5'] # K hydraul
        else : dicK = {'name':[]}
        s = gmeshString(self.core,dicD,dicK)
        os.chdir(self.core.fileDir)
        f1 = open('gm_in.txt','w');f1.write(s);f1.close()
        bindir = self.core.baseDir+os.sep+'bin'+os.sep
        os.system(bindir+'gmsh gm_in.txt -2 -o gm_out.msh')
        #os.chdir(self.core.fileDir)
        f1 = open('gm_out.msh','r');s = f1.read();f1.close()
        nodes,elements = readGmshOut(s)
        self.nodes,nnod0 = nodes,len(nodes)
        self.nnod = nnod0*1
        s = self.arr2string(nodes); self.nodestring = s.replace('.  ','  ')
        nel,nc = shape(elements)
        elements[:,0]=arange(nel)
        elements[:,1]=elements[:,2] # this is the material number, which starts from 1 in gmsh and 0 in opgeo
        self.nel,self.nlay = nel,1
        self.elements = elements
        self.elx = nodes[elements[:,3:],1]
        self.ely = nodes[elements[:,3:],2]
        elements[:,2] = -100 # before the nodes number
        s = self.arr2string(elements); self.elementstring = s.replace('-100',' -1 tri')
        self = createTriangul(self)
        if self.core.addin.getDim() == '3D':
            self.threeD = True
            self.nodes = self.nodes3D(nodes)
            s = self.arr2string(self.nodes); self.nodestring = s.replace('.  ','  ');#print self.nodestring
            self.elements = self.elts3D(elements,nnod0)
            s = self.arr2string(self.elements); self.elementstring = s.replace('-100',' -1 pris')            
        return 
        
    def getCenters(self):    
        if self.core.dicval['OpgeoFlow']['domn.1'][0]==0:
            return getXYmeshCenters(self.core,'Z',0)
        else :
            return self.elcenters
            
    def getNumber(self,typ):
        if typ=='elements': return self.nel
        else : return self.nnod
            
    def arr2string(self,arr):
        s=''
        nr,nc = shape(arr)
        for i in range(nr):
            s += str(int(arr[i,0]))+' '
            for j in range(1,nc):
                s += str(arr[i,j])+' '
            s += '\n'
        return s
        
    def arr2string1(self,arr):
        s=''
        nr,nc = shape(arr)
        for i in range(nr):
            s += str(arr[i,0])+' '
            for j in range(1,nc):
                s += str(arr[i,j])+' '
            s += '\n'
        return s
    
    def nodes3D(self,nodes):
        '''builds 3d nodes from 2d ones'''
        zb = makeZblock(self.core)
        self.core.Zblock = zb
        nlay,nnod = shape(zb)
        self.nlay = nlay
        nod1 = self.nodes*1
        nod2 = zeros((nlay*nnod,4))
        for il in range(nlay):
            nod2[nnod*il:nnod*(il+1),:] = c_[nod1[:,:3],zb[il:il+1].T]
        nod2[:,0] = range(nnod*nlay)
        self.nnod = nnod*nlay
        return nod2
        
    def elts3D(self,elements,nnod0):
        '''builds 3D elements from 2d ones
        col 0 : nb, 1 : material, 2 : -100 to be replaced, 3...node nb
        up to now nlay = nmedia!!!!!
        there are node layers (=phys layer +1) and elts layers (=physical layers)
        seem to be correct'''
        nlay,nel = self.nlay,self.nel;print nel, nnod0
        elts1 = self.elements
        elts2 = zeros(((nlay-1)*self.nel,9)).astype('int')
        elts2[:,2] = -100
        for il in range(nlay-1):
            elts2[nel*il:nel*(il+1),3:] = c_[elts1[:,3:]+nnod0*il,elts1[:,3:]+nnod0*(il+1)]
        nel3 = nel*(nlay-1)
        elts2[:,0] = range(nel3)
        # need to retrieve the true material nb per layer as gmsh built it in 2D
        for il in range(nlay-1):
            zindx = zone2mesh(self.core,'OpgeoFlow','flow.5',il,iper=0,loc='elements',val='nb')
            elts2[nel*il:nel*(il+1),1] = zindx
        self.nel = nel3
        return elts2