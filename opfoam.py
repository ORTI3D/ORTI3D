# -*- coding: utf-8 -*-
"""
Created on Sun May 17 10:55:23 2020

@author: olivier
geometry of the openfoam mesh, we use here mesh from modflowUsg
Simply said upper-triangular order is that the order of faces corresponds to the order of the cells they connect.
- take all (higher numbered) cells connected to a cell.
- sort them according to their cell label
- the faces connecting to those cells should now also be ordered.
renumberMesh will do this for you. (IT Works)

The algorithm to assign indices to a cell is as follows:
* Partition faces (yes, faces) in boundaries
* Maintain a cell index counter, start it at zero.
* Go through all these faces: find the face its owner (i.e. the cell that owns that face).
If that cell has not been assigned a cell index yet, assign it a cell index.
Increase that cell index counter. 
"""
# a file with point coordinates nbofpts ((x0 y0 z0) (x1 y1 z1)...)
# faces file nboffaces (nbpts (pt0 pt1 pt2..) nbpts (pt0 pt1...) ...)
# owner file  al ist of the cell that owns the considered face nb (0 1 2 0)
# neighbour 
from scipy import zeros,ones,array,arange,r_,c_,around,argsort,unique,cumsum,where,shape,\
    amin,amax,mod,lexsort
from geometry import *
from geometryMesh import *

import os,time
from opfoamKeywords import OpF
from opfoamKeywords import OpT

class opfoam(unstructured):
    
    def __init__(self,core):
        self.core = core
        # unstructured provides the function buildMesh0 = usg
       
    def buildMesh(self,opt):
        self.buildMesh0('OpenFlow',opt)
        mshType = self.core.getValueFromName('OpenFlow','MshType')
        if mshType == 0: # regular cells
            self.opfRect()
        if mshType == 1: # nested
            #msh = myRect(self,nodes,elements)
            #msh.calc()
            pass
        if mshType == 2: # triangle
            pass
        if mshType == 3: # voronoi
            msh = myVor(self)
            msh.transformVor(self.points,self.elts,self.dcoo1,self.dicD,self.dicFeats) # OA 15/8/20
            xn,yn = self.nodes[:,0],self.nodes[:,1];print('voronoi made')
            self.trg = mptri.Triangulation(xn,yn) #,triangles=elts) 
            self.makeBC('OpenFlow')
            self.getPointsFaces()
        self.addMeshVects()
        Zblock = makeZblock(self.core)
        nlay = getNlayers(self.core)
        if self.core.addin.getDim()=='3D': 
            thick = Zblock[:-1]-Zblock[1:]
            self.core.addin.get3D();print('start 3D')
        else :
            thick = Zblock[0]-Zblock[-1]
            #self.fahl = [array(lg)*thick[i] for i,lg in enumerate(self.fahl)]

    
    def opfRect(self):
        grd = self.core.addin.getFullGrid()
        dx,dy = array(grd['dx']), array(grd['dy']);
        nx,ny,x0,y0 = grd['nx'],grd['ny'],grd['x0'],grd['y0']
        xv,yv = r_[x0,x0+cumsum(dx)],r_[y0,y0+cumsum(dy)]
        npt = (nx+1)*(ny+1)
        l = [];fc = [];cn=[]
        # build the outside boundaries (faces without neighbour)
        for i in range(nx): l.append([i,i+1,i,-3])# botm
        for j in range(ny): l.append([j*(nx+1)+nx,(j+1)*(nx+1)+nx,(j+1)*nx-1,-2]);# right
        for i in range(nx,0,-1): l.append([ny*(nx+1)+i,ny*(nx+1)+i-1,nx*(ny-1)+i-1,-1]);#top
        for j in range(ny,0,-1): l.append([j*(nx+1),(j-1)*(nx+1),(j-1)*nx,-4]);# left
        # then the faces in the center
        for j in range(ny):
            for i in range(nx-1):
                ic=j*nx+i;
                l.append([j*(nx+1)+i+1,(j+1)*(nx+1)+i+1,ic,ic+1]) # right
        for j in range(ny-1):
            for i in range(nx):      
                ic=j*nx+i;
                l.append([(j+1)*(nx+1)+i+1,(j+1)*(nx+1)+i,ic,ic+nx]) # top
        # set the upper face for all cells
        for j in range(ny):
            for i in range(nx):
                fc.append([j*(nx+1)+i,j*(nx+1)+i+1,(j+1)*(nx+1)+i+1,(j+1)*(nx+1)+i,j*(nx+1)+i])
        fc0 = array(l,ndmin=2);fcup = array(fc,ndmin=2)
        faces = fc0[fc0[:,3]>=0];bfaces = fc0[fc0[:,3]<0]
        ag = argsort(bfaces[:,3])
        self.bfaces = bfaces[ag[-1::-1]]
        xm,ym = meshgrid(xv,yv)
        points = c_[reshape(xm,(npt,1)),reshape(ym,(npt,1))]
        dxm,dym = meshgrid(dx,dy)
        self.nbc=4
        if self.core.addin.getDim() in ['2D','3D']:
            self.carea = ravel(dxm*dym);
        else : #xsection
            dz = self.core.dicval['OpenFlow']['dis.6'][0]-self.core.dicval['OpenFlow']['dis.7'][0]
            self.carea = ravel(dxm*dz)
        self.points,self.faces,self.fcup = points,faces,fcup
        self.elx,self.ely=points[fcup,0],points[fcup,1]
        return points,faces,bfaces,fcup

        '''
        for i in range(len(faces)): 
            plot(points[faces[i,:2],0],points[faces[i,:2],1])
        ic=57;fc=faces[faces[:,2]==ic,:2] 
        for a in fc: plot(points[a,0],points[a,1])
        for j in range(len(fcup)): 
            plot(points[fcup[j],0],points[fcup[j],1])
        '''
            
    def opfDictFromUsg(self):
        '''get some values form usg to openfoam dicts'''
        l0 = ['disu.2','disu.3','disu.9','bas.5']
        l1 = ['dis.1','dis.2','dis.5','head.1']
        for i in range(len(l0)):
            self.core.dicval['OpenFlow'][l1[i]] = self.core.dicval['Modflow'][l0[i]]
            self.core.dictype['OpenFlow'][l1[i]] = self.core.dictype['Modflow'][l0[i]]
        #domain
        self.core.diczone['OpenFlow'].dic['dis.2'] = self.core.diczone['Modflow'].dic['disu.3']        
        # fixed heads are for the zones of bas.5
        if 'bas.5' in self.core.diczone['Modflow'].dic.keys():
            self.core.diczone['OpenFlow'].dic['head.2'] = self.core.diczone['Modflow'].dic['bas.5']        
        l0 = ['disu.7','disu.8','lpf.8','lpf.9','wel.1','drn.1','ghb.1','rch.2']
        l1 = ['dis.6','dis.7','khy.2','khy.3','wel','drn','ghb','rch']
        for i in range(len(l0)):
            self.core.dicval['OpenFlow'][l1[i]] = self.core.dicval['Modflow'][l0[i]]
            self.core.dictype['OpenFlow'][l1[i]] = self.core.dictype['Modflow'][l0[i]]
            if l0[i] in self.core.diczone['Modflow'].dic.keys():
                self.core.diczone['OpenFlow'].dic[l1[i]] = self.core.diczone['Modflow'].dic[l0[i]]
            if l0[i] in self.core.dicarray['Modflow'].keys():
                self.core.dicarray['OpenFlow'][l1[i]] = self.core.dicarray['Modflow'][l0[i]]
        # transport
        self.core.dicval['OpenTrans'],self.core.dictype['Optrans'] = {},{}
        #self.core.diczone['OpenTrans'] = dicZone(self.core,'Optrans')
        l0 = ['bct.2','bct.3','bct.20','pcb.2']
        l1 = ['cactiv','poro','cinit','cfix']
        for i in range(len(l0)):
            self.core.dicval['OpenTrans'][l1[i]] = self.core.dicval['MfUsgTrans'][l0[i]]
            self.core.dictype['OpenTrans'][l1[i]] = self.core.dictype['MfUsgTrans'][l0[i]]
            if l0[i] in self.core.diczone['MfUsgTrans'].dic.keys():
                self.core.diczone['OpenTrans'].dic[l1[i]] = self.core.diczone['MfUsgTrans'].dic[l0[i]]
        l0=['ph.3','ph.4','ph.5']
        l1=['sinit','sfix','srch']
        for i in range(len(l0)):
            self.core.dicval['OpenChem'][l1[i]] = self.core.dicval['Pht3d'][l0[i]]
            self.core.dictype['OpenChem'][l1[i]] = self.core.dictype['Pht3d'][l0[i]]
            if l0[i] in self.core.diczone['MfUsgTrans'].dic.keys():
                self.core.diczone['OpenChem'].dic[l1[i]] = self.core.diczone['Pht3d'].dic[l0[i]]

    def findSpecies(self,core):
        '''find the components and species for phreeqc'''
        listE = core.addin.chem.getDictSpecies()
        listS = listE['i'];listS.extend(listE['k']);listS.extend(listE['kim'])
        listS.sort()
        nbs,ncomp,gcomp,lcomp,lspec,i = len(listS),0,len(listE['g']),[],[],0
        while i<nbs:
            if listE['i'][i] == 'O(0)': lspec.append('O(0)')
            if listE['i'][i] in ['H(1)','H(+1)']: lspec.append('H(1)')
            elif '(' in listE['i'][i]: 
                lspec.append(listE['i'][i])
                lcomp.append(listE['i'][i]);ncomp +=1
                if '(' in listE['i'][i+1]: 
                    lspec.append(listE['i'][i+1]);i +=1
            else :
                lcomp.append(listE['i'][i]);ncomp += 1
            i += 1
        ncomp += 2 # there is pH and pe but we need to have H20,H,0,charge and not pH,pe
        return ncomp,gcomp,lcomp,listE['g'],lspec
    
# **************** potential link with USG faces, not used
    '''
    def facesLinkUsg(self):
        #here we make the link between faces of opfreak and Usg in order
        #to recover the fluxes if necessary, it needs to find the Usg face number
        #in the fcup freak list
        pts = self.points
        for ic in range(self.ncell):
            p0 = pts[fcup[ic]][0];#n,a = shape(p0)
            x,y = mesh.elx[ic],mesh.ely[ic]
            i0 = where((x==p0[0])&(y==p0[1]))[0][0]
    '''
