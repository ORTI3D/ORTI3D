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
from .geometry import *
from .geometryMesh import *

import os,time
from .opfoamKeywords import OpF
from .opfoamKeywords import OpT

class opfoam(unstructured):
    
    def __init__(self,core):
        self.core = core
        # unstructured provides the function buildMesh0 = usg
        
    def makeBC(self):
        '''
        create an array bcindx that gather the cell number and the 
        '''
        if self.core.getValueFromName('OpenFlow','MshType')>0:
            self.ncell_lay = len(self.carea)
            lbc,lnb = [],[]
            for k in self.dicFeats.keys():
                if k[:7]=='bc_cell':
                    arr = self.dicFeats[k]
                    lbc.extend(list(arr));lnb.extend([int(k[7:])]*len(arr))
            bcindx=zeros((len(lbc),2));bcindx[:,0]=lbc;bcindx[:,1]=array(lnb)+1
            self.nbc = len(unique(bcindx[:,1]))
            self.bcindx = bcindx  
        else: 
            bcindx = None
            grd = self.core.addin.getFullGrid()
            nx,ny,dx,dy = grd['nx'],grd['ny'],grd['dx'],grd['dy']
            self.ncell_lay = nx*ny
            dxm,dym = meshgrid(dx,dy)
            self.carea = ravel(dxm*dym)
            self.nbc = 4
       
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
            self.makeBC()
            self.opfMesh2Faces()
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
        self.carea = ravel(dxm*dym);self.nbc=4
        self.points,self.faces,self.fcup = points,faces,fcup
        self.elx,self.ely=points[fcup,0],points[fcup,1]
        return points,faces,bfaces,fcup

    def opfMesh2Faces(self):
        '''
        uses the "usg" type mesh to build the points and faces necessary for opf
        in M :itriangle, icell, ineighb,x,y
        '''
        M = self.M
        elx,ely,indx = self.elx,self.ely,self.indx
        u,id0 = unique(M[:,0],return_index=True)
        points = M[id0,3:5];npts = len(points);
        self.points=points
        nfc = len(M);nco = self.ncorn;nbc=int(amax(self.bcindx[:,0]))
        faces,i0c,i,ic = [],0,0,0
        fcup=[0]*self.ncell_lay
        while ic<nbc+1: # corners+bdy we keep only the internal faces
            ic = M[i,1]
            if ic not in self.bcindx[:,0]: 
                if M[i+1,1]>ic: i0c = i+1;
                i += 1
                continue
            if M[i+1,1]>ic: # last face in cell
                if M[i,2]>nbc:
                    faces.append([M[i,0],M[i0c,0],M[i,1],M[i,2],M[i,-1]])
                i0c = i+1;#print(ic,i0c)
            else :
                if M[i,2]>nbc:
                    faces.append([M[i,0],M[i+1,0],M[i,1],M[i,2],M[i,-1]])
            i += 1
        # now all the faces        
        i0c = i
        for i in range(nfc-1):
            ifc,ic,inb = M[i,:3].astype('int')
            if ic in self.bcindx[:,0]: 
                if M[i+1,1]>ic: i0c = i+1;
                continue
            if M[i+1,1]>ic : # last face of the current cell ic
                if (inb>ic)&(ifc!=M[i0c,0]): # store last face only if not done already
                    faces.append([M[i,0],M[i0c,0],ic,inb,M[i,-1]])
                # store fcup
                a = r_[M[i0c:i+1,0],M[i0c,0]]
                fcup[ic] = a.astype('int')
                i0c = i+1 # will be the starting point of th enext cell
            else:
                if (inb>ic)&(ifc!=M[i+1,0])&(inb>nbc): # store face only if not done already
                    faces.append([M[i,0],M[i+1,0],ic,inb,M[i,-1]])
        # last cell
        i+=1;inb=M[i,2]
        if (inb>ic)&(ifc!=M[i0c,0]): # store last face only if not done already
            faces.append([M[i,0],M[i0c,0],ic,inb,M[i,-1]])
        a = r_[M[i0c:i+1,0],M[i0c,0]]
        fcup[-1] = a.astype('int')
        faces = array(faces)
        lastcpt=[];lastsz=0
        #build the corners (from elx as they were done that way in voronoi)
        for i in range(nco):
            lc = self.bcindx[self.bcindx[:,1]==i+1,0].astype('int')
            pt0 = npts*1
            b = faces[faces[:,2]==i];
            b =b[b[:,-1].argsort(),:2].astype('int') # to keep the angle trigo round
            # a trick to get the pt where the loop is broken 
            if lastsz==1: ipt1 = lastcpt[i-1]
            else : ipt1 = pt0+2
            il=-1;u=list(b[0])
            for j in range(1,len(b)):
                if b[j,0]!=b[j-1,1] : 
                    il = j
                    for k in range(3-lastsz):
                        points = r_[points,c_[elx[i][j+1+k],ely[i][j+1+k]]];
                        u.append(npts);npts+=1
                    if lastsz==1: u.append(ipt1)
                    u.extend(b[j])
                else: u.append(b[j,1])
            if il==-1: 
                il=len(b)-1
                istart=where((elx[i]==points[b[il,1],0])&(ely[i]==points[b[il,1],1]))[0][0]
                for k in range(3-lastsz):
                    points = r_[points,c_[elx[i][istart+1+k],ely[i][istart+1+k]]];
                    u.append(npts);npts+=1
                if lastsz==1: u.append(ipt1)
                u.append(b[0,0])
            fcup[i] = array(u)
            if len(lc)>2: inb=lc[lc>nco-1][0]
            else : inb=i+1
            fc1 = array([u[u.index(pt0)-1],pt0,i,inb,0],ndmin=2)
            fc2 = array([pt0,pt0+1,i,-i-1,0],ndmin=2)
            if i==0: inb=-nco
            else :inb = -i
            fc3 = array([pt0+1,ipt1,i,inb,0],ndmin=2)
            faces = r_[faces,fc1,fc2,fc3]
            lastcpt.append(pt0)
            if len(lc)==2:lastsz=1
            else : lastsz=0
                
        # then bc lines
        for i in range(nco):
            lc = list(self.bcindx[self.bcindx[:,1]==i+1,0].astype('int'))
            cnt = 0
            for i0,j in enumerate(lc): # go along the boundary
                if j<nco: continue
                pt0 = npts*1
                b = faces[faces[:,2]==j];
                b =b[b[:,-1].argsort(),:2].astype('int') # to keep the angle trigo round
                # a trick to get the pt where the loop is broken 
                il=-1;u=list(b[0])
                for j1 in range(1,len(b)):
                    if b[j1,0]!=b[j1-1,1] : 
                        il = j1;u.extend([-2,-1]);u.extend(b[j1])
                        if i0!=len(lc)-1: 
                            points = r_[points,c_[elx[j][j1+1],ely[j][j1+1]]]
                            npts+=1
                    else:
                        u.append(b[j1,1])
                if il==-1: 
                    il=len(b)-1;u.extend([-2,-1]);u.append(b[0,0])
                    istart=where((elx[j]==points[b[-1,1],0])&(ely[j]==points[b[-1,1],1]))[0][0]
                    if i0!=len(lc)-1: 
                        points = r_[points,c_[elx[j][istart+1],ely[j][istart+1]]]
                        npts+=1
                    
                if (i0==len(lc)-1)&(i==nco-1): ipt1,inb1=lastcpt[0]+2,0 # the last cell back to 0
                elif i0==len(lc)-1: ipt1,inb1=lastcpt[i+1]+2,i+1; # lat cell in the row
                else: ipt1,inb1 = npts-1,j+1
                if cnt == 0: ipt2 = lastcpt[i];cnt=1
                else : ipt2 = pt0-1
                    
                fc1 = array([u[u.index(-2)-1],ipt1,j,inb1,0],ndmin=2)
                fc2 = array([ipt1,ipt2,j,-i-1,0],ndmin=2)
                faces = r_[faces,fc1,fc2]
                u[u.index(-2)]=ipt1;u[u.index(-1)]=ipt2
                fcup[j] = array(u)
                cnt=1
        #plot(points[fcup[j],0],points[fcup[j],1],'k')     
        #sort faces to have BC at the end
        faces = faces[:,:-1].astype('int')
        face1 = faces[faces[:,3]>=0]
        idr = where(face1[:,3]<face1[:,2])[0] # neighbour should not be smaller than owner
        face1[idr,:2]=face1[idr,1::-1]
        face1[idr,2:]=face1[idr,3:1:-1]
        face1 = face1[lexsort((face1[:,3],face1[:,2]))] #in lexsort col in reverse
        bfaces = faces[faces[:,3]<0]
        bfaces = bfaces[lexsort((-bfaces[:,3],bfaces[:,2]))] #in lexsort col in reverse        
        self.points,self.faces,self.bfaces,self.fcup = points,face1,bfaces,fcup
        return points,face1,bfaces,fcup
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
        if len(lspec)==0:
            lspec = [listE['i'][0]] # we need at least one species for sel out
        return ncomp,gcomp,lcomp,lspec
    
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
