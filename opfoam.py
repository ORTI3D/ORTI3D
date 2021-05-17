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
    amin,amax,mod
from ilibq.geometry import *

import os

class opfoam:
    def __init__(self,md,mesh):
        self.mesh = mesh
        if md.getValueFromName('Modflow','MshType')>0:
            self.ncell_lay = len(mesh.carea)
            lbc,lnb = [],[]
            for k in mesh.dicFeats.keys():
                if k[:7]=='bc_cell':
                    arr = mesh.dicFeats[k]
                    lbc.extend(list(arr));lnb.extend([int(k[-1])]*len(arr))
            bcindx=zeros((len(lbc),2));bcindx[:,0]=lbc;bcindx[:,1]=array(lnb)+1
            self.area = mesh.carea
            self.nbc = len(unique(bcindx[:,-1]))
        else: 
            bcindx = None
            grd = md.addin.getFullGrid()
            nx,ny,dx,dy = grd['nx'],grd['ny'],grd['dx'],grd['dy']
            self.ncell_lay = nx*ny
            dxm,dym = meshgrid(dx,dy)
            self.area = ravel(dxm*dym)
            self.nbc = 4
        self.bcindx = bcindx
        
    def opfMeshFromUsg(self,md):
        # determiner le num des bcs
        mesh,bcindx,dim = self.mesh,self.bcindx,md.addin.getDim()
        # creer un array de faces à partir des cells en conservant ds colonne 2 le num de cellule
        # le num de cellule est -nb pour les line nb de BC
        '''
        for ic in mesh.dicFeats['bc_cell0']: plot(mesh.elx[ic],mesh.ely[ic],'g')
        '''
        lnghb=[];ncell=self.ncell_lay # we use only the first layer
        #each array has columns first pt, 2nd pt, cell nb, neighbour, x, y 
        for ic in range(ncell):
            nf=len(mesh.elx[ic])-1
            if ic in bcindx[:,0]: # in BC
                cn1 = zeros(nf);idx=mesh.indx[ic]; # idx for BC
                bcv=bcindx[bcindx[:,0]==ic,1] # for BC
                try :
                    len(idx)
                    if ic==0:
                        cn1[idx[0]],cn1[idx[1]] = -min(bcv),-max(bcv)
                    else :
                        cn1[idx[0]],cn1[idx[1]] = -max(bcv),-min(bcv)
                except TypeError:
                    cn1[idx] = -bcv
                k = 0 # adding the neighb if it is not a bdy 
                for j in range(nf):
                    if cn1[j]==0: 
                        cn1[j] = mesh.cneighb[ic][k];k+=1            
            else :
                if dim == '3D': cn1 = mesh.cneighb[ic][:-1] # dont take the bottom
                else : cn1 = mesh.cneighb[ic]
            lnghb.append(cn1)

        # specific run for the first cell (no neighbour)
        ic=0
        nf = len(lnghb[0])
        fc = zeros((nf,4))
        fc[:,0],fc[:-1,1] = list(range(nf)),list(range(1,nf))
        fc[:,2],fc[:,3] = [0]*nf,lnghb[0]
        points = zeros((nf,2))
        points[:,0],points[:,1] = mesh.elx[0][:-1],mesh.ely[0][:-1]
        pcnt = nf-1
        a = list(range(nf));a.append(0);fcup = [a]

        for ic in range(1,ncell):
            ng = lnghb[ic];nf = len(ng)
            m=fc[fc[:,3]==ic]
            lown = m[:,2];fup=[]
            # look at the first point
            if ng[-1] in lown:
                p_prec = m[lown==ng[-1],0][0] # 0 becasue faces was done in reverse order
            else :
                if ng[0] not in lown: # don't add point if next face was done
                    pcnt += 1; p_prec = pcnt;
                    points = r_[points,array([mesh.elx[ic][0],mesh.ely[ic][0]],ndmin=2)]
                else :
                    p_prec = m[lown==ng[0],1][0]
            fup.append(int(p_prec))
            for jf in range(nf) : # loop on faces    
                if ng[jf] in lown : # the face has been done
                    p_prec = m[lown==ng[jf],0][0] # nothing else for faces
                    fup.append(int(p_prec))
                else :
                    if jf<nf-1 and ng[jf+1] in lown: # next face has been done
                        pt2 = m[lown==ng[jf+1],1][0];fup.append(int(pt2))
                        addf = [p_prec,pt2,ic,ng[jf]]
                        fc = r_[fc,array(addf,ndmin=2)]
                    elif jf==nf-1: #next face (0) has been done
                        if ng[0] in lown:
                            pt2 = m[lown==ng[0],1][0];fup.append(int(pt2))
                            addf = [p_prec,pt2,ic,ng[jf]]
                            fc = r_[fc,array(addf,ndmin=2)]   
                        else : # new face, the last one so don't add point
                            addf = [p_prec,fup[0],ic,ng[jf]];fup.append(fup[0])
                            fc = r_[fc,array(addf,ndmin=2)]   
                            p_prec = fup[0]                             
                    else : # fully new face, not the last one
                        pcnt +=1
                        addf = [p_prec,pcnt,ic,ng[jf]];fup.append(pcnt)
                        fc = r_[fc,array(addf,ndmin=2)]            
                        points = r_[points,array([mesh.elx[ic][jf+1],mesh.ely[ic][jf+1]],ndmin=2)]
                        p_prec = pcnt
            fcup.append(fup)
            #print(ic,amax(faces[:,:2]),shape(points)[0])
        #sort faces to have BC at the end
        faces = fc[fc[:,3]>=0]
        bfaces = fc[fc[:,3]<0]
        ag = argsort(bfaces[:,3])
        bfaces = bfaces[ag[-1::-1]]
        '''
        nf=shape(faces)[0]
        for kf in range(nf): 
            plot(points[faces[kf,:2].astype('int'),0],points[faces[kf,:2].astype('int'),1])
        nf=shape(bfaces)[0]
        for kf in range(nf): 
            plot(points[bfaces[kf,:2].astype('int'),0],points[bfaces[kf,:2].astype('int'),1])
        for ic in range(ncell):
            plot(points[fcup[ic],0],points[fcup[ic],1])
        '''
        self.points,self.faces,self.fcup = points,faces,fcup
        return points,faces,bfaces,fcup
    
    def opfMeshReg(self,md):
        grd = md.addin.getFullGrid()
        dx,dy = array(grd['dx']), array(grd['dy']);
        nx,ny,x0,y0 = grd['nx'],grd['ny'],grd['x0'],grd['y0']
        xv,yv = r_[x0,x0+cumsum(dx)],r_[y0,y0+cumsum(dy)]
        ncell = nx*ny;npt = (nx+1)*(ny+1)
        ic= 0
        l = [];fc = []
        for j in range(ny):
            for i in range(nx):
                if j==0: l.append([j*(nx+1)+i,j*(nx+1)+i+1,ic,-3]) # botm
                if i<nx-1: l.append([j*(nx+1)+i+1,(j+1)*(nx+1)+i+1,ic,ic+1]) # right
                if i==nx-1: l.append([j*(nx+1)+i+1,(j+1)*(nx+1)+i+1,ic,-2]) # right
                if j<ny-1: l.append([(j+1)*(nx+1)+i+1,(j+1)*(nx+1)+i,ic,ic+nx]) # top
                if j==ny-1: l.append([(j+1)*(nx+1)+i+1,(j+1)*(nx+1)+i,ic,-1]) # top
                if i==0: l.append([(j+1)*(nx+1)+i,j*(nx+1)+i,ic,-4])
                fc.append([j*(nx+1)+i,j*(nx+1)+i+1,(j+1)*(nx+1)+i+1,(j+1)*(nx+1)+i,j*(nx+1)+i])
                ic += 1
        fc0 = array(l,ndmin=2);fcup = array(fc,ndmin=2)
        faces = fc0[fc0[:,3]>=0];bfaces = fc0[fc0[:,3]<0]
        ag = argsort(bfaces[:,3])
        bfaces = bfaces[ag[-1::-1]]
        xm,ym = meshgrid(xv,yv)
        points = c_[reshape(xm,(npt,1)),reshape(ym,(npt,1))]
        return points,faces,bfaces,fcup
    
# **************** potential link with USG faces, not used
    def facesLinkUsg(self):
        '''here we make the link between faces of opfreak and Usg in order
        to recover the fluxes if necessary, it needs to find the Usg face number
        in the fcup freak list
        '''
        pts = self.points
        for ic in range(self.ncell):
            p0 = pts[fcup[ic]][0];#n,a = shape(p0)
            x,y = mesh.elx[ic],mesh.ely[ic]
            i0 = where((x==p0[0])&(y==p0[1]))[0][0]

#***************************************************************
#************************* in geometry  ***********************

def cellsUnderPoly_new(core,dicz,media,iz):
    ''' a new version fo cellunder poly because the old one was too slow
    it starts from the first points of poly, search in the neighbouring cells
    which one is under the currnet segment and followup to next segment
    and then end of poly
    '''
    mesh = core.addin.mesh
    xc,yc,idcell = mesh.elcenters[:,0],mesh.elcenters[:,1],mesh.idc
    elxa, elya = array(mesh.elxa),array(mesh.elya)
    poly = dicz['coords'][iz]
    if len(poly)==1: #OA 18/12/20
        x,y = poly[0];dst=(x-xc)**2+(y-yc)**2
        indx=where(dst==amin(dst))[0]
        return indx,dicz['value'][iz]
    lcoefs=lcoefsFromPoly(poly)
    z = 0 # to bemodifiedfurther
    x,y = zip(*coo)
    d = sqrt((x[0]-xc)**2+(y[0]-yc)**2)
    ic = where(amin(d)==d)[0][0]
    iseg = 0 # index of the segments
    # find in the neighb cells which one is below the segments
    lpt = [ic]
    while True:
        print(iseg)
        a,b = lcoefs[0,iseg],lcoefs[1,iseg]
        cnb = list(mod(mesh.cneighb[ic],mesh.ncell_lay));cnb.remove(ic)
        for ing in cnb:
            pos = sign(a*mesh.elx[ing]+b*mesh.ely[ing]-1) # posit. vs the line
            seg0 = sign(b*xc[ing] - a*yc[ing] - b*x[iseg] +a*y[iseg])#bx’-ay’-bx0+ay0
            seg1 = sign(b*xc[ing] - a*yc[ing] - b*x[iseg+1] +a*y[iseg+1])#bx’-ay’-bx0+ay0
            if (max(pos)!=min(pos))*(seg0 != seg1)*1: break
        if ing==cnb[-1] : 
            if iseg == len(x)-2: break
            else : iseg += 1
        else : 
            lpt.append(ing) ;ic = ing