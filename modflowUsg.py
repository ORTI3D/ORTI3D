# -*- coding: utf-8 -*-
"""
Created on Sun Sep 17 07:08:20 2017

@author: oatteia
"""
from .config import *
from .geometry import *
from .geometryMesh import *

class modflowUsg(unstructured):
    def __init__(self,core):
        self.core=core
        
    def buildMesh(self,opt):
        self.buildMesh0('Modflow',opt)
        mshType = self.core.getValueFromName('Modflow','MshType')
        if mshType == 0: # regular cells
            msh = usgRect(self)
            msh.calc()
        if mshType == 1: # nested
            #msh = myRect(self,nodes,elements)
            #msh.calc()
            pass
        if mshType == 2: # triangle
            msh = usgTrg(self,nodes,elements)
            msh.calc()
            xe,ye = self.elcenters[:,0],self.elcenters[:,1]
            self.trg = mptri.Triangulation(xe,ye) #
        if mshType == 3: # voronoi
            msh = myVor(self)
            msh.transformVor(self.points,self.elts,self.dcoo1,self.dicD,self.dicFeats) # OA 15/8/20
            xn,yn = self.nodes[:,0],self.nodes[:,1];print('voronoi made')
            self.trg = mptri.Triangulation(xn,yn) #,triangles=elts) 
        self.addMeshVects()
        Zblock = makeZblock(self.core)
        self.nlay = getNlayers(self.core)
        thick=[0]*self.nlay
        if self.core.addin.getDim()=='3D': 
            for i in range(self.nlay): thick[i] = ravel(Zblock[i]-Zblock[i+1])
            thick = array(thick)
            self.core.addin.get3D();print('start 3D')
            self.add3d(self.nlay,thick);print('3D done')
        else :
            thick = Zblock[0]-Zblock[-1]
            self.fahl = [array(lg)*thick[i] for i,lg in enumerate(self.fahl)]
                            
    def writeDisu(self):
        #6. Area(NDSLAY) - U1DREL, red only once for vertically consistent grid
        s = 'INTERNAL 1.0 (FREE)  0 #AREA \n'
        s += self.lst2lines(self.carea,'float')
        #7. IAC(NODES) - U1DINT
        s += 'INTERNAL 1 (FREE)   0 #IAC \n' # nb of connections
        s += self.lst2lines([len(x)+1 for x in self.cneighb],'int')
        # 8. JA(NJAG) - U1DINT # cell nb starts from 1 in Modflow!!
        s += 'INTERNAL 1 (FREE)   0 #JA \n' # list of neighbors for each cell
        s += '\n'.join([str(i+1)+' '+' '.join([str(x+1) for x in self.cneighb[i]]) for i in range(self.ncell)])
        s += '\n'
        #9. IVC(NJAG) - U1DINT # to be red only for vertical subdivision =0
        #10a. CL1(NJAGS) - U1DREL – for symertric input
        #10b. CL2(NJAGS) - U1DREL – for symertric input
        #11. CL12(NJAG) - U1DREL - for non symertric input
        s += 'INTERNAL 1.0 (FREE)  0  #CL12 dist node to face\n'
        s += '\n'.join(['0 '+' '.join(['%9.4e '%x for x in self.cdist[i]]) for i in range(self.ncell)])
        s += '\n'
        #12. FAHL(NJAG/NJAGS) - U1DREL -
        s += 'INTERNAL 1.0 (FREE)  -1  #FAHL face area \n'
        s += '\n'.join(['0 '+' '.join(['%9.4e '%x for x in self.fahl[i]]) for i in range(self.ncell)])
        s += '\n'
        #13. PEfahl NSTP TSMULT Ss/Tr
        #s +='1.0 1 1.0  SS \n' written in the modflow writer
        return s.replace('\n\n','\n')
        
    def lst2lines(self,lst,fmt='float'):
        s,l0 = '',len(lst)
        if fmt=='float': fmt1 = '%9.4e '
        if fmt=='int': fmt1 = '%5i '
        if l0<50:
            s += ' '.join([fmt1%a for a in lst])+ '\n'
        else :
            for i in range(int(l0/50)):
                s += ' '.join([fmt1%a for a in lst[i*50:(i+1)*50]])+'\n'
            s += ' '.join([fmt1%a for a in lst[(i+1)*50:]])+ '\n'
        return s
        
    def add3d(self,nlay,thk):
        '''transforms the 2D mesh in a 3D ones. nlay-1 cells are added
        area are not changed, but the connections are modified to represent the 3D,
        the length and distances are the same for each layer, except the one added for 3d
        carea1 is calculated for consistency in searching area'''
        ncell2d = self.ncell*1
        cneighb0,cdist0,carea0,fahl0,angl0 = self.cneighb,self.cdist,self.carea,self.fahl,self.angl
        cneighb1, cdist1, carea1, fahl1, angl1 = [],[],[],[],[]
        for il in range(nlay):
            print(il)
            if il==0: # add the cell below
                cneighb1.extend([r_[array(cn),i+ncell2d] for i,cn in enumerate(cneighb0)])
                cdist1.extend([r_[array(cd), thk[0,i]/2] for i,cd in enumerate(cdist0)])
                fahl1.extend([r_[array(fl)*thk[0,i],carea0[i]] for i,fl in enumerate(fahl0)]) 
                angl1.extend([r_[array(ag),0] for i,ag in enumerate(angl0)]) 
            elif il==nlay-1: # add the cell on top
                cneighb1.extend([r_[array(cn)+il*ncell2d,i+(il-1)*ncell2d] for i,cn in enumerate(cneighb0)])
                cdist1.extend([r_[array(cd), thk[-1,i]/2] for i,cd in enumerate(cdist0)])
                fahl1.extend([r_[array(fl)*thk[-1,i],carea0[i]] for i,fl in enumerate(fahl0)]) 
                angl1.extend([r_[array(ag),0] for i,ag in enumerate(angl0)]) 
            else : # add both cell on top and blow
                cneighb1.extend([r_[array(cn)+il*ncell2d,i+(il-1)*ncell2d,i+(il+1)*ncell2d] for i,cn in enumerate(cneighb0)])
                cdist1.extend([r_[array(cd), thk[il-1,i]/2, thk[il+1,i]/2] for i,cd in enumerate(cdist0)])
                fahl1.extend([r_[array(fl)*thk[il,i],carea0[i],carea0[i]] for i,fl in enumerate(fahl0)]) 
                angl1.extend([r_[array(ag),0,0] for i,ag in enumerate(angl0)]) 
            carea1.extend(carea0)
        self.cneighb,self.cdist,self.fahl,self.carea1,self.angl = cneighb1,cdist1,fahl1,carea1,angl1
        self.ncell =  ncell2d*nlay
        self.ncell_lay =  ncell2d
        self.nconnect = self.nconnect*nlay+ncell2d*2*(nlay-1)

class usgRect:
    def __init__(self,parent):
        ''' the grid is made of rectangular cells, include variables width or height'''
        self.parent = parent

    def calc(self):
        p = self.parent
        p.carea,p.fahl,p.cneighb,p.cdist=[],[],[],[]
        grd = self.parent.core.addin.getFullGrid()
        dx,dy = array(grd['dx']), array(grd['dy']);
        nx,ny,x0,y0 = grd['nx'],grd['ny'],grd['x0'],grd['y0']
        xv,yv = r_[x0,x0+cumsum(dx)],r_[y0,y0+cumsum(dy)]
        ncell = nx*ny;p.ncell = ncell
        p.carea = zeros(ncell)
        p.nconnect = (nx-2)*(ny-2)*4+(nx-2)*3*2+(ny-2)*3*2+4*2
        cn,fa,ds,agl,ex,ey=[],[],[],[],[],[]
        for j in range(ny):
            for i in range(nx):
                a,b,c,d = [],[],[],[]
                if i>0: a.append(j*nx+i-1);b.append(dy[j]);c.append(dx[i]/2);d.append(pi) # left
                if j>0: a.append((j-1)*nx+i);b.append(dx[i]);c.append(dy[j]/2);d.append(-pi/2) # bottom
                if i<nx-1: a.append(j*nx+i+1);b.append(dy[j]);c.append(dx[i]/2);d.append(0) # right
                if j<ny-1: a.append((j+1)*nx+i);b.append(dx[i]);c.append(dy[j]/2);d.append(pi/2) # top
                cn.append(a);fa.append(b);ds.append(c);agl.append(d)
                ex.append(array([xv[i],xv[i+1],xv[i+1],xv[i]]))
                ey.append(array([yv[j],yv[j],yv[j+1],yv[j+1]]))
                p.carea[j*nx+i] = dy[j]*dx[i]
        p.cneighb, p.fahl,p.cdist,p.angl,p.elx,p.ely = cn,fa,ds,agl,ex,ey   
        xc,yc = array([mean(a) for a in p.elx]),array([mean(a) for a in p.ely])
        p.elcenters = c_[xc,yc]
        p.nodes = p.elcenters
        
class usgTrg:
    def __init__(self,parent, nodes, elements):
        self.parent,self.nodes,self.elements = parent,nodes, elements

    def calc(self):
        p = self.parent
        p.carea,p.fahl,p.cneighb,p.lcoefs,p.cdist,p.angles=[],[],[],[],[],[]
        p.nodes,p.elements = self.nodes,self.elements
        elt = p.elements[:,3:]
        ncell,a =shape(elt)
        p.ncell,p.ncell_lay,p.nconnect = ncell,ncell,0
        p.nd_elt = [[] for i in range(len(self.nodes))]
        p.nd_dst = [[] for i in range(len(self.nodes))]
        for ic in range(ncell):
            inod = elt[ic,:];#print inod
            inod1 = r_[inod,inod[0]]
            x,y = p.nodes[inod1,1],p.nodes[inod1,2]
            xp,yp = mean(x[:-1]),mean(y[:-1])
            for j in inod: 
                p.nd_elt[j].append(ic)
                p.nd_dst[j].append(sqrt((p.nodes[j,1]-xp)**2+(p.nodes[j,2]-yp)**2))
            p.carea.append(abs(sum(x[:-1]*y[1:]-x[1:]*y[:-1]))/2)
            neighb,llcoef,dlist,llist,alist = [],[],[],[],[]
            for i in range(3):
                #find cells neihgboring the current one
                il1,il2 = list(where(elt==inod1[i])[0]),list(where(elt==inod1[i+1])[0]);
                il3 = list(set(il1)&set(il2))
                il3.remove(ic)
                if len(il3)>0:
                    neighb.extend(il3)
                    # finds the distance from pt to face, using the eq of le triangle sides
                    c = 1
                    xym = array([[x[i],y[i]],[x[i+1],y[i+1]]])
                    if sum(abs(xym[:,0]))==0: 
                        lcoef =array([1,0]);c=0
                    elif sum(abs(xym[:,1]))==0: 
                        lcoef =array([0,1]);c=0
                    else : 
                        lcoef = solve(xym,ones((2,1))[:,0])
                    llcoef.append(lcoef)
                    dst = abs(dot(array([xp,yp]),lcoef)-c)/sqrt(sum(lcoef**2))
                    dlist.append(dst)
                    # find the length of the vertex and the angle
                    llist.append(sqrt((x[i]-x[i+1])**2+(y[i]-y[i+1])**2))
                    angl = arctan((y[i+1]-y[i])/(x[i+1]-x[i]))-pi/2
                    if x[i+1]<x[i] : angl += pi
                    alist.append(angl)
            p.cneighb.append(neighb)
            p.lcoefs.append(llcoef)
            p.cdist.append(dlist)
            p.fahl.append(llist)
            p.angles.append(alist)
            p.nconnect += len(dlist)
        p.carea1 = p.carea*1
        p.elx = p.nodes[p.elements[:,3:],1]
        p.ely = p.nodes[p.elements[:,3:],2]
        xc,yc = mean(p.elx,axis=1),mean(p.ely,axis=1)
        p.elcenters = c_[xc,yc]
