# -*- coding: utf-8 -*-
"""
Created on Sun Sep 17 07:08:20 2017

@author: oatteia
"""
#from scipy.spatial import Voronoi
#from shapely.geometry import Polygon
from scipy import pi
from operator import itemgetter
from .config import *
from .geometry import *
from .myInterpol import *
import matplotlib.tri as mptri

class modflowUsg:
    def __init__(self,core):
        self.core=core
    
    def buildMesh(self,opt):
        '''opt=old the .msh file is used if it is found, new : it is rewritten'''
        self.core.addin.mesh = self
        dct = self.core.diczone['Modflow'].dic
        fmsh = self.core.fileName+'_out.msh'
        os.chdir(self.core.fileDir)
        if fmsh not in os.listdir(os.getcwd()): opt = 'new'
        dicD = dct['disu.3'] # where the domain boundaries are 
        idom = dicD['name'].index('domain')
        dcoo = dicD['coords'][idom]
        dcoo = self.verifyDomain(dcoo)
        # if 'lpf.8' in dct: dicK = dct['lpf.8'] # K hydraul
        # else : 
        dicK = {'name':[]}
        if opt=='new':
            s = gmeshString(self.core,dicD,dicK)
            f1 = open('gm_in.txt','w');f1.write(s);f1.close()
            bindir = self.core.baseDir+os.sep+'bin'+os.sep
            #clmax = self.core.dicval['Modflow']['disu.3'][0];print('start gmesh')
            #os.system(bindir+'gmsh gm_in.txt -2 -smooth 5 -clmax '+str(clmax)+' -o '+fmsh)
            #os.system(bindir+'gmsh gm_in.txt -2 -smooth 5 -o '+fmsh)
            os.system(bindir+'gmsh gm_in.txt -2 -smooth 5 -algo meshadapt -o '+fmsh)
        f1 = open(fmsh,'r');s = f1.read();f1.close();print('gmesh file read')
        nodes,elements,line_dom = readGmshOut(s,outline=True);
        points,elts = nodes[:,1:3],elements[:,-3:]
        points,elts,line_dom = self.verifyMesh(points,elts,line_dom)
        npts = len(points)
        #bbox = [amin(points[:,0]),amax(points[:,0]),amin(points[:,1]),amax(points[:,1])]
        dns = float(dicD['value'][dicD['name'].index('domain')])
        mshType = self.core.getValueFromName('Modflow','MshType')
        if mshType == 0: # triangle
            msh = myTrg(self,nodes,elements)
            msh.calc()
            xe,ye = self.elcenters[:,0],self.elcenters[:,1]
            self.trg = mptri.Triangulation(xe,ye) #
        if mshType == 1: # voronoi
            msh = myVor(self)
            msh.transformVor(points,elts,dcoo,line_dom)
            xn,yn = self.nodes[:,0],self.nodes[:,1]
            self.trg = mptri.Triangulation(xn,yn) #,triangles=elts) 
        self.addMeshVects()
        Zblock = makeZblock(self.core)
        nlay = getNlayers(self.core)
        if self.core.addin.getDim()=='3D': 
            thick = Zblock[:-1]-Zblock[1:]
            self.core.addin.get3D()
            self.add3d(nlay,thick) # !!! media = layer up to now
        else :
            thick = Zblock[0]-Zblock[-1]
            self.fahl = [array(lg)*thick[i] for i,lg in enumerate(self.fahl)]

    def getCenters(self): return self.nodes[:,0],self.nodes[:,1]
    def getNumber(self,typ): 
        if typ == 'nodes': return len(self.nodes)
        if typ == 'elements': return self.ncell_lay
        if typ == 'tot_elements': return self.ncell
        
    def verifyDomain(self,dcoo):
        dcout = [dcoo[0]]
        for i in range(1,len(dcoo)):
            if dcoo[i]!=dcoo[i-1]: dcout.append(dcoo[i])
        return dcout
        
    def verifyMesh(self,pts,elts,line_dom):
        '''try to eliminate double points'''
        ndom = len(line_dom)
        dpts = sum(abs(pts[1:,:]-pts[:-1,:]),axis=1);dpts=r_[dpts,1]
        ptout = pts[dpts!=0]
        ldout = line_dom[dpts[:ndom]!=0]
        indx = where(dpts==0)[0]
        for i in list(indx):
            i1 = where(elts==i)[0]
            i2 = where(elts[i1]==i-1)[0]
            elts = delete(elts,i1[i2],axis=0)
            elts[elts>=i] -= 1
        return ptout,elts,ldout
        
    def addMeshVects(self):
        self.ncell,self.ncell_lay = len(self.elx),len(self.elx)
        ic,self.idc=0,[]
        for a in self.elx: self.idc.append([ic,ic+len(a)]);ic+=len(a)
        self.elxa,self.elya=[],[]
        for a in self.elx: self.elxa.extend(a)
        for a in self.ely: self.elya.extend(a)
                        
    def writeDisu(self):
        mh = self
        #6. Area(NDSLAY) - U1DREL, red only once for vertically consistent grid
        s = 'INTERNAL 1.0 (FREE)  0 #AREA \n'
        s += self.lst2lines(mh.carea,'float')
        #7. IAC(NODES) - U1DINT
        s += 'INTERNAL 1 (FREE)   0 #IAC \n' # nb of connections
        s += self.lst2lines([len(x)+1 for x in mh.cneighb],'int')
        # 8. JA(NJAG) - U1DINT # cell nb starts from 1 in Modflow!!
        s += 'INTERNAL 1 (FREE)   0 #JA \n' # list of neighbors for each cell
        s += '\n'.join([str(i+1)+' '+' '.join([str(x+1) for x in mh.cneighb[i]]) for i in range(mh.ncell)])
        s += '\n'
        #9. IVC(NJAG) - U1DINT # to be red only for vertical subdivision =0
        #10a. CL1(NJAGS) - U1DREL – for symertric input
        #10b. CL2(NJAGS) - U1DREL – for symertric input
        #11. CL12(NJAG) - U1DREL - for non symertric input
        s += 'INTERNAL 1.0 (FREE)  0  #CL12 dist node to face\n'
        s += '\n'.join(['0 '+' '.join(['%9.4e '%x for x in mh.cdist[i]]) for i in range(mh.ncell)])
        s += '\n'
        #12. FAHL(NJAG/NJAGS) - U1DREL -
        s += 'INTERNAL 1.0 (FREE)  -1  #FAHL face area \n'
        s += '\n'.join(['0 '+' '.join(['%9.4e '%x for x in mh.fahl[i]]) for i in range(mh.ncell)])
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
        # 1st layer
        cneighb0, cdist0, carea0, fahl0 = self.cneighb, self.cdist, self.carea,self.fahl
        cneighb1, cdist1, carea1, fahl1 = [],[],[],[]
        for il in range(nlay):
            if il==0: # add the cell below
                cneighb1.extend([r_[array(cn),i+ncell2d] for i,cn in enumerate(cneighb0)])
                cdist1.extend([r_[array(cd), thk[0,i]/2] for i,cd in enumerate(cdist0)])
                fahl1.extend([r_[array(fl)*thk[0,i],carea0[i]] for i,fl in enumerate(fahl0)]) 
            elif il==nlay-1: # add the cell on top
                cneighb1.extend([r_[array(cn)+il*ncell2d,i+(il-1)*ncell2d] for i,cn in enumerate(cneighb0)])
                cdist1.extend([r_[array(cd), thk[-1,i]/2] for i,cd in enumerate(cdist0)])
                fahl1.extend([r_[array(fl)*thk[-1,i],carea0[i]] for i,fl in enumerate(fahl0)]) 
            else : # add both cell on top and blow
                cneighb1.extend([r_[array(cn)+il*ncell2d,i+(il-1)*ncell2d,i+(il+1)*ncell2d] for i,cn in enumerate(cneighb0)])
                cdist1.extend([r_[array(cd), thk[il-1,i]/2, thk[il+1,i]/2] for i,cd in enumerate(cdist0)])
                fahl1.extend([r_[array(fl)*thk[il,i],carea0[i],carea0[i]] for i,fl in enumerate(fahl0)]) 
            carea1.extend(carea0)
        self.cneighb, self.cdist, self.fahl,self.carea1 = cneighb1, cdist1, fahl1,carea1
        self.ncell =  ncell2d*nlay
        self.ncell_lay =  ncell2d
        self.nconnect = self.nconnect*nlay+ncell2d*2*(nlay-1)

class myTrg:
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
                    angl = arctan((y[i+1]-y[i])/(x[i+1]-x[i]))
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
                        
class myVor:
    def __init__(self,parent):
        self.parent = parent
        
    def transformVor(self,points,elts,dcoo,line_dom):
        '''reads the voronoi to create several objects:
        cells ordered as the initial points, removing the ones that contain -1
        - nodes : the inital points
        - xyverts : an array of vertices coordinates
        - cverts: for each cell a list of n vertices (points of the cell)
        - rarray: same but in a form of an array with last point till the max size
        - carea: for each cell its area ((x1y2-x2y1)+(x2y1-x3y2)..)/2
        - (for each cell a list of n-1 ridges (the sides))
        - fahl : for each cell a list of ridges length
        - cneighb: for each cell a list of n-1 connecting other cell
        - cdist: for each cell a list of n-1 distances from node to face
        note : cdist, cneighb, fahl only include the links to cells that are inside 
        the domain an thus contain elements
        dcoo: coordinates of the domain
        '''
        p = self.parent
        npts,a = shape(points);nelts,a = shape(elts)
        l0 = []
        p.carea,p.fahl,p.cneighb,p.cdist,p.elx,p.ely=l0*1,l0*1,l0*1,l0*1,l0*1,l0*1
        p.nconnect,p.nodes = 0,[]
        line_dom = list(line_dom);nd=len(line_dom)
        xdom,ydom = zip(*dcoo)
        dmax = (max(xdom)-min(xdom)+max(ydom)-min(ydom))/40;eps=dmax/1e6
        # new method with arrays
        it = array(arange(nelts),ndmin=2).T
        xc,yc = self.calcCentre(points,elts,eps)
        xp,yp = points[elts,0],points[elts,1] # triangle points in matrices nelts,3
        xmid = c_[(xp[:,0:2]+xp[:,1:3])/2,(xp[:,0]+xp[:,2])/2] # mid point in the triangle face
        ymid = c_[(yp[:,0:2]+yp[:,1:3])/2,(yp[:,0]+yp[:,2])/2]
        dx = xp[:,1:3]-xp[:,0:2]+eps; dx = c_[dx,xp[:,0]-xp[:,2]+eps]
        dy = yp[:,1:3]-yp[:,0:2]; dy = c_[dy,yp[:,0]-yp[:,2]]
        theta = arctan(dy/dx)  # angle from point 1 to 2
        theta[dx<0] += pi
        # to calculate the mid point and length
        xv,yv = c_[xp[:,0],xmid[:,0],xc,xmid[:,2],xp[:,0]],c_[yp[:,0],ymid[:,0],yc,ymid[:,2],yp[:,0]]
        l1 = sqrt((xv[:,1:4]-xv[:,:3])**2+(yv[:,1:4]-yv[:,:3])**2);
        #pd1,l1a,l1b = l1[:,0],l1[:,1],l1[:,2]
        xv,yv = c_[xp[:,1],xmid[:,1],xc,xmid[:,0],xp[:,1]],c_[yp[:,1],ymid[:,1],yc,ymid[:,0],yp[:,1]]
        l2 = sqrt((xv[:,1:4]-xv[:,:3])**2+(yv[:,1:4]-yv[:,:3])**2);
        #pd2,l2a,l2b = l2[:,0],l2[:,1],l2[:,2]
        xv,yv = c_[xp[:,2],xmid[:,2],xc,xmid[:,1],xp[:,2]],c_[yp[:,2],ymid[:,2],yc,ymid[:,0],yp[:,2]]
        l3 = sqrt((xv[:,1:4]-xv[:,:3])**2+(yv[:,1:4]-yv[:,:3])**2);
        #pd3,l3a,l3b = l3[:,0],l3[:,1],l3[:,2]
        # creating the three matrices of segments (p1,p2), (p2,p3),(p3,p1)
        M1 = c_[it,elts[:,0:2],xc,yc,l1[:,0],xmid[:,0],ymid[:,0],theta[:,0]]
        M2 = c_[it,elts[:,1:3],xc,yc,l2[:,0],xmid[:,1],ymid[:,1],theta[:,1]]
        M3 = c_[it,elts[:,2],elts[:,0],xc,yc,l3[:,0],xmid[:,2],ymid[:,2],theta[:,2]]
        M4 = r_[M1,M2,M3] # struct it,pt1,pt2,xc,yc,mid btw p1&p2, area
        # sort the segments by pt1 and then angle
        indx = argsort(M4[:,1]*10+M4[:,-1])
        M = M4[indx,:] # order by pt1 and angle
        # add the 2nd segments when they are at the boundary (there is only one in M now)
        ndom = max(line_dom) # the boundary segments are the first one
        M0 = M[(M[:,1]<=ndom)&(M[:,2]<=ndom),:] # where the segments are both in bdy
        M1 = M0*1
        M1[:,1],M1[:,2] = M0[:,2],M0[:,1] # invert the points
        M1[:,-1] = mod(M1[:,-1]+pi,2*pi) # add pi to the angle
        M2 = r_[M1,M]
        u,indx = unique(M2[:,1]*1e4+M2[:,2],return_index=True)# some segmts were already there
        M = M2[indx,:]
        indx = argsort(M[:,1]*10+M[:,-1]) # sort again
        M = M[indx,:]
        # store the data as list where the list index is the cell nb
        ip1 = M[:,1].astype('int')
        indsplit = where(ip1[1:]-ip1[:-1]==1)[0]+1 # this is the nb of points for each elt
        p.nconnect,a = shape(M)
        p.cneighb = split(M[:,2].astype('int'),indsplit);nc=len(p.cneighb)
        lvx = split(M[:,3],indsplit);p.elx = [r_[v,v[0]] for v in lvx]
        lvy = split(M[:,4],indsplit);p.ely = [r_[v,v[0]] for v in lvy]
        p.cdist = split(M[:,5],indsplit)
        xmd1,ymd1 = split(M[:,6],indsplit),split(M[:,7],indsplit)
        # redo differently elx,ely and length for the poly at the boundaries (the rest is ok)
        indx = []
        for ip1 in range(ndom):
            ip2 = p.cneighb[ip1]
            if (ip2[0]<=ndom)&(ip2[-1]<=ndom):# the 2 pts at the end of the poly are in line_dom
                p.elx[ip1][0] = xmd1[ip1][0]
                p.ely[ip1][0] = ymd1[ip1][0]
                p.elx[ip1][-1] = xmd1[ip1][-1]
                p.ely[ip1][-1] = ymd1[ip1][-1]
                p.elx[ip1] = r_[p.elx[ip1],p.elx[ip1][0]]
                p.ely[ip1] = r_[p.ely[ip1],p.ely[ip1][0]]
                indx.append(-1)
            else : # the 2 pts in line_dom are in the middle of the poly
                ia=where(ip2<=ndom)[0]
                p.elx[ip1] = r_[p.elx[ip1][:ia[0]+1],xmd1[ip1][ia],p.elx[ip1][ia[-1]+1:]]
                p.ely[ip1] = r_[p.ely[ip1][:ia[0]+1],ymd1[ip1][ia],p.ely[ip1][ia[-1]+1:]]
                indx.append(ia[-1])
        p.fahl = [sqrt((p.elx[i][1:]-p.elx[i][:-1])**2+(p.ely[i][1:]-p.ely[i][:-1])**2) for i in range(nc)]
        for ip1 in range(ndom):
            p.fahl[ip1] = delete(p.fahl[ip1],indx[ip1])
        la = [abs(sum(p.elx[i][:-1]*p.ely[i][1:]-p.elx[i][1:]*p.ely[i][:-1]))/2 for i in range(nc)]
        p.carea = array(la)
        p.nodes = array(points,ndmin=2)
        p.elcenters = p.nodes
        
    def verifyDomain(self,line_dom,points):
        lout=[line_dom[0]]
        for i in range(1,len(line_dom)):
            if sum(abs(points[line_dom[i]]-points[line_dom[i-1]]))!=0:
                lout.append(line_dom[i])
        return lout
        
    def calcCentre(self,points,elts,eps):
        ''' at https://cral.univ-lyon1.fr/labo/fc/Ateliers_archives/ateliers_2005-06/cercle_3pts.pdf
        (one sign was wrong)'''
        i1,i2,i3 = elts[:,0],elts[:,1],elts[:,2]
        x1,x2,x3 = points[i1,0],points[i2,0],points[i3,0]
        y1,y2,y3 = points[i1,1],points[i2,1],points[i3,1]
        # line equations
        dx1,dy1 = x2-x1, y2-y1+eps;a1 = -dx1/dy1
        dxy1 = x2**2-x1**2+y2**2-y1**2; b1 = dxy1/2/dy1
        dx2,dy2 = x3-x2, y3-y2+eps;a2 = -dx2/dy2
        dxy2 = x3**2-x2**2+y3**2-y2**2;b2 = dxy2/2/dy2
        # centers
        xc = (b2-b1)/(a1-a2);yc = a1*xc+b1
        return xc,yc