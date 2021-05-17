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
from numpy import in1d,nonzero,insert # OA 15/8/20

class modflowUsg:
    def __init__(self,core):
        self.core=core
    
    def buildMesh(self,opt):
        '''opt=old the .msh file is used if it is found, new : it is rewritten'''
        self.core.addin.mesh = self
        dct = self.core.diczone['Modflow'].dic
        mshType = self.core.getValueFromName('Modflow','MshType')
        if mshType >1: # case of true unstructured grid built through gmesh
            fmsh = self.core.fileName+'_out.msh'
            os.chdir(self.core.fileDir)
            if fmsh not in os.listdir(os.getcwd()): opt = 'new'
            dicD = dct['disu.3'] # where the domain boundaries are 
            idom = dicD['name'].index('domain')
            dcoo = dicD['coords'][idom]
            dcoo1 = self.verifyDomain(dcoo)
            npts = len([b for b in dicD['name'] if b[:5]=='point']) # OA 15/8/20
            liNames = [b for b in dicD['name'] if b[:4] in ['line','faul']]
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
                os.system(bindir+'gmsh gm_in.txt -2 -smooth 5 -algo front2d -o '+fmsh) #meshadapt
            f1 = open(fmsh,'r');s = f1.read();f1.close();print('gmesh file read')
            nbdy = len(dcoo1) # 15/8/20
            nodes,elements,self.dicFeats = readGmshOut(s,outline=True,nbdy=nbdy,liNames=liNames); # OA 15/8/20
            self.dicFeats['points']=list(range(nbdy,nbdy+npts)) # OA 15/8/20
            points,elts = nodes[:,1:3],elements[:,-3:]
            #points,elts,lines = self.verifyMesh(points,elts,lines)
            #npts = len(points)
            #ndom = where(array(lines,ndmin=2)[:,2]==0)[0][0] # this is the place where we get back to pt 0
            #line_dom = unique(lines[:ndom])
            #dns = float(dicD['value'][dicD['name'].index('domain')])
        if mshType == 0: # regular cells
            msh = myRect(self)
            msh.calc()
        if mshType == 1: # nested
            #msh = myRect(self,nodes,elements)
            #msh.calc()
            pass
        if mshType == 2: # triangle
            msh = myTrg(self,nodes,elements)
            msh.calc()
            xe,ye = self.elcenters[:,0],self.elcenters[:,1]
            self.trg = mptri.Triangulation(xe,ye) #
        if mshType == 3: # voronoi
            msh = myVor(self)
            msh.transformVor(points,elts,dcoo1,dicD,self.dicFeats) # OA 15/8/20
            xn,yn = self.nodes[:,0],self.nodes[:,1];print('voronoi made')
            self.trg = mptri.Triangulation(xn,yn) #,triangles=elts) 
        self.addMeshVects()
        Zblock = makeZblock(self.core)
        nlay = getNlayers(self.core)
        if self.core.addin.getDim()=='3D': 
            thick = Zblock[:-1]-Zblock[1:]
            self.core.addin.get3D();print('start 3D')
            self.add3d(nlay,thick);print('3D done')
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
        #xd,yd = zip(*dcoo);eps = (max(xd)-min(xd)+max(yd)-min(yd))/100
        for i in range(1,len(dcoo)-1):
            if dcoo[i]!=dcoo[i-1]: dcout.append(dcoo[i])
        return dcout
    
#    def getBoundaries(self,dcoo,nodes,lines):
#        poly=dcoo*1; poly.append(poly[0])
#        lcpoly = lcoefsFromPoly(poly)
#        line_dom = []
#        for i in range(len(lines)):
#            p0 = [(nodes[lines[i,0],1],nodes[lines[i,0],2]),(nodes[lines[i,1],1],nodes[lines[i,1],2])]
#            lc = lcoefsFromPoly(p0)
#            if any(sum(lc-lcpoly,axis=0)==0): line_dom.append(i)
#        return line_dom
        
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
        
    def connectZones(self):
        '''this function creates a local dictionnary that details teh connection between transp/pht3d
        zones and the flow ones. 1st from each point in a zone search for the corresp cell.
        2nd compare the cells with all of each modflow zones
        '''
class myRect:
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
                if i>0: a.append(j*nx+i-1);b.append(dy[j]);c.append(dx[i]/2);d.append(pi/2) # left
                if j>0: a.append((j-1)*nx+i);b.append(dx[i]);c.append(dy[j]/2);d.append(0) # top
                if i<nx-1: a.append(j*nx+i+1);b.append(dy[j]);c.append(dx[i]/2);d.append(pi/2) # right
                if j<ny-1: a.append((j+1)*nx+i);b.append(dx[i]);c.append(dy[j]/2);d.append(0) # bottom
                cn.append(a);fa.append(b);ds.append(c);agl.append(d)
                ex.append(array([xv[i],xv[i+1],xv[i+1],xv[i]]))
                ey.append(array([yv[j],yv[j],yv[j+1],yv[j+1]]))
                p.carea[j*nx+i] = dy[j]*dx[i]
        p.cneighb, p.fahl,p.cdist,p.angl,p.elx,p.ely = cn,fa,ds,agl,ex,ey   
        xc,yc = array([mean(a) for a in p.elx]),array([mean(a) for a in p.ely])
        p.elcenters = c_[xc,yc]
        p.nodes = p.elcenters
        
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
        
    def transformVor(self,points,elts,dcoo,dicD,dicFeats): # OA 15/8/20 dicFeats
        '''reads the voronoi to create several objects:
        cells ordered as the initial points, removing the ones that contain -1
        - nodes : the inital points (on the triangular gmesh grid)
        - carea: for each cell its area ((x1y2-x2y1)+(x2y1-x3y2)..)/2
        - (for each cell a list of n-1 ridges (the sides))
        - fahl : for each cell a list of ridges length
        - cneighb: for each cell a list of n-1 connecting other cell
        - cdist: for each cell a list of n-1 distances from node to face
        - elx and ely he coordinates of the elements nodes
        note : cdist, cneighb, fahl only include the links to cells that are inside 
        the domain an thus contain elements
        dcoo: coordinates of the domain
        '''
        #line_dom = self.getLineDomains(dcoo,lines)
        p = self.parent
        npts,a = shape(points);nelts,a = shape(elts)
        l0 = []
        p.carea,p.cneighb,p.cdist,p.elx,p.ely=l0*1,l0*1,l0*1,l0*1,l0*1
        p.nconnect,p.nodes = 0,[]
        line_dom = dicFeats['bc0'] # OA 15/8/20 + 2 lines bleow 
        for i in range(1,len(dcoo)): line_dom = r_[line_dom,dicFeats['bc'+str(i)]]# 15/8/20
        line_dom = unique(line_dom) 
        xdom,ydom = zip(*dcoo)
        dmax = (max(xdom)-min(xdom)+max(ydom)-min(ydom))/40;eps=dmax/1e9
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
        indx = nonzero(in1d(M[:,1],line_dom)&in1d(M[:,2],line_dom))[0] # OA 15/8/20
        M0 = M[indx,:] # where the segments are both in bdy
        M1 = M0*1
        M1[:,1],M1[:,2] = M0[:,2],M0[:,1] # invert the points
        M1[:,-1] = mod(M1[:,-1]+pi,2*pi) # add pi to the angle
        M2 = r_[M1,M]
        u,indx = unique(M2[:,1]*1e6+M2[:,2],return_index=True)# OA 15/4/21 1e4-> 1e6
        M = M2[indx,:]
        indx = argsort(M[:,1]*10+M[:,-1]) # sort again
        M = M[indx,:]
        # store the data as list where the list index is the cell nb
        M = self.putInDomain(elts,dcoo,M,dicFeats) # OA adede 25/4/20
        ip1 = M[:,1].astype('int')
        indsplit = where(ip1[1:]-ip1[:-1]==1)[0]+1 # this is the nb of points for each elt
        p.nconnect = shape(M)[0]
        p.cneighb = split(M[:,2].astype('int'),indsplit);nc=len(p.cneighb)
        lvx = split(M[:,3],indsplit);p.elx = [r_[v,v[0]] for v in lvx]
        lvy = split(M[:,4],indsplit);p.ely = [r_[v,v[0]] for v in lvy]
        p.cdist = split(M[:,5],indsplit)
        xmd1,ymd1 = split(M[:,6],indsplit),split(M[:,7],indsplit)
        # redo differently elx,ely and length for the poly at the boundaries (the rest is ok)
        indx = [0]*(max(line_dom)+1) # to remove the face that are not connected
        p.ncorn = len(dcoo)
        for ip1 in range(p.ncorn): # points on angles
            ip2 = p.cneighb[ip1]
            if (ip2[0] in line_dom)&(ip2[-1] in line_dom):# the 2 pts at the end of the poly are in line_dom
                p.elx[ip1][0] = xmd1[ip1][0]
                p.ely[ip1][0] = ymd1[ip1][0]
                p.elx[ip1][-1] = xmd1[ip1][-1]
                p.ely[ip1][-1] = ymd1[ip1][-1] 
                p.elx[ip1] = r_[p.elx[ip1],points[ip1,0],p.elx[ip1][0]]
                p.ely[ip1] = r_[p.ely[ip1],points[ip1,1],p.ely[ip1][0]]
                indx[ip1] = [len(ip2),len(ip2)+1] # OA 24/4/20
            else : # the 2 pts in line_dom are in the middle of the poly
                ia = nonzero(in1d(ip2,line_dom))[0]
                p.elx[ip1] = r_[p.elx[ip1][:ia[0]+1],xmd1[ip1][ia[0]],points[ip1,0],xmd1[ip1][ia[1]],p.elx[ip1][ia[-1]+1:]]
                p.ely[ip1] = r_[p.ely[ip1][:ia[0]+1],ymd1[ip1][ia[0]],points[ip1,1],ymd1[ip1][ia[1]],p.ely[ip1][ia[-1]+1:]]
                indx[ip1] = ia+1 # OA 24/4/20
        line_dom1 = line_dom[p.ncorn:]
        for ip1 in line_dom1: #range(line_dom[0],ndom+1): # points on lines
            ip2 = p.cneighb[ip1]
            if (ip2[0] in line_dom)&(ip2[-1] in line_dom):# the 2 pts at the end of the poly are in line_dom
                p.elx[ip1][0] = xmd1[ip1][0]
                p.ely[ip1][0] = ymd1[ip1][0]
                p.elx[ip1][-1] = xmd1[ip1][-1]
                p.ely[ip1][-1] = ymd1[ip1][-1] 
                p.elx[ip1] = r_[p.elx[ip1],p.elx[ip1][0]]
                p.ely[ip1] = r_[p.ely[ip1],p.ely[ip1][0]]
                indx[ip1] = -1 #len(ip2))
            else : # the 2 pts in line_dom are in the middle of the poly
                ia = nonzero(in1d(ip2,line_dom))[0]
                p.elx[ip1] = r_[p.elx[ip1][:ia[0]+1],xmd1[ip1][ia],p.elx[ip1][ia[-1]+1:]]
                p.ely[ip1] = r_[p.ely[ip1][:ia[0]+1],ymd1[ip1][ia],p.ely[ip1][ia[-1]+1:]]
                indx[ip1] = ia[-1]
        # make the faults
        p.nodes = array(points,ndmin=2)
        lfaults = [b for b in dicFeats if b[:5]=='fault']
        for fau in lfaults:
            lcell = unique(dicFeats[fau])
            iz = dicD['name'].index(fau) # index of fault zone
            lcell,ilines = self.sortFaultCells(p,lcell,len(dicD['coords'][iz])-1)
            lcoefs = lcoefsFromPoly(dicD['coords'][iz]);
            for ic in range(len(lcell)):
                p = self.splitCell(p,ic,lcell,lcoefs,ilines)
        nc= len(p.elx) # nb of cells has changed
        p.fahl = [sqrt((p.elx[i][1:]-p.elx[i][:-1])**2+(p.ely[i][1:]-p.ely[i][:-1])**2) for i in range(nc)]
        p.angl = []
        for i in range(nc) :
            angl = arctan((p.ely[i][1:]-p.ely[i][:-1])/(p.elx[i][1:]-p.elx[i][:-1]))
            angl[p.elx[i][1:] < p.elx[i][:-1]] += pi
            angl[angl<0]=2*pi+angl[angl<0]
            p.angl.append(angl)
        for ip1 in range(p.ncorn): # OA 26/4/20 changed to indx
            p.fahl[ip1] = delete(p.fahl[ip1],indx[ip1])
            p.angl[ip1] = delete(p.angl[ip1],indx[ip1])
        for ip1 in line_dom1: # OA 26/4/20 changed to indx
            p.fahl[ip1] = delete(p.fahl[ip1],indx[ip1])
            p.angl[ip1] = delete(p.angl[ip1],indx[ip1])
        la = [abs(sum(p.elx[i][:-1]*p.ely[i][1:]-p.elx[i][1:]*p.ely[i][:-1]))/2 for i in range(nc)]
        p.carea = array(la);p.indx = indx;
        p.elcenters = p.nodes*1 # OA 9/1/21
        for i in line_dom: p.elcenters[i] = (mean(p.elx[i]),mean(p.ely[i]))
        l = []
        for k in p.dicFeats.keys():
            if k[:2]=='bc': l.append(unique(p.dicFeats[k]))
        for i in range(p.ncorn): p.dicFeats['bc_cell'+str(i)]=l[i]

    def putInDomain(self,elts,dcoo,M,dicFeats):  # added 25/4/20 some poits can be outside domain
        poly=dcoo*1;poly.append(dcoo[0]);npl=len(poly)
        lcoefs=lcoefsFromPoly(poly)
        iIn=array(pointsInPoly(M[:,3],M[:,4],poly,lcoefs))
        lout = where(iIn==0)[0];nl=len(lout)
        lc1 = [];poly=array(poly,ndmin=2)
        # determine the line of interest
        lines = [unique(dicFeats[a]) for a in dicFeats.keys() if a[:2]=='bc']
        it0 = unique(M[lout,0]).astype('int') # index of the triangle on which is the outlier
        # on whcich bdy is the triangle
        for i in it0:
            for ip in range(npl-1):
                if sum(in1d(lines[ip],elts[i]))==2: 
                    lc1.append(lcoefs[:,ip])
        # set the point symetric to the line (ax+by+c=0) ici ax+by=1
        # y’=(-2abx-y(b^2-a^2)-2bc)/(b^2+a^2) et x’=x-a(y-y’)/b
        for i,it in enumerate(it0):
            a,b = lc1[i];c=-1
            x,y = M[M[:,0]==it,3],M[M[:,0]==it,4]
            yn = (-2*a*b*x-y*(b**2-a**2)-2*b*c)/(b**2+a**2+1e-6) #OA 15/4/21 added 1e-6
            M[M[:,0]==it,4] = yn
            M[M[:,0]==it,3] = x-a/(b+1e-6)*(y-yn) #OA 15/4/21 added 1e-6
#        for i in range(nl):
#            a,b = lc1[i,:];c=-1;i1 = lout[i]; x,y = ptx[i1],pty[i1]
#            pty1[i1] = (-2*a*b*x-y*(b**2-a**2)-2*b*c)/(b**2+a**2)
#            ptx1[i1] = x-a/b*(y-pty1[i1])
        return M
    
    def sortFaultCells(self,p,lcell,nlines): # OA added 17/8/20
        nc = len(lcell);lcold=list(lcell);cnew = lcold[0];lcnew = [cnew]
        ilines = [];count=1
        for i in range(nc-1):
            cnew = lcold[nonzero(in1d(lcold,p.cneighb[cnew]))[0][0]]
            lcold.remove(lcnew[-1])
            lcnew.append(cnew)
            ilines.append(count-1)
            if cnew==lcell[count]: count +=1
        ilines.append(ilines[-1])
        return lcnew,ilines
        
    def splitCell(self,p,ic,lcell,lcoefs,ilines) : # OA added 16/8/20...
        '''
        split the cells in the middle where the fault line passes (not at the position of the fault)
        the cell that has a -1 sign is kept at cll position, the sgn=1 is added
        at the end of p
        '''
        nctot = len(p.elx);cll=lcell[ic];il0=ilines[max(0,ic-1)];il1=ilines[ic]
        sgn0 = sign(p.elx[cll]*lcoefs[0,il0]+p.ely[cll]*lcoefs[1,il0]-1) # this sign gives nfo about the side of the fualt ine
        sgn = sign(p.elx[cll]*lcoefs[0,il1]+p.ely[cll]*lcoefs[1,il1]-1) # this sign gives nfo about the side of the fualt ine
        imod0 = where(abs(sgn0[1:]-sgn0[:-1])!=0)[0] # where change ccur
        imod = where(abs(sgn[1:]-sgn[:-1])!=0)[0] # where change ccur
        if il1 != il0: imod[1] = imod0[1]
        xn = (p.elx[cll][imod]+p.elx[cll][imod+1])/2
        yn = (p.ely[cll][imod]+p.ely[cll][imod+1])/2
        elx1 = r_[p.elx[cll][:imod[0]+1],xn,p.elx[cll][imod[1]+1:]]
        ely1 = r_[p.ely[cll][:imod[0]+1],yn,p.ely[cll][imod[1]+1:]]
        elx2 = r_[xn[-1::-1],p.elx[cll][imod[0]+1:imod[1]+1],xn[1]]
        ely2 = r_[yn[-1::-1],p.ely[cll][imod[0]+1:imod[1]+1],yn[1]]
        cdm = mean(p.cdist[cll])
        if sgn[0]==-1: nb1,nb2=nctot,cll
        else : nb1,nb2=cll,nctot
        cn1 = r_[p.cneighb[cll][:imod[0]+1],nb1,p.cneighb[cll][imod[1]:]]
        cd1 = r_[p.cdist[cll][:imod[0]+1],cdm,p.cdist[cll][imod[1]:]]
        cn2 = r_[nb2,p.cneighb[cll][imod[0]:imod[1]+1]]
        cd2 = r_[cdm,p.cdist[cll][imod[0]:imod[1]+1]]
        # change the elx,ely and neighbors for the neighb of the current cell
        cngh = p.cneighb[cll]
        if (ic>0)&(ic<len(lcell)-1):
            cn2[where(cn2==lcell[ic-1])[0]]=nctot-1 # the previous valeu was not corected
        else : # first and last cell
            i0 = 1- nonzero(in1d(cngh[imod],lcell))[0][0] # indx of the neighb not in lcell
            ic0 = cngh[imod[i0]];ic1 = list(p.cneighb[ic0]).index(cll)
            p.cneighb[ic0] = insert(p.cneighb[ic0],ic1+1,nctot)
            p.cdist[ic0] = insert(p.cdist[ic0],ic1+1,p.cdist[ic0][ic1])
            p.elx[ic0] = insert(p.elx[ic0],ic1+1,xn[i0])
            p.ely[ic0] = insert(p.ely[ic0],ic1+1,yn[i0])
        # save to p
        if sgn[0]==-1:
            p.elx[cll]=elx1;p.ely[cll]=ely1;p.elx.append(elx2);p.ely.append(ely2)
            p.cneighb[cll]=cn1;p.cneighb.append(cn2)
            p.cdist[cll]=cd1;p.cdist.append(cd2)
            for c in cn2 : 
                if c not in lcell : p.cneighb[c][p.cneighb[c]==cll]=nctot
        else :
            p.elx[cll]=elx2;p.ely[cll]=ely2;p.elx.append(elx1);p.ely.append(ely1)
            p.cneighb[cll]=cn2;p.cneighb.append(cn1)
            p.cdist[cll]=cd2;p.cdist.append(cd1)
            for c in cn1 : 
                if c not in lcell : p.cneighb[c][p.cneighb[c]==cll]=nctot
        # area will be calc later 
        p.nodes = r_[p.nodes,p.nodes[cll:cll+1]]
        return p
       
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