# -*- coding: utf-8 -*-
"""
Created on Sun Sep 17 07:08:20 2017

@author: oatteia
"""
#from scipy.spatial import Voronoi
from scipy import pi
from operator import itemgetter
from config import *
from geometry import *
from myInterpol import *
import matplotlib.tri as mptri

class modflowUNS:
    def __init__(self,core):
        self.core=core
    
    def buildMesh(self,opt):
        '''opt=old the .msh file is used if it is found, new : it is rewritten'''
        self.core.addin.mesh = self
        dct = self.core.diczone['Modflow'].dic
        fmsh = self.core.fileDir+os.sep+self.core.fileName+'_out.msh'
        os.chdir(self.core.fileDir)
        if fmsh not in os.listdir(os.getcwd()): opt = 'new'
        dicD = dct['disu.3'] # where the domain boundaries are + line .. + points (only pts and line will be used)
        idom = dicD['name'].index('domain')
        dcoo = dicD['coords'][idom]
        if dcoo[-1]==dcoo[0]: dcoo=dcoo[:-1]
        dicD['coords'][idom] = dcoo
        if dct.has_key('lpf.8'): dicK = dct['lpf.8'] # K hydraul
        else : dicK = {'name':[]}
        if opt=='new':
            s = gmeshString(self.core,dicD,dicK)
            f1 = open('gm_in.txt','w');f1.write(s);f1.close()
            bindir = self.core.baseDir+os.sep+'bin'+os.sep
            clmax = self.core.dicval['Modflow']['disu.3'][0]
            #os.system(bindir+'gmsh gm_in.txt -2 -smooth 5 -o '+fmsh)
            #os.system(bindir+'gmsh gm_in.txt -2 -clmax '+str(clmax)+' -o '+fmsh)
            os.system(bindir+'gmsh gm_in.txt -2 -smooth 5 -clmax '+str(clmax)+' -o '+fmsh)
        f1 = open(fmsh,'r');s = f1.read();f1.close()
        nodes,elements = readGmshOut(s);#print 'mfu 32',shape(nodes),shape(elements)
        points = nodes[:,1:3]
        npts0 = len(points)
        bbox = [amin(points[:,0]),amax(points[:,0]),amin(points[:,1]),amax(points[:,1])]
        dns = float(dicD['value'][dicD['name'].index('domain')])
        points = self.addSidesPts(points,bbox,dns);
        mshType = self.core.getValueFromName('Modflow','MshType')
        if mshType == 0: # triangle
            msh = myTrg(self,nodes,elements)
            msh.calc()
            #xn,yn = msh.nodes[:,1],msh.nodes[:,2]
            #self.trg = mptri.Triangulation(xn,yn,triangles=elements[:,-3:]) # be careful the numbering can be differnet from gmesh
            xe,ye = self.elcenters[:,0],self.elcenters[:,1]
            self.trg = mptri.Triangulation(xe,ye) #,triangles=elements[:,-3:]) # be careful the numbering can be differnet from gmesh
        if mshType == 1: # voronoi
            vor = Voronoi(points)
            msh = myVor(self,vor)
            msh.transformVor(bbox,npts0)
            xn,yn = self.nodes[:,0],self.nodes[:,1]
            self.trg = mptri.Triangulation(xn,yn,triangles=elements[:,-3:]) # be careful the numbering can be differnet from gmesh
        if self.core.addin.getDim()=='3D': 
            med = self.core.dicaddin['3D']
            top = [float(x) for x in med['topMedia']]
            bot = top[1:]; bot.append(float(med['zmin']))
            thick = [top[i]-bot[i] for i in range(len(top))]
            self.core.addin.get3D()
            self.add3d(getNmedia(self.core),thick) # !!! media = layer up to now

    def getCenters(self): return self.nodes[:,0],self.nodes[:,1]
    def getNumber(self,typ): 
        if typ == 'nodes': return len(self.nodes)
        if typ == 'elements': return self.ncell_lay
        if typ == 'tot_elements': return self.ncell
        
    def addSidesPts(self,pts,bbox,dens):
        #left
        iside = where(pts[:,0]==bbox[0])[0]
        p_out = pts[iside]*1;p_out[:,0] -= dens/3
        pts[iside,0] += dens/2
        pts = r_[pts,p_out]
        #rigth
        iside = where(pts[:,0]==bbox[1])[0]
        p_out = pts[iside]*1;p_out[:,0] += dens/3
        pts[iside,0] -= dens/2
        pts = r_[pts,p_out]
        #bottom
        iside = where(pts[:,1]==bbox[2])[0]
        p_out = pts[iside]*1;p_out[:,1] -= dens/3
        pts[iside,1] += dens/2
        pts = r_[pts,p_out]
        #top
        iside = where(pts[:,1]==bbox[3])[0]
        p_out = pts[iside]*1;p_out[:,1] += dens/3
        pts[iside,1] -= dens/2
        pts = r_[pts,p_out]
        return pts
        
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
            for i in range(l0/50):
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
            if il==0: 
                cneighb1.extend([r_[array(cn),i+ncell2d] for i,cn in enumerate(cneighb0)])
                cdist1.extend([r_[array(cd), thk[0]/2] for cd in cdist0])
                fahl1.extend([r_[array(fl)*thk[0],carea0[i]] for i,fl in enumerate(fahl0)]) 
            elif il==nlay-1: 
                cneighb1.extend([r_[array(cn)+il*ncell2d,i+(il-1)*ncell2d] for i,cn in enumerate(cneighb0)])
                cdist1.extend([r_[array(cd), thk[-1]/2] for cd in cdist0])
                fahl1.extend([r_[array(fl)*thk[-1],carea0[i]] for i,fl in enumerate(fahl0)]) 
            else :
                cneighb1.extend([r_[array(cn)+il*ncell2d,i+(il-1)*ncell2d,i+(il+1)*ncell2d] for i,cn in enumerate(cneighb0)])
                cdist1.extend([r_[array(cd), thk[il-1]/2, thk[il+1]/2] for cd in cdist0])
                fahl1.extend([r_[array(fl)*thk[il],carea0[i],carea0[i]] for i,fl in enumerate(fahl0)]) 
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
    def __init__(self,parent,vor):
        self.parent,self.nodes,self.xyverts,self.vor = parent,vor.points,vor.vertices,vor
        
    def transformVor(self,bbox,npts):
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
        not : cdist, cneighb, fahl only include the links to cells that are inside 
        the domain an thus contain les elements
        '''
        p = self.parent
        p.cverts,p.carea,p.fahl,p.cneighb,p.cdist=[],[],[],[],[]
        p.ncell,p.ncell_lay,p.xyverts = npts,npts,self.xyverts;#print npts
        p.nodes = self.nodes[:npts,:];#print shape(self.nodes)
        lnp = [len(x) for x in self.vor.regions]
        mx_lnp = max(lnp)
        p.rarray = zeros((npts,mx_lnp)).astype('int')
        riverts = array(self.vor.ridge_vertices,ndmin=2)
        p.nconnect = 0
        for ip in range(npts):
            ir = self.vor.point_region[ip]
            xp,yp = p.nodes[ip]
            iverts = self.vor.regions[ir]
            nv = len(iverts)
            #print ip,xp,yp,iverts,nv
            p.rarray[ip,:nv] = iverts;
            p.rarray[ip,nv:] = iverts[-1]
            if iverts[0]!=iverts[-1]:iverts.append(iverts[0])
            nverts=len(iverts)
            p.cverts.append(iverts)
            # calculates area
            x,y = p.xyverts[iverts,0],p.xyverts[iverts,1]
            area = abs(sum(x[:-1]*y[1:]-x[1:]*y[:-1]))/2
            p.carea.append(area)
            # find the points connected to ip
            idx = list(where(self.vor.ridge_points==ip)[0])
            riv1, rip1 = riverts[idx,:],self.vor.ridge_points[idx,:]
            #print riv1,rip1
            plist, p_in = [],[] # p_in is a list of points that are inside the domain
            for i in range(nverts-1):
                v1,v2 = iverts[i:i+2];#print v1,v2
                r1,r2 = list(where(riv1==v1)[0]), list(where(riv1==v2)[0]);#print r1,r2
                ind =list(set(r1)&set(r2))[0]
                pt = list(rip1[ind]);pt.remove(ip)
                if pt[0]<npts:
                    p_in.append(i)
                    plist.append(pt[0])
            # length of ridge
            ln0 = sqrt((x[1:]-x[:-1])**2+(y[1:]-y[:-1])**2) # length of face
            ln0 = ln0[p_in]
            # finds the distance from pt to face
            dlist = []
            for i in range(nverts-1):
                if i not in p_in: continue
                c = 1
                xym = array([[x[i],y[i]],[x[i+1],y[i+1]]])
                if sum(abs(xym[:,0]))==0: 
                    lcoef =array([1,0]);c=0
                elif sum(abs(xym[:,1]))==0: 
                    lcoef =array([0,1]);c=0
                else : 
                    lcoef = solve(xym,ones((2,1))[:,0])
                dst = abs(dot(array([xp,yp]),lcoef)-c)/sqrt(sum(lcoef**2))
                dlist.append(dst)
            p.cneighb.append(plist)    
            p.fahl.append(ln0)
            p.cdist.append(dlist)
            p.nconnect += len(dlist)
        p.xyverts[:,0] = clip(p.xyverts[:,0],bbox[0],bbox[1])
        p.xyverts[:,1] = clip(p.xyverts[:,1],bbox[2],bbox[3])
        p.elx = [p.xyverts[vts,0] for vts in p.cverts]
        p.ely = [p.xyverts[vts,1] for vts in p.cverts]
        p.carea1 = p.carea*1
        p.elcenters = array(p.nodes)

        #print 'len cv',len(self.cverts);#'nconnect',self.nconnect
