# -*- coding: utf-8 -*-
"""
Created on Mon Jul  4 07:11:43 2022

@author: olivi
"""
import os
from config import *
import matplotlib.tri as mptri
import numpy as np
from numpy import pi,in1d,nonzero,insert
from geometry import *
from scipy import zeros,ones,array,arange,r_,c_,around,argsort,unique,cumsum,where,shape,\
    amin,amax,mod,lexsort


##########################   GMESH    ########################
##########################################################

def gmeshString(core,dicD,dicM):
    """creates a gmesh string to use in gmesh to generate the mesh
    dicD is the domain zones, containing points, lines and domain
    dicM is the media zones dict, they will be used to build the grid"""
    s,p_list,p_link = '',[],[]
    i_domn = dicD['name'].index('domain') # search for the line called domain
    dcoords = dicD['coords'][i_domn][:-1]
    ddens = float(dicD['value'][i_domn])
    ddens = corrMeshDens(dcoords,ddens)
    if isPolyTrigo(dcoords) : dcoords = dcoords[-1::-1] # poly mus tbe clockwise
    mgroup = core.addin.getModelGroup()
    if mgroup=='Opgeo': dens =  str(core.dicval[mgroup+'Flow']['domn.4'][0])
    ldom=len(dcoords)
    # domain points
    for i,pt in enumerate(dcoords):
        s+='Point('+str(i+1)+')={'+str(pt[0])+','+str(pt[1])+',0,'+str(ddens[i])+'}; \n'
        p_list.append(pt)
    npt = i+1;lli = list(range(i+1)) # OA 26/4
    # domain lines
    s+= stringLine(p_link,list(range(ldom)),1,opt='close')
    #other points: present in domain, name shall start by point
    indx = [dicD['name'].index(b) for b in dicD['name'] if b[:5]=='point']
    p_coord = [dicD['coords'][iz][0] for iz in indx] # a list of all coords in zone points... aded [0]
    p_dens = [dicD['value'][iz] for iz in indx] # same for densities
    #add the fault points OA 16/8/20
#    indx = [dicD['name'].index(b) for b in dicD['name'] if b[:5]=='fault']
#    if len(indx)>0:
#        for iz in indx:
#            pts1,pts2 = ptForFaults(dicD['coords'][iz],float(dicD['value'][iz]))
#            p_coord.extend(pts1);p_coord.extend(pts2)
#            p_dens.extend([dicD['value'][iz]]*len(pts1)*2)
    if len(p_coord)>0:
        p_list,p_link,spt,sa = stringPoints(p_list,p_coord,p_dens,npt)
        s+=sa
        npt += len(p_coord)
    # other lines present in the domain
    s1 = ''
    for iz,n in enumerate(dicD['name']):
        if n[:4] in ['line','faul']: # OA 16/8/20 correct for line
            p_coord = dicD['coords'][iz]
            p_dens = dicD['value'][iz]
            p_list,p_link,spt,sa = stringPoints(p_list,p_coord,p_dens,npt)
            s1+=sa
            p_range = list(range(npt,npt+len(p_coord)))
            s1+= stringLine(p_link,p_range,npt)
            npt += len(p_coord)
            lli.extend(p_range) # OA 26/4 npt->nli
    # add the material zones (opgeo), but don't add points twice!!
    isurf,s2 = 2,''
    for iz in range(len(dicM['name'])):
        if (dicM['name'][iz] == 'domain')  or (dicM['name'][iz][0]=='_'): continue
        p_coord = dicM['coords'][iz];#print 'geom 517',iz,p_coord
        x,y = list(zip(*p_coord)); x,y = array(x),array(y)
        #dens = str(min(max(x)-min(x),max(y)-min(y))/1.5)[:6]
        d = sqrt((x[1:]-x[:-1])**2+(y[1:]-y[:-1])**2)
        dens = str(min(d)*2)[:6];#print 'geom 599 dens',dicM['name'],dens
        if len(p_coord)>2: # don't take points, just lines
            p_list,p_link,spt,ss = stringPoints(p_list,p_coord,dens,npt)
            s2+=ss
            p_range = list(range(npt,npt+len(p_coord)))
            s2+= stringLine(p_link,p_range,isurf,'close')
            npt += len(p_coord)
            s2+='Plane Surface('+str(isurf)+')={'+str(isurf)+'}; \n';
            isurf+=1
    a=''
    if isurf >2 : a = ','+','.join([str(x) for x in range(2,isurf)])
    s = s+s1+s2+'Plane Surface(1)={1'+a+'}; \n';
    #s+='Physical Surface(1)={1}; \n';
    for ip in range(npt):
        s+='Point{'+str(ip+1)+'} In Surface{1}; \n'
    for il in lli: # OA 26/4
        s+='Line{'+str(il+1)+'} In Surface{1}; \n'
    return s#+'Recombine Surface{1}; \n'
    
def corrMeshDens(coord,dens): # OA added 28/9/19
    npts=len(coord)
    try: len(dens)
    except TypeError: dens = [dens]*npts
    coo,dns = array(coord,ndmin=2),array(dens)
    dst = sqrt((coo[1:,0]-coo[:-1,0])**2+(coo[1:,1]-coo[:-1,1])**2)
    dstmin = r_[dns[0],minimum(dst[:-1],dst[1:]),dns[-1]]
    dns[dns>dstmin] = dstmin[dns>dstmin]
    return dns
                
def stringPoints(p_list,p_coord,p_dens,istart):    
    '''creates lines of points from coordinates and returns the point list'''
    def isCooInList(coo,lst,eps):
        x,y = list(zip(*lst));x,y = array(x),array(y)
        #d = sqrt((x[1:]-x[:-1])**2+(y[1:]-y[:-1])**2)
        #eps = min(d)/20
        for i,c1 in enumerate(lst):
            if abs(coo[0]-c1[0])<eps and abs(coo[1]-c1[1])<eps: return i
        return -1
        
    s,spt,p_link='','',[]
    x,y = list(zip(*p_list))
    eps = min(max(x)-min(x),max(y)-min(y))/500
    if type(p_dens) != type([5]): 
        p_dens = [p_dens]*len(p_coord)
    for ip,coo in enumerate(p_coord):
        prf,idx ='',isCooInList(coo,p_list,eps)
        if idx>=0: 
            p_link.append([istart+ip,idx]) # make the link btw the present pt and the existing one with same coords
            prf='//' # don't add twice the same point
        p_list.append(coo) # store the point
        if type(coo[0])==type((5,6)): coo = coo[0]
        dens = p_dens[ip]
        s+=prf+'Point('+str(istart+ip+1)+')={'+str(coo[0])+','+str(coo[1])+',0,'+dens+'}; \n'
        spt+= str(istart+ip+1)+','
    return p_list,p_link,spt,s
    
def stringLine(p_link,p_range,il,opt='None'):
    '''creates a line string from a list of points number
    it has to consider pre-existing points'''
    s,pnew,pold='',[],[]
    if len(p_link)>0: pnew,pold = list(zip(*p_link)) # get the new and old ref of the same points
    s+='// p '+str(pnew)+'  '+str(pold)+'\n'
    l_link = []
    #print pnew,pold
    for i in p_range[:-1]:
        prf='';#print i
        if (i in pnew) and (i+1 in pnew):
            if pold[pnew.index(i+1)]==pold[pnew.index(i)]+1: # the line must be in the same order
                l_link.append((i,pold[pnew.index(i)]))
                prf='//' # don't write an existing line
        a,b = ptreplace(pold,pnew,i,i+1)
        s+=prf+'Line('+str(i+1)+')={'+str(a)+','+str(b)+'}; \n'
    if opt=='close':
        a,b = ptreplace(pold,pnew,p_range[-1],p_range[0])
        s+='Line('+str(p_range[-1]+1)+')={'+str(a)+','+str(b)+'}; \n'
    if opt=='None': p_range = p_range[:-1]
    v,lnew,lold='',[],[]
    if len(l_link)>0:lnew,lold = list(zip(*l_link))
    s+='// l'+str(lnew)+'  '+str(lold)+'\n'
    for ip in p_range:
        if ip in lnew: v+=str(lold[lnew.index(ip)]+1)+','
        else : v+=str(ip+1)+','
    s+='Line Loop('+str(il)+')={'+v[:-1]+'};\n'
    return s
    
def ptreplace(pold,pnew,i,j):
    a = i+1
    if i in pnew: a = pold[pnew.index(i)]+1
    b = j+1
    if j in pnew: b = pold[pnew.index(j)]+1
    return a,b

def ptForFaults(poly,dns): # !!! not used finally
    '''create new points around a fault to build the fault later. form a poly
    we create points in the perpendicular direction at a dens distance. In the middle
    of the poly the points are at a dens distance of the two line. Then the interval
    is divided on both sides to have the same nb of points
    '''
    x,y = zip(*poly)
    lc0 = lcoefsFromPoly(poly);c=-1
    a,b = lc0[:1,:],lc0[1:,:]
    A = sqrt(dns**2/(1+a**2/b**2))
    yn1 = y[:-1]-A;xn1=x[:-1]-a/b*A # starting point
    yn1 = r_[yn1,y[:-1]+A];xn1=r_[xn1,x[:-1]+a/b*A]
    yn2 = y[1:]-A;xn2=x[1:]-a/b*A # ending point
    yn2 = r_[yn2,y[1:]+A];xn2=r_[xn2,x[1:]+a/b*A]
    xn1 = c_[xn1,xn2[:,-1]];xn2=c_[xn1[:,0],xn2];xn=(xn1+xn2)/2 # mke the average
    yn1 = c_[yn1,yn2[:,-1]];yn2=c_[yn1[:,0],yn2];yn=(yn1+yn2)/2
    d1 = sqrt((x-xn)**2+(y-yn)**2) # correct for distance for pt sin the middle
    xn=x+(xn-x)*dns/d1;yn=y+(yn-y)*dns/d1
    pts0,pts1=[(xn[0,0],yn[0,0])],[(xn[1,0],yn[1,0])]
    for i in range(len(x)-1):
        d=sqrt((xn[:,i+1]-xn[:,i])**2+(yn[:,i+1]-yn[:,i])**2) # dst from 1 pt to the next
        dv = int(around(mean(d)/dns))
        lx,ly = linspace(xn[0,i],xn[0,i+1],dv),linspace(yn[0,i],yn[0,i+1],dv)
        pts0.extend(zip(lx[1:],ly[1:]))
        lx,ly = linspace(xn[1,i],xn[1,i+1],dv),linspace(yn[1,i],yn[1,i+1],dv)
        pts1.extend(zip(lx[1:],ly[1:]))
    return pts0,pts1
        
def readGmshOut(s,outline=False,nbdy=0,npts=0,liNames=[]): # OA 28/7/20 added nbdy
    '''reads the inside of a gmesh mesh file from the string s
    and returns nnod : nb of nodes, nodes : nodes coordinates, 
    nel : nb of elements, el
    '''
    b = s.split('$Nodes')[1]
    c = b.split('$EndNodes')[0]
    c1 = c.split('\n')       
    c2 = [x.split() for x in c1[2:-1]];#print len(c2)
    nodes = array(c2,dtype='float')
    nodes[:,0] = nodes[:,0]-1 # node number start from 0
    # elements
    b = s.split('$Elements')[1];#print len(b)
    c = b.split('$EndElements')[0];#print len(c)
    c1 = c.split('\n')
    c2 = [x.split() for x in c1[2:] if len(x.split())==8];#elements have a length 8
    elements = array(c2,dtype='int');#print 'ogModel 48',arr
    elements = elements-1
    elements = elements[:,[0,3,4,5,6,7]] # 2nd column is useless
    if outline: 
        dcout = {} # OA 15/8/20 modif to dict
        c3 = array([x.split() for x in c1[2:] if len(x.split())==7],dtype='int') # OA 28/7/20
        lines = unique(c3[:,4]); count=0;liout = []
        for l in lines: 
            if l<=nbdy: dcout['bc'+str(l-1)] = c3[c3[:,4]==l,-2:]-1 #domain
            else : #true lines
                arr = c3[c3[:,4]==l,-2:]-1
                if len(liout)==0: liout.append(arr)
                else :
                    if liout[-1][-1,1] == arr[0,0]: # this is the same poly
                        liout[-1] = r_[liout[-1],arr]
                    else :
                        liout.append(arr)
                        count += 1
        for i in range(len(liout)):
            dcout[liNames[i]] = liout[i]
        return nodes,elements,dcout
    else :
        return nodes,elements

def createTriangul(obj):
    '''get the triangulation for matplotlib representation'''
    xnds,ynds = obj.nodes[:,1],obj.nodes[:,2]
    trg = mptri.Triangulation(xnds,ynds,triangles=obj.elements[:,-3:]) # be careful the numbering can be differnet from gmesh
    obj.trg = trg
    nel = len(trg.triangles)
    obj.elx = trg.x[trg.triangles]
    obj.ely = trg.y[trg.triangles]
    obj.elcenters = zeros((nel,2))
    obj.elcenters[:,0] = mean(obj.elx,axis=1)
    obj.elcenters[:,1] = mean(obj.ely,axis=1)
    return obj

class unstructured:
    
    def __init__(self,core):
        self.core=core
        
    def buildMesh0(self,modName,opt):
        '''opt=old the .msh file is used if it is found, new : it is rewritten'''
        self.core.addin.mesh = self
        self.nlay = getNlayers(self.core)
        dct = self.core.diczone[modName].dic
        if modName=='OpenFlow' and 'dis.2' not in dct.keys(): 
            self.core.dicval[modName]['dis.1'][0]=0;
            #self.points,self.elts,self.dicD,self.dcoo1,self.dicFeats = [],[],{},[],{}
        mshType = self.core.getValueFromName(modName,'MshType')
        if mshType>1: # case of true unstructured grid built through gmesh
            fmsh = self.core.fileName+'_out.msh'
            os.chdir(self.core.fileDir)
            if fmsh not in os.listdir(os.getcwd()): opt = 'new'
            if modName=='Modflow' : dicD = dct['disu.3'] # where the domain boundaries are 
            elif modName=='OpenFlow' : dicD = dct['dis.2'] # where the domain boundaries are 
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
                os.system(bindir+'gmsh gm_in.txt -2 -smooth 5 -algo front2d -o '+fmsh) #meshadapt
            f1 = open(fmsh,'r');s = f1.read();f1.close();print('gmesh file read')
            nbdy = len(dcoo1) # 15/8/20
            nodes,elements,self.dicFeats = readGmshOut(s,outline=True,nbdy=nbdy,liNames=liNames); # OA 15/8/20
            self.dicFeats['points']=list(range(nbdy,nbdy+npts)) # OA 15/8/20
            self.points,self.elts = nodes[:,1:3],elements[:,-3:]
            self.dicD,self.dcoo1 = dicD,dcoo1
            
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
        
    def addMeshVects(self):
        ''' add elxa and elya, a long list format for elx and ely '''
        self.ncell,self.ncell_lay = len(self.elx),len(self.elx)
        ic,self.idc=0,[]
        for a in self.elx: self.idc.append([ic,ic+len(a)]);ic+=len(a)
        self.elxa,self.elya=[],[]
        for a in self.elx: self.elxa.extend(a)
        for a in self.ely: self.elya.extend(a)
        
    def makeBC(self,modName):
        '''
        create an array bcindx that gather the cell number and the 
        '''
        if self.core.getValueFromName(modName,'MshType')>0:
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
        
    def getPointsFaces(self):
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
    
    def pointsHexaForVTK(self,modName,points,fcup):
        '''
        core is the model,fDir the major openfoam folder
        '''
        core,mesh,nlay = self.core,self,self.nlay
        nplay = shape(points)[0] # nb of points for each layer
        ncell = len(fcup) # nb of cells per layer
        # below approximative, just for flat layers
        if core.addin.getDim()=='3D':
            if 'importArray' in core.dictype[modName]['dis.6']:
                zlist = self.getZfromPoints(modName,points)
            else :
                zlist = [float(a) for a in core.addin.get3D()['topMedia']]
                zlist.append(float(core.addin.get3D()['zmin']))
            zlist = zlist[-1::-1] # layers are reversed, start from 0 in OpF
        else : 
            zlist = [core.dicval[modName]['dis.7'][0],core.dicval[modName]['dis.6'][0]]
        # points 
        nbp = shape(points)[0]*(nlay+1)
        sp = ''
        x0,y0 = amin(points[:,0]),amin(points[:,1])
        strt = 0
        for i,z in enumerate(zlist):
            if type(z)==type(0.): z = [z]*nplay
            coo = zeros((shape(points)[0],3))
            coo[:,0],coo[:,1],coo[:,2] = points[:,0]-x0,points[:,1]-y0,z
            sp += '\n'.join([' '.join(x) for x in coo.astype('str')]) 
            sp += '\n'
        
        lh = [];lnh=[];idx=[];cnt=0
        # internal faces for 1st nlay-1 # This part is quite slow!! because top face added there
        for il in range(nlay): 
            for i in range(ncell):
                #sh += str(len(fcup[i])-1)+'('+' '.join([str(ip+nplay*(il+1)) for ip in fcup[i][:-1]])+ ')\n'
                #lh.append([ip+nplay*(il+1) for ip in fcup[i][:-1]])
                f = fcup[i]
                if len(fcup[i])<8: 
                    lh.append(r_[nplay*il+f[:-1],nplay*(il+1)+f[:-1]])
                    lnh.append((len(f)-1)*2)
                    idx.append(cnt);cnt+=1
                else : # we need to make two poly (VTK does not support more than 6 points)
                    lh.append(r_[nplay*il+f[:5],nplay*(il+1)+f[:5]])
                    lh.append(r_[nplay*il+f[4:],nplay*(il+1)+f[4:]])
                    lnh.extend([10,(len(fcup[i])-4)*2])
                    idx.extend([cnt,cnt]);cnt+=1
        return (nlay+1)*nplay,sp,len(lnh),lnh,lh,idx
    
    def getZfromPoints(self,modName,points):
        '''
        get the z coords of the points grid, cannot be done through getValueLong or zblock
        as the points are not at the center of the cell
        this is called only when arrays are present (interpolated case not treated)
        as z will be reversed, here the list of points for each layer starts at the top
        '''
        core = self.core;lzout=[];zb=core.Zblock;dzmin=(amax(zb)-amin(zb))/500
        core.lcellInterp = [] # to reset the values where to search (default cell centers)
        grd = core.addin.getFullGrid();intp,ysign,zdx,zdy=False,0,None,None
        xx,yy,intp,z0 = points[:,0],points[:,1],False,1e6
        fDir = self.core.fileDir
        for i in range(self.nlay): 
            fNameExt = core.dicarray[modName]['dis.6'][i] #tops
            if '.' not in fNameExt: # not an array
                z1 = (xx*0+1)*float(core.dicval[modName]['dis.6'][i])
            else:
                if fNameExt[-3:] == 'var' : ysign,zdx,zdy,zgrd = core.importGridVar(self.core.fileDir,fNameExt) # OA 13/6/20 add ysign
                elif fNameExt[-3:] == 'asc' : zgrd =core.importAscii(fDir,fNameExt)
                else : zgrd = loadtxt(fDir+fNameExt)
                if ysign == -1 : 
                    zgrd = zgrd[-1::-1];zdy = zdy[-1::-1]
                z1 = linIntpFromGrid(core,grd,zgrd,xx,yy,intp,zdx,zdy)
            lzout.append(minimum(z1,z0-dzmin)) # to avoid negative thickness
            z0 = z1*1
        # for the bottom
        if 'dis.7' in core.dicarray[modName].keys(): #bottom
            fNameExt = core.dicarray[modName]['dis.7'][-1]
            if fNameExt[-3:] == 'var' : ysign,zdx,zdy,zgrd = core.importGridVar(self.core.fileDir,fNameExt) # OA 13/6/20 add ysign
            elif fNameExt[-3:] == 'asc' : zgrd =core.importAscii(fDir,fNameExt)
            else : zgrd = loadtxt(fDir+fNameExt)
            if ysign == -1 : 
                zgrd = zgrd[-1::-1];zdy = zdy[-1::-1]
            z1 =linIntpFromGrid(core,grd,zgrd,xx,yy,intp,zdx,zdy)
        else :
            z1 = (xx*0+1)*float(core.dicval[modName]['dis.7'][-1])
        lzout.append(minimum(z1,z0-dzmin))# removed [::-1]
        core.lcellInterp = [] # to reset the values where to search (default cell centers)
        return lzout
    
    def stringVTKgeom(self):
        ''' returns a string for the geometry part of the VTK file, with triangles'''
        def fmtlist(a):
            return str(a).replace('[','').replace(']','').replace(',','')
        pts=self.elcenters
        zb=self.core.Zblock
        zmid=(zb[:-1]+zb[1:])/2;nlay,npts=shape(zmid)
        # point coord list
        sp = ''
        for il in range(nlay):
            sp += fmtlist(c_[pts,zmid[il]])+'\n'
        #wedge cells
        tri = self.trg.triangles
        ntri = len(tri)
        cells = c_[ones((ntri,1))*6,tri,tri+npts]
        for i in range(1,nlay-1):
            c1 = c_[ones((ntri,1))*6,tri+i*npts,tri+(i+1)*npts]
            cells = r_[cells,c1]
        ncells = len(cells)  
        s='# vtk DataFile Version 2.0\n'
        s+='Unstructured Grid Example\nASCII\nDATASET UNSTRUCTURED_GRID\n'
        s+='POINTS '+str(npts)+' float\n'+sp
        s+='\nCELLS '+str(ncells)+' '+str(ncells*3)+'\n'
        s+= fmtlist(cells)+'\n'
        s+= '\nCELL_TYPES '+str(ncells)+'\n' 
        s+= '\n'.join((ones(ncells)*6).astype('int').astype('str'))
        s+='\n'
        return s
    
    def writeVTKgeom_back(self,npt,sp,nh,lnh,lh):
        ''' returns a string for the geometry part of the VTK file, old one with polygons'''
        def fmtlist(a):
            return str(a).replace('[','').replace(']','').replace(',','')
        #4 pts/fcup: 12; 3:13; 5: 15; 6: 16
        lnh2=array(lnh)
        lnh2[lnh2==12]=16;lnh2[lnh2==6]=13;lnh2[lnh2==8]=12;lnh2[lnh2==10]=15;
            
        s='# vtk DataFile Version 2.0\n'
        s+='Unstructured Grid Example\nASCII\nDATASET UNSTRUCTURED_GRID\n'
        s+='POINTS '+str(npt)+' float\n'+sp
        s+='\nCELLS '+str(nh)+' '+str(sum(lnh)+nh)+'\n'
        for i in range(nh):
            s+=str(len(lh[i]))+' '+fmtlist(lh[i])+'\n'
        s += '\nCELL_TYPES '+str(nh)+'\n' 
        for i in range(nh): 
            s+= str(lnh2[i])+'\n'
        s+='\n'
        return s

    def writeVTKgeomXML(self):
        pass
        '''
        s='<?xml version="1.0"?>\n'
        s+='<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian" header_type="UInt32" compressor="vtkZLibDataCompressor">\n'
        s+='<UnstructuredGrid><Piece NumberOfPoints="'+str(np*2)+'" NumberOfCells="'+str(nh)+'\">\n'
        #s+= <PointData></PointData><CellData></CellData>
        s += '<Points>\n<DataArray type="Int32" Name="Points" NumberOfComponents="3" format="ascii" >\n'
        s += sp
        s += '</DataArray>\n  </Points>\n  <Cells>\n  '
        s += '<DataArray type="Int64" Name="connectivity" format="ascii" >\n'
        s += '\n'.join([fmtlist(a) for a in lh2])
        s += '</DataArray>\n'
        s += '<DataArray type="Int64" Name="offsets" format="ascii">\n'
        csum = cumsum(lnh)
        s += ' '.join([str(a*2)+' ' for a in csum])+'\n'
        s += '</DataArray>\n'
        s += '<DataArray type="UInt8" Name="types" format="ascii" RangeMin="5" RangeMax="7">\n'
        # 5 faces 15 penta prism, 6 faces 16 hexa prism
        s += ' '.join([str(a+10)+' ' for a in lnh])+'\n'
        s += '</DataArray>\n </Cells>\n </Piece>\n </UnstructuredGrid>\n </VTKFile>'
        f1=open('test.vtu','w');f1.write(s);f1.close()
        '''
        
                        
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
        theta = arctan(dy/dx) # angle from point 1 to 2
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
        p.M = M
        ip1 = M[:,1].astype('int')
        indsplit = where(ip1[1:]-ip1[:-1]==1)[0]+1 # this is the nb of points for each elt
        p.nconnect = shape(M)[0]
        p.cneighb = split(M[:,2].astype('int'),indsplit);nc=len(p.cneighb)
        lvx = split(M[:,3],indsplit);p.elx = [r_[v,v[0]] for v in lvx]
        lvy = split(M[:,4],indsplit);p.ely = [r_[v,v[0]] for v in lvy]
        p.cdist = split(M[:,5],indsplit)
        xmd1,ymd1 = split(M[:,6],indsplit),split(M[:,7],indsplit)
        # redo differently elx,ely and length for the poly at the boundaries (the rest is ok)
        indx = [0]*(max(line_dom)+1) # indx contains the face number at the bdy in line_dom
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
            angl -= pi/2
            #angl[angl<0]=2*pi+angl[angl<0]
            p.angl.append(angl)
        for ip1 in range(p.ncorn): # OA 26/4/20 changed to indx
            p.fahl[ip1] = delete(p.fahl[ip1],indx[ip1])
            p.angl[ip1] = delete(p.angl[ip1],indx[ip1])
        for ip1 in line_dom1: # OA 26/4/20 changed to indx
            p.fahl[ip1] = delete(p.fahl[ip1],indx[ip1])
            p.angl[ip1] = delete(p.angl[ip1],indx[ip1])
        la = [abs(sum(p.elx[i][:-1]*p.ely[i][1:]-p.elx[i][1:]*p.ely[i][:-1]))/2 for i in range(nc)]
        p.carea = array(la);p.indx = indx;p.ncell=len(p.elx)
        p.elcenters = p.nodes*1 # OA 9/1/21
        for i in line_dom: p.elcenters[i] = (mean(p.elx[i]),mean(p.ely[i]))
        l = []
        for k in p.dicFeats.keys():
            if k[:2]=='bc': l.append(unique(p.dicFeats[k]))
        for i in range(p.ncorn): p.dicFeats['bc_cell'+str(i)]=l[i]

    def putInDomain(self,elts,dcoo,M,dicFeats):  # added 25/4/20 some poits can be outside domain
        ''' put points out of the domain back into it '''
        poly=dcoo*1;poly.append(dcoo[0]);npl=len(poly)
        lcoefs=lcoefsFromPoly(poly)
        iIn=array(pointsInPoly(M[:,3],M[:,4],poly,lcoefs))
        lout = where(iIn==0)[0];nl=len(lout)
        # determine the line of interest
        lines=[]
        for k in dicFeats.keys():
            if (k[:2]=='bc')&(k[:4]!='bc_c'):
                lines.append(unique(dicFeats[k]))
        it0 = unique(M[lout,0]).astype('int') # index of the triangle on which is the outlier
        # on whcich bdy is the triangle
        lc1=[]
        for i in it0:
            for ip in range(npl-1):
                if sum(in1d(lines[ip],elts[i]))==2: 
                    lc1.append(lcoefs[:,ip])
        # set the point symetric to the line (ax+by+c=0) ici ax+by=1
        # y’=(-2abx-y(b^2-a^2)-2bc)/(b^2+a^2) et x’=x-a(y-y’)/b
        for i,it in enumerate(it0):
            a,b = lc1[i];c=-1
            x,y = M[M[:,0]==it,3],M[M[:,0]==it,4]
            yn = (-2*a*b*x-y*(b**2-a**2)-2*b*c)/(b**2+a**2+1e-12) #OA 15/4/21 added 1e-6
            M[M[:,0]==it,4] = yn
            M[M[:,0]==it,3] = x-a/(b+1e-12)*(y-yn) #OA 15/4/21 added 1e-6
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
    
