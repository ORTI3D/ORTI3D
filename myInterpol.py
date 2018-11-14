# -*- coding: utf-8 -*-
"""
Created on Thu Jan 26 05:19:46 2017

@author: oatteia
"""

#!/usr/bin/env python
import numpy as np
import os
from numpy.linalg import solve
from scipy.stats import tvar as variance
from scipy.spatial import Delaunay
from .geometry import *
from .config import *
#from .pykrige.ok import OrdinaryKriging
#import pykrige.kriging_tools as kt


def linIntpFromGrid(core_grd,z_grid,xx,yy): # added OA 22/5/17
    '''linear interpolation from a regular grid
    that has exactly the same size as the current domain
    xx,yy are the position of the point where the z vlaue is searched for
    x0,y0 are the origin of the domain
    dx,dy are the cell width of the regular grid'''
    x0,x1,y0,y1 = core_grd['x0'],core_grd['x1'],core_grd['y0'],core_grd['y1']
    nr,nc = shape(z_grid)
    dx,dy = (x1-x0)/nc,(y1-y0)/nr
    nr_out =0
    if len(shape(xx))> 1 : 
        nr_out,nc_out = shape(xx)
        xx,yy = ravel(xx),ravel(yy)
    vi,vj = floor((yy-(y0+dy/2))/dy),floor((xx-(x0+dx/2))/dx)
    vi = clip(vi,0,nr-1); vj=clip(vj,0,nc-1)
    dxx, dyy = (xx-(x0+dx/2+dx*vj))/dx,(yy-(y0+dy/2+dy*vi))/dy
    vi,vj,vi1,vj1 = list(vi),list(vj),list(clip(vi+1,0,nr-1)),list(clip(vj+1,0,nc-1))
    dzx = z_grid[vi,vj1]-z_grid[vi,vj]
    dzy = z_grid[vi1,vj]-z_grid[vi,vj]
    dzxy = z_grid[vi,vj]+z_grid[vi1,vj1]-z_grid[vi,vj1]-z_grid[vi1,vj]
    zz = z_grid[vi,vj] + dzx*dxx + dzy*dyy + dzxy*dxx*dyy
    if nr_out>0: zz=reshape(zz,((nr_out,nc_out)))
    return clip(zz,amin(z_grid),amax(z_grid))

def invDistance(xpt,ypt,zpt,x,y,power=1.):
    """inverse distance interpolation on x,y points, and reference points
    xpt,wpt,zpt, with a given power"""
    n=len(x);z0=x*0.;#l=len(xpt)
    for i in range(n):
        d=maximum(sqrt((x[i]-xpt)**2.+(y[i]-ypt)**2.),1e-4)
        lb=1./(d**power)
        lb=lb/sum(lb)
        z0[i]=sum(zpt*lb)
    return z0
    
######################   kriging  ############################ 
def krige(xpt,ypt,zpt,rg,x,y,vtype='spher'):
    """ krige function to interpolate over a vector of points of x,y coords
    using the base points xpt ypt and the vario distance rg (range)"""
    #print 'geom kr l 522',len(xpt),rg
    n = len(x)
    z0=zeros(n)
    vari=0.5*variance(zpt);#print 'vari',vari
    def gamma(vtype,d,rg): # OA 25/10/18
        if vtype in ['spher',None]: gam = 3/2.*d/rg-1/2.*(d/rg)**3.;gam[d>rg]=1
        elif vtype=='gauss': gam=exp(-(d/rg)**2)
        elif vtype=='gauss3': gam=exp(-(d/rg)**3)
        return gam
    for i in range(n):
        x0=x[i];y0=y[i];
        d=sqrt((x0-xpt)**2+(y0-ypt)**2);
        ind=argsort(d)
        x1,y1,z1=xpt[ind[:16]],ypt[ind[:16]],zpt[ind[:16]];# select closest 16 points
        xm1,ym1=meshgrid(x1,y1)
        d=sqrt((xm1-transpose(xm1))**2.+(ym1-transpose(ym1))**2.)
        gam=gamma(vtype,d,rg) # OA 25/10/18
        l1=len(x1)+1
        A=ones((l1,l1)); #starting the A matrix
        A[:-1,:-1]=gam*vari
        A[-1,-1]=0.
        d = sqrt((x0-x1)**2.+(y0-y1)**2.)
        gam=gamma(vtype,d,rg) # OA 25/10/18
        B = ones((l1,1))
        B[:-1,0]=gam*vari 
        lb=solve(A,B)
        z0[i]=sum(z1*transpose(lb[:-1]))
    #z0 = clip(z0,min(zpt)*0.9,max(zpt)*1.1)
    return z0
    
def krige_OK(xpt,ypt,zpt,rg,x,y,vtype='spher'):
    vparms = [0.5*variance(zpt),rg,0] # sill, range, nugget
    OK = OrdinaryKriging(xpt,ypt,zpt, variogram_model='exponential',
        variogram_parameters=vparms, verbose=False, enable_plotting= False)#'spherical'
    zvalues, sigmasq = OK.execute('points', x, y)
    return zvalues
    
import itertools
from numpy.linalg import lstsq
def polyfit2d(x, y, z, order=3):
    ncols = (order + 1)**2
    G = zeros((x.size, ncols))
    ij = itertools.product(list(range(order+1)), list(range(order+1)))
    for k, (i,j) in enumerate(ij):
        G[:,k] = x**i * y**j
    m, _, _, _ = lstsq(G, z)
    return m

def polyval2d(x, y, m):
    order = int(sqrt(len(m))) - 1
    ij = itertools.product(list(range(order+1)), list(range(order+1)))
    z = zeros_like(x)
    for a, (i,j) in zip(m, ij):
        z += a * x**i * y**j
    return z
###################################  simulation using GSLIB  ################
def makeSimul(core,data,asearch,zbounds=(-3,3)):
    '''runs a simulation with sgsim in the simul folder, using initial heads
    data cols : x,y,head,conc
    zbounds does not seem to change anything
    '''
    fdir = core.fileDir+os.sep
    x,y,Z,asearch = data[:,0],data[:,1],data[:,2],str(asearch)
    zb0,zb1 = zbounds;zb0,zb1 = str(zb0)+' ',str(zb1)+' '
    npts = len(x)
    meanZ,stdZ = mean(Z),np.std(Z)
    v = (Z-meanZ)/stdZ
    # write data.txt
    s='title \n3 \nX \nY \nVar \n'
    a = ['%10.5e %10.5e %10.5e '%(x[i],y[i],v[i]) for i in range(npts)]
    s += '\n'.join(a)
    f1 = open(fdir+'sim_data.txt','w');f1.write(s); f1.close()
    #write sgsim parameters
    grd = core.dicaddin['Grid']
    x0,x1,dx = float(grd['x0']),float(grd['x1']),float(grd['dx'][0])
    y0,y1,dy = float(grd['y0']),float(grd['y1']),float(grd['dy'][0])
    nx,ny = int((x1-x0)/dx),int((y1-y0)/dy)
    s = 'START OF PARAMETERS\n'+fdir+'sim_data.txt \n' #input file\n'   # datafl: the input data file
    s += '1 2 0 3 0 0\n' #column indexes 
    s += '-4.0  4.0  \n' #limit to ignore 
    s += '0      \n' #transformation? 
    s += fdir+'sim_transf.out \n' # output file for transfo
    s += '0      \n' #source of transformation 
    s += fdir+'sim_transf.in \n' #file with data for transfo
    s += '0 0    \n' #cols in that file 
    s += zb0+' '+zb1+'   \n' #zmin and zmax 
    s += '1 '+zb0+' \n' #low tail \n'  
    s += '1 '+zb1+'  \n' #high tail \n'  
    s += '2      \n' #debug level\n'  
    s += fdir+'sim_debug.txt \n' #   dbgfl: the file for the debugging output.
    s += fdir+'sim_results.txt \n' #  outfl: the output grid 
    s += '1      \n' #nb of simulations \n' #  nsim: the number of simulations to generate.
    s += str(nx)+' '+str(x0)+' '+str(dx)+' \n' #nx,xmin,dx  
    s += str(ny)+' '+str(y0)+' '+str(dy)+' \n' #ny,ymin,dy  
    s += '1 0. 1.   \n' #nz,zmin,dz  
    s += str(int(rand()*1e7))+'   \n' #random seed 
    s += '0 20      \n' #min and max nb points 
    s += '5         \n' #nb of max prev simulated to use .
    s += '0         \n' #two part search flag 
    s += '1 3       \n' #multi grid \n'  
    #s += '0 \n' #  nmult: the number of multiple grid refinements
    s += '0         \n' # nb of octants 
    s += asearch+' '+asearch+' 1.0 \n' #search radii 
    s += '0.0 0.0 0.0 \n' #search angles 
    s += asearch+' '+asearch+' 1.0 \n' #search covar \n'  
    s += '0         \n' # krige type 
    #s += '0 \n' #   rho: correlation coefficient 
    s += 'sim_var.txt  \n' #file for local var mean 
    s += '0         \n' #column in this file\n'
    s += '1 0.02    \n' #nbstruct and nugget
    s += '1 0.98 0.0 0.0 0.0  \n' #stru nb, sill, angles  
    s += asearch+' '+asearch+' 1.0    \n' #ranges in 3 directions \n' #
    f1=open(fdir+'sgsim.par','w');f1.write(s);f1.close()
    # run the simulation
    os.system(core.baseDir+os.sep+'bin'+os.sep+'sgsim.exe '+fdir+'sgsim.par')
    m0 = loadtxt(fdir+'sim_results.txt',skiprows=3)
    m = reshape(m0,(ny,nx))
    Z1 = m*stdZ+meanZ
    return Z1

###################################  VORONOI ##################""
def getPolyList(listP,xb,yb):
    polylist = []
    #add points on sides
    y=arange(yb[0],yb[1],yb[1]/4);x=zeros(4)+xb[0] # left
    x=r_[x,arange(xb[0],xb[1],xb[1]/4)];y=r_[y,zeros(4)+yb[1]] # top
    y=r_[y,arange(yb[1],yb[0],-yb[1]/4)];x=r_[x,zeros(4)+xb[1]] # right
    x=r_[x,arange(xb[1],xb[0],-xb[1]/4)];y=r_[y,zeros(4)+yb[0]] # botm
    listP.extend(list(zip(x,y)))
    delauny = Delaunay(listP)
    segments = voronoi(delauny)
    for i in range(len(listP[:-16])): 
        p = findPoly(delauny,segments,listP,i,xb,yb);#print i,p
        polylist.append(p)
    return polylist
    

def findPoly(delauny,segments,listP,ipt,xb,yb):
    """finds a polygon around a point identified in the vt series (vertices)"""
    vt = delauny.vertices
    ltri = where(sum(vt==ipt,axis=1))[0]
    lseg = [] # will contain arrays (1,4) for each segment : x,y 1st pt and 2nd pt
    for tri in ltri:
        nghb = delauny.neighbors[tri]
        if -1 in nghb: continue
        for i,n in enumerate(nghb):
            if ipt not in vt[n]: continue
            s = segments[tri*3+i]
            seg=r_[s[0],s[1]]
#            if any(seg[:,[0,2]]>xb[1]) or any(seg[:,[0,2]]<xb[0]) or any(seg[:,[1,3]]>yb[1]) or any(seg[:,[1,3]]<yb[0]):
#                return None
            if any(seg[[0,2]]>xb[1]) or any(seg[[0,2]]<xb[0]) or any(seg[[1,3]]>yb[1]) or any(seg[[1,3]]<yb[0]):
                return None
            flag = 0
            for k in range(len(lseg)): 
                if all(sort(lseg[k])==sort(seg)): flag = 1# already here
            if flag == 0: lseg.append(seg)
    #find 1st point
    npt = len(lseg)
    if len(lseg)==0: return None
    poly=[lseg[0][:2],lseg[0][2:]];pt=poly[1];current = 0;nit = 0
    while (len(poly)<npt) and (nit<npt):
        nit += 1;
        for i in range(npt):
            if all(lseg[i][:2]==pt) and i!=current:
                pt=lseg[i][2:];poly.append(pt);current=i;
            if all(lseg[i][2:]==pt) and i!=current:
                pt=lseg[i][:2];poly.append(pt);current=i;
    if len(poly)<3: return None
    poly.append(poly[0])
    return poly

def voronoi(delauny):
    triangles = delauny.points[delauny.vertices]
    circum_centers = np.array([triangle_csc(tri) for tri in triangles])
    segments = []
    for i, triangle in enumerate(triangles):
        circum_center = circum_centers[i]
        for j, neighbor in enumerate(delauny.neighbors[i]):
            if neighbor != -1:
                segments.append((circum_center, circum_centers[neighbor]))
            else:
                ps = triangle[(j+1)%3] - triangle[(j-1)%3]
                ps = np.array((ps[1], -ps[0]))

                middle = (triangle[(j+1)%3] + triangle[(j-1)%3]) * 0.5
                di = middle - triangle[j]

                ps /= np.linalg.norm(ps)
                di /= np.linalg.norm(di)

                if np.dot(di, ps) < 0.0:
                    ps *= -1000.0
                else:
                    ps *= 1000.0
                segments.append((circum_center, circum_center + ps))
    return segments

def triangle_csc(pts):
    rows, cols = pts.shape

    A = np.bmat([[2 * np.dot(pts, pts.T), np.ones((rows, 1))],
                 [np.ones((1, rows)), np.zeros((1, 1))]])

    b = np.hstack((np.sum(pts * pts, axis=1), np.ones((1))))
    x = np.linalg.solve(A,b)
    bary_coords = x[:-1]
    return np.sum(pts * np.tile(bary_coords.reshape((pts.shape[0], 1)), (1, pts.shape[1])), axis=0)

def irreg2mat(core,XP,YP,Z):
    """ transform data from an irregular matrix XP, YP,
    to a regular one using the model grid
    """
#if 3>2:
    grd = core.addin.getFullGrid()
    x0,dx,nx = grd['x0'],grd['dx'][0],grd['nx']
    y0,dy,ny = grd['y0'],grd['dy'][0],grd['ny']
    Z1 = zeros((ny,nx))+0.
    # etendre un peu la grille pour eviter les trous
    [l,c] = shape(XP)
    x1=XP[:,1:]/2+XP[:,:-1]/2;y1=YP[:,1:]/2+YP[:,:-1]/2;z1=Z[:,1:]/2+Z[:,:-1]/2
    x2=zeros((l,2*c-1))+0.;y2=x2*0.;z2=x2*0.
    x2[:,::2]=XP;x2[:,1::2]=x1;y2[:,::2]=YP;y2[:,1::2]=y1;z2[:,::2]=Z;z2[:,1::2]=z1
    x3=x2[:,1:]/2+x2[:,:-1]/2;y3=y2[:,1:]/2+y2[:,:-1]/2;z3=z2[:,1:]/2+z2[:,:-1]/2
    x4=zeros((l,4*c-3))+0.;y4=x4*0.;z4=x4*0.
    x4[:,::2]=x2;x4[:,1::2]=x3;y4[:,::2]=y2;y4[:,1::2]=y3;z4[:,::2]=z2;z4[:,1::2]=z3
    ix=around(x4/dx);ix=clip(ix,0,nx-1);ix=ix.astype(int)
    iy=around(y4/dy);iy=clip(iy,0,ny-1);iy=iy.astype(int)
    # mettre dans la matrice puis boucher les derniers trous
    put(Z1,iy*nx+ix,z4)
    Z1b=ravel(Z1)
    ind = compress(Z1b==0.,arange(nx*ny))
    iy0 = floor(ind/nx);iy0=iy0.astype(int)
    ix0 = mod(ind,nx);ix0=ix0.astype(int)
    ind1=(iy0-1)*nx+ix0;ind2=(iy0+1)*nx+ix0;ind3=iy0*nx+ix0-1;ind4=iy0*nx+ix0+1
    ind1=clip(ind1,0,nx*ny-1);ind2=clip(ind2,0,nx*ny-1)
    ind3=clip(ind3,0,nx*ny-1);ind4=clip(ind4,0,nx*ny-1)
    nb=sign(take(Z1b,ind1))+sign(take(Z1b,ind2))+sign(take(Z1b,ind3))+sign(take(Z1b,ind4))
    nb=clip(nb,1,4)
    put(Z1,ind,(take(Z1b,ind1)+take(Z1b,ind2)+take(Z1b,ind3)+take(Z1b,ind4))/nb)
    return Z1
