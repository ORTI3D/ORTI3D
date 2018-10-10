if 3>2:
    import ctypes as C
    import numpy
    from scipy import *

    _lib=numpy.ctypeslib.load_library('lib/rflowC_dll','.')
    ti=numpy.ctypeslib.ndpointer(dtype = numpy.int_)
    tf=numpy.ctypeslib.ndpointer(dtype = numpy.float_)
    _lib.calcTube.argtypes = [tf,tf,C.c_int,C.c_int,tf,tf,C.c_int,C.c_int,tf,tf,tf,\
                          tf,tf,tf,tf,tf,ti]
    _lib.calcTube.restype = C.c_void_p

    _lib.ptsLigne.argtypes = [ti,tf,tf,tf,tf,tf,tf,tf,\
                              tf,tf,tf,tf,tf]
    _lib.ptsLigne.restype = C.c_void_p

    _lib.calcPart.argtypes = [ti,tf,tf,tf,tf,\
                          tf,tf,tf,tf,ti]
    _lib.calcPart.restype = C.c_void_p

    _lib.ptsLineP.argtypes = [ti,tf,tf,tf,tf,tf,tf,tf,tf,tf,ti]
    _lib.ptsLineP.restype = C.c_void_p

    _lib.interp2d.argtypes = [ti,tf,tf,tf,tf,tf]
    _lib.interp2d.restype = C.c_void_p

def calcTube1(data,alpha0,vx,vy,xpin,ypin,clim):
#if 3>2:
    for n in ['data','alpha0','vx','vy','xpin','ypin','clim']:
        exec('tmp = numpy.asarray('+n+');'+n+'i = tmp.astype(numpy.double)')
    ny, nx = vx.shape;
    np=len(xpin);nt=int(max(nx,ny)*1.7);ntmx=array([0]);
    tmp = numpy.asarray(ntmx);ntmxo=tmp.astype(numpy.int_);
    for n in ['xp','yp','tp','cu1','cu2']:
        exec(n+'o = numpy.empty([nt,np], dtype=numpy.float_)')
    _lib.calcTube(datai,alpha0i,nx,ny,vxi,vyi,nt,np,xpini,ypini,climi,\
              xpo,ypo,tpo,cu1o,cu2o,ntmxo)
    nto=ntmxo[0]-1;#print 'nto',nto,xpo[:nto,0],ypo[:nto,0]
    return xpo[:nto,:],ypo[:nto,:],tpo[:nto,:],cu1o[:nto,:],cu2o[:nto,:]

def ptsLigne(ind,tet,dc,xp,yp,tp,cu1,cu2):
#if 3>2:
    tmp = numpy.asarray(ind);indi=tmp.astype(numpy.int_);
    for n in ['tet','dc','xp','yp','tp','cu1','cu2']:
        exec('tmp = numpy.asarray('+n+',dtype=numpy.double);'+n+'i = tmp.astype(numpy.float_)')
    ntmx,nt,np=ind;
    for n in ['xp','yp','tp','cu1','cu2']:
        exec(n+'o = numpy.empty([ntmx,np], dtype=numpy.float_)')
    _lib.ptsLigne(indi,teti,dci,xpi,ypi,tpi,cu1i,cu2i,\
              xpo,ypo,tpo,cu1o,cu2o)
    return xpo,ypo,tpo,cu1o,cu2o

def calcParticule(data,vx,vy,clim):
    for n in ['data','vx','vy','clim']:
        exec('tmp = numpy.asarray('+n+');'+n+'i = tmp.astype(numpy.float_)')
    ny, nx = vx.shape;it=1;ndx=1;ndy=1; # ndx pour grille variable
    nt=int(max(nx,ny)*1.7);it=array([1]);
    dims = array([nx,ny,ndx,ndy,nt]);#print data,dims,dx,dy
    tmp = numpy.asarray(dims);dimsi=tmp.astype(numpy.int_);
    tmp = numpy.asarray(it);ito=tmp.astype(numpy.int_);
    for n in ['xp','yp','tp','cu1']:
        exec(n+'o = numpy.empty([nt], dtype=numpy.float_)')
    _lib.calcPart(dimsi,datai,vxi,vyi,climi,\
	xpo,ypo,tpo,cu1o,ito)
    #print 'rfloC part',ito,xpo[:ito],ypo[:ito]
    return xpo[:ito],ypo[:ito],tpo[:ito],cu1o[:ito]

def ptsLigne2(xp,yp,tp,dt,zp=None):
    """ juste pour une ligne"""
#if 3>2:
    rt=1
    if zp==None:
        zp=xp*1;rt=0
    dims=array([len(xp)]);dt=array([dt])
    tmp=numpy.asarray(dims);dimsi=tmp.astype(numpy.int_)
    for n in ['xp','yp','tp','dt','zp']:
        exec('tmp = numpy.asarray('+n+');'+n+'i = tmp.astype(numpy.float_)')
    nt=100;it=ones(100,dtype=int)-2
    for n in ['xp','yp','tp','zp']:
        exec(n+'o = numpy.empty([nt], dtype=numpy.float_)')
    tmp = numpy.asarray(it);ito=tmp.astype(numpy.int_);
    _lib.ptsLineP(dimsi,xpi,ypi,tpi,dti,zpi,xpo,ypo,tpo,zpo,ito)
    itmx=amax(where(ito>-1))
    if rt==0: return xpo[:itmx],ypo[:itmx],tpo[:itmx],ito[:itmx]
    else : return xpo[:itmx],ypo[:itmx],tpo[:itmx],zpo[:itmx],ito[:itmx]

def interp2d(dims,data,xp,yp,zp):
    for n in ['data','xp','yp','zp']:
        exec('tmp = numpy.asarray('+n+');'+n+'i = tmp.astype(numpy.float_)')
    tmp = numpy.asarray(dims);dimsi=tmp.astype(numpy.int_);
    [nx,ny,np]=dims
    zpo=numpy.empty([ny,nx], dtype=numpy.float_)
    _lib.interp2d(dimsi,datai,xpi,ypi,zpi,zpo)
    return zpo

