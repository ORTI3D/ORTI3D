from .config import *
#from matplotlib import pylab as pl
from .myInterpol import *
import matplotlib.tri as mptri
import numpy as np

"""all geometrical operations are performed here
all coordinates are in real world ones. so matrices have index 0 for rows
at x=0. As modflow work the other way round, the transformation are done
by the modflow writer"""

def arr2string1(arr):
    '''only for coordinates'''
    nr,nc = shape(arr)
    #ndigit=4-int(log10(abs(amax(arr,axis=0))));#print 'nd',ndigit
    #arr = around(arr,ndigit)
    a='\n'.join([' '.join(['%+11.4e' %float(arr[j,i]) for i in range(nc)]) for j in range(nr)])
    return a
    
def arr2string1b(arr):
    '''only for coordinates'''
    a='\n'.join([str(list(x)) for x in arr])
    a = a.replace('[','')
    a = a.replace(']','')
    return a
    
def niceformat(val):
    if float(val)<1e-3 or float(val)>1e3: 
        return '%5.1e '%float(val)
    else :
        return str(val)[:5]
    
def makeGrid(core,gIn):
    """ from a dictionnary returned by the addin dialog, creates a grid
    linear or irregular, x0, x1, y0,y1 are kept from the original
    """
    grd={}
    li = ['x0','x1','y0','y1']
    x0,x1,y0,y1=[float(gIn[a]) for a in li]
    grd['x0'],grd['x1'],grd['y0'],grd['y1']=x0,x1,y0,y1;
    dx, dy = calcDx(core,gIn,x0,x1,'dx'),calcDx(core,gIn,y0,y1,'dy')
    grd['dx'],grd['dy'],grd['nx'],grd['ny']=dx,dy,len(dx),len(dy)
    grd['epsilon'] = min((y1-y0),(x1-x0))/1000 # this is the epsilon to snap nodes
    return grd

def calcDx(core,gIn,x0,x1,name): 
    """calculates the grid fixed or variable on one dimension (always called dx in there"""     
    if name == 'dx' and gIn[name]=='fixed': dx=core.dicval['Modflow']['dis.4']
    elif name == 'dy' and gIn[name]=='fixed': dx=core.dicval['Modflow']['dis.5']
    elif len(gIn[name])==1 or gIn[name][1] in ['',' ']: # linear case
        dx = float(gIn[name][0])
        nx = int(round((x1-x0)/dx))
        dx = ones((nx))*dx
    else:
    # variable list dxi values at point xi
        dxIn = gIn[name] # a list of values
        ldx=linspace(0,0,0) # vect vides
        xi,dxi=[],[]
        for i in range(len(dxIn)):
            sxi=dxIn[i].split()
            if len(sxi)>1: 
                xi.append(float(sxi[0]));dxi.append(float(sxi[1]))
        for i in range(1,len(xi)):
            g2 = calcGriVar(xi[i-1],xi[i],dxi[i-1],dxi[i]);#print shape(ldx),shape(g2),g2
            ldx=concatenate([ldx,g2],axis=0)
        dx = ldx*1.
    return dx
    
def calcGriVar(x0,x1,dx0,dx1):
    """calculates variables grid between two points"""
    a=logspace(log10(dx0),log10(dx1),100);#print 'aq calcgriv',a
    n=round(100*(x1-x0)/sum(a))
    if n<=1:
        return [x1-x0]
    else :
        dxv=logspace(log10(dx0),log10(dx1),int(n));#print x0,x1,dx0,dx1,dxv
        if abs((dxv[0]-(x1-x0))) < dxv[0]/1000: return dxv[:1] # added 28/3/17 oa for fixed grids
        dxv[dxv==max(dxv)]+=x1-x0-sum(dxv) # to fit the total
        return dxv
    
def getNmedia(core):
    return len(core.addin.get3D()['topMedia'])
    
def makeLayers(core):
    """returns a list containing for each media the number of layers
    the list layers can have 3 possiblities : 
        - one value the same nb of layers in each media
        - a series of numbers : the nb of layers in each media
        - or two numbers for each media : the relative thickness of the first and last layer
    all oriented from bottom to top"""
    #zmin = core.addin.get3D()['zmin']
    litop = core.addin.get3D()['topMedia']
    nbM = len(litop) # nb of media may be string list
    lilay = core.addin.get3D()['listLayers']*1 # may be string list
    if core.addin.getDim()=='3D':
        if len(lilay)== 1:  
            lilay = [lilay[0]]*getNmedia(core)
        dzL = []
        for im in range(nbM):
            dzL.append([])
            if type(lilay[im]) in [type('a'),type('a')]: 
                if len(lilay[im].split())>1:
                    dz0, dz1 = lilay[im].split(); #[:2]
                    dzL[im] = calcGriVar(0,1,float(dz0),float(dz1))
                    lilay[im] = len(dzL[im])
                else : 
                    lilay[im]=int(lilay[im]); dzL[im] = [1./lilay[im]]*lilay[im]
            else :
                dzL[im] = [1./lilay[im]]*lilay[im]
    elif core.addin.getDim()=='2D':
        lilay=[1]
        dzL = [1.]
    else : # xsection radial
        nx,ny,xvect,yvect = getXYvects(core)
        lilay=[ny];dzL=[1.]
    #print 'geom 87 dzL',lilay,dzL
    return sum(lilay), lilay, dzL  #lilay : nb of layer per media, dzL relative thickness of a layer in a media 

def getNlayersPerMedia(core):    
    nbL, lilay, dzL = makeLayers(core)
    return lilay
def getNlayers(core):
    nbL, lilay, dzL = makeLayers(core)
    return nbL
    
def media2layers(core,media):
    nbL, lilay, dzL = makeLayers(core); #print 'geom', nbL,lilay,dzL,media
    if type(media)==type([5]):
        if len(media)>1:
            lay1 = sum(lilay[:media[0]])
            lay2 = sum(lilay[:media[-1]+1])
        else :# just one media in list
            m = int(media[0])
            lay1 = sum(lilay[:m])
            lay2 = sum(lilay[:m+1]) 
    else : # jsut one media (not in list)
        lay1 = sum(lilay[:media])
        lay2 = sum(lilay[:media+1]) 
    return arange(lay1,lay2).astype(int)   
    
def mediaInlayers(core):
    nbL, lilay, dzL = makeLayers(core)
    lim = [0]
    for nb in lilay: lim.append(lim[-1]+nb)
    return lim[:-1]

def makeZblock(core):
    """to create the block of cells in 3D,
    for modflow it is a 3D matrix
    for unstruct it is 2D"""
    nlay, lilay, dzL = makeLayers(core)
    nx,ny,xvect,yvect = getXYvects(core); 
    mgroup = core.dicaddin['Model']['group']
    #nbM,optionT,optionB = len(lilay),None,None # EV 19/02/20
    nbM = len(lilay)
    #### choise of line depending the type of model # EV 10/02/20
    mgroup = core.dicaddin['Model']['group']
    modName,line = mgroup.split()[0],'dis.6'
    if modName == 'Opgeo': modName += 'Flow'
    if mgroup in ['Modflow series']: lineTop='dis.6'; lineBot='dis.7' # oa 2/2/20 removed usg_rect
    elif mgroup == 'Modflow USG': lineTop='disu.7'; lineBot='disu.8' # OA 1/8/19 changed UNS to USG
    elif mgroup == 'Min3p': 
        modName = 'Min3pFlow'; lineTop='spat.7'; lineBot='spat.8'
    elif mgroup == 'Opgeo': lineTop='domn.5' ; lineBot='domn.6'
    #### stucture of Zblock
    if core.addin.getDim() not in ['Radial','Xsection']: # 2 or 3D case
        if core.addin.mesh==None: # OA 18/7/19 pb for min3p
            Zblock = zeros((nlay+1,ny,nx)) 
        else : #unstructured grids
            if mgroup[0]=='M': #Min3p amd Modfow USG
                Zblock = zeros((nlay+1,core.addin.mesh.getNumber('elements')))
            elif mgroup[0]=='O': # opengeosys
                Zblock = zeros((nlay+1,core.addin.mesh.nnod)) # opengoesys   
        #### extend the number of type at the number of media # EV 10/02/20
        vtypeTop=core.dictype[modName][lineTop]
        if len(vtypeTop)<nbM:vtypeTop.extend([vtypeTop[0]]*(nbM-len(vtypeTop)))
        vtypeBot=core.dictype[modName][lineBot]
        if len(vtypeBot)<nbM:vtypeBot.extend([vtypeBot[0]]*(nbM-len(vtypeBot)))
        #### get the type of model parameter for each media # EV 10/02/20
        intp=[]
        for im in range(nbM): ## im is the media number
            if im<nbM-1:line=lineTop
            else : line=lineBot
            vtype=core.dictype[modName][line][im] ## get the type of entrance data (zone, array, interpolation...) 
            if vtype=='importArray': intp.append(4) 
            elif vtype in ['one_value','zone']: intp.append(3)
            elif vtype=='interpolate' : 
                intp.append(1) 
            else: ## formula 
                formula = core.dicformula[modName][line][im]
                value=core.formExec(formula)
                return value
        #### creates the first media
        i= 0
        top = getTopBotm(core,modName,lineTop,intp[0],0,refer=None,mat='top') # EV 10/02/20 # EV 19/02/20
        #top = getTopBotm(core,modName,lineTop,intp[0],0,optionT,refer=None,mat='top')
        Zblock[0] = top
        #### create the other media and the layers
        for im in range(nbM):
            if im<nbM-1: ## bottom is in top matrix
                line=lineTop
                botm = getTopBotm(core,modName,line,intp[im+1],im+1,refer=top,mat='top')# EV 19/02/20
            else : ## last one, we take the bottom
                line=lineBot
                botm = getTopBotm(core,modName,line,intp[im],im,refer=top,mat='botm') # EV 19/02/20
            ## creates the sublayers
            if lilay[im]==1: ## just one sublayer
                dff = top-botm;
                limit = mean(dff)/500 ## limit the thick not to be 0
                #if type(dzL[im])==type([5]): dzL[im]=dzL[im][0] ## don't understand why i need that
                Zblock[i+1] = top - maximum(dff,limit) #*dzL[im]
                i+=1
            else: ## dzL[im] is a list
                dz = array(dzL[im]);
                for il in range(lilay[im]):
                    dff = top-botm
                    limit = mean(dff)/500
                    Zblock[i+1] = top - maximum(dff,limit)*sum(dz[:il+1])
                    i += 1
            top = botm
    else : ## radial and Xsection cases
        Zblock = ones((ny+1,1,nx));#print yvect
        for i in range(len(yvect)):
            Zblock[i] = yvect[ny-i];#print i,yvect[ny-i]
    #print Zblock
    return Zblock
    
def getTopBotm(core,modName,line,intp,im,refer,mat):#,optionT # EV 19/02/20
    if intp==1: # EV 11/02/20
        parms = core.dicinterp[modName][line][im] # EV 19/02/20
        z, mess = zone2interp(core,modName,line,im,parms,refer,iper=0) # EV 19/02/20
    elif intp==3 :# EV 11/02/20
        if  core.addin.mesh==None: # OA 19/4/20
            z = zone2grid(core,modName,line,im) # OA 18/7/19 modif type mesh for min"p
        else : 
            if modName == 'Opgeo': z = zone2mesh(core,modName,line,im,loc='nodes') #OA 17/2/20 replace mgroup
            else : z = zone2mesh(core,modName,line,im,loc='elements')
    elif intp==4 : # EV 11/02/20
        z = zone2array(core,modName,line,im) # EV 20/02/20
        if z.size == 0 : #EV 01/04/20
            z = zone2grid(core,modName,line,im)
            core.dictype[modName][line][im]='one_value'
    if shape(refer)==shape(z) : return minimum(z,refer-1)
    else : return z
    
def getXYvects(core):
    g = core.addin.getFullGrid()
    xvect=concatenate([array(g['x0'],ndmin=1),g['x0']+cumsum(g['dx'])],axis=0)
    yvect=concatenate([array(g['y0'],ndmin=1),g['y0']+cumsum(g['dy'])],axis=0)
    #print 'geom 180 grid',g, xvect
    return g['nx'],g['ny'],xvect,yvect

def getXYmeshCenters(core,plane,section):
    '''returns the center of each cell for a given plane and orientation'''
    a,b,xv0,yv0 = getXYvects(core);#print 'geom l 141',plane,layer
    xv1, yv1 = (xv0[1:]+xv0[:-1])/2, (yv0[1:]+yv0[:-1])/2
    nlay = getNlayers(core)
    zb = core.Zblock
    if core.addin.getDim() not in ['Radial','Xsection']: # 2 or 3D case
        if plane=='Z': return meshgrid(xv1,yv1)
        elif plane == 'Y':
            return xv1*ones((nlay,1)),zb[:,section,:]
        elif plane == 'X':
            return yv1*ones((nlay,1)),zb[:,:,section]
    else : # radial or cross section
        return meshgrid(xv1,yv1)

def getXYmeshSides(core,plane,section):
    '''returns the sides of the cells for a given plane and orientation'''
    a,b,xv1,yv1 = getXYvects(core);#print 'geom l 141',plane,layer
    nlay = getNlayers(core)
    zb = core.Zblock
    if core.addin.getDim() not in ['Radial','Xsection']: # 2 or 3D case
        if plane=='Z': return meshgrid(xv1,yv1)
        elif plane == 'Y':
            zb1 = zb[:,section,:]
            return xv1*ones((nlay+1,1)), c_[zb1[:,:1],(zb1[:,:-1]+zb1[:,1:])/2, zb1[:,-1:]]
        elif plane == 'X':
            zb1 = zb[:,:,section]
            return yv1*ones((nlay+1,1)), c_[zb1[:,:1],(zb1[:,:-1]+zb1[:,1:])/2, zb1[:,-1:]]
    else : # radial or cross section
        return meshgrid(xv1,yv1)
        
def getMesh3D(core):
    '''returns the coordinates of the cell sides in 3D, only for structured'''
    x1,y1=getXYmeshSides(core,'Z',0);ny,nx=shape(x1)
    z1=core.Zblock;nz,a,b=shape(z1)
    z2=concatenate([z1[:,:,:1],z1],axis=2)
    z2=concatenate([z2[:,:1,:],z2],axis=1)
    x2=x1*ones((nz,ny,nx))
    y2=y1*ones((nz,ny,nx))
    return x2,y2,z2
    
def getMesh3Dcenters(core):
    '''returns the coordinates of the cell centers in 3D, only for structured'''
    x1,y1=getXYmeshCenters(core,'Z',0)
    z1=core.Zblock
    z2=z1[1:,:,:]/2+z1[:-1,:,:]/2
    x2=x1*ones(shape(z2))
    y2=x2*ones(shape(z2))
    return x2,y2,z2
    
def block(core,modName,line,intp=False,opt=None,iper=0):
    #print("block",core.addin.mesh)
    if core.addin.mesh == None: # OA 29/2/20 added mstType rect
        return blockRegular(core,modName,line,intp,opt,iper)
    else : #mfUnstruct or modName[:5]=='Opgeo':
        return blockUnstruct(core,modName,line,intp,opt,iper)
    
def block1(core,modName,line,intp=False,opt=None,nvar=2):
    #a special block for n variables (2 for drn,ghb or chd, 3 for riv)
    iper=0
    if core.addin.mesh == None: # OA 29/2/20 added mstType rect
        return blockRegular(core,modName,line,intp,opt,iper)
    else : #mfUnstruct or modName[:5]=='Opgeo':
        blk = ones(nvar,(getNlayers(core),core.addin.mesh.getNumber('elements')))
        return blockUnstruct(core,modName,line,intp,opt,iper)
        
def blockUnstruct(core,modName,line,intp,opt,iper):
    '''returns data for a 3D unstructured block'''
#    if 1 in intp : # OA 17/2/20
#        parms = core.dicinterp[modName][line][0]
#        m, mess = zone2interp(core,modName,line,0,parms,iper=iper);#print 'geom 282',amax(m) # EV 19/02/20
#        m = array(m,ndmin=2) # OA 1/5/20 this and three below added
#        for im in range(1,getNmedia(core)):
#            m = r_[m,array(zone2interp(core,modName,line,0,parms,iper=iper),ndmin=2)]
#        return m
#    else: # changed from 0 to []
#        m = array(zone2mesh(core,modName,line,0),ndmin=2) # just one layer
#        for im in range(1,getNmedia(core)):
#            m = r_[m,array(zone2mesh(core,modName,line,im),ndmin=2)]
#        return m
    m0 = ones((getNlayers(core),core.addin.mesh.getNumber('elements')))
    nmedia = getNmedia(core)
    if type(intp) != type([5,6]) : intp=[0]*nmedia # OA added 28/10/20
    lilay = getNlayersPerMedia(core)
    lay = 0
    for im in range(nmedia): # 3D case, includes 2D
        if intp[im]==1 :
            parms = core.dicinterp[modName][line][im] # EV 19/02/20
            a,mess = zone2interp(core,modName,line,im,parms,iper=iper) # EV 19/02/20
        elif intp[im]==3 :
            a = zone2mesh(core,modName,line,im)
        elif intp[im]==4 :
            #try : 
            a = zone2array(core,modName,line,im) 
            if a.size == 0 : 
                a = zone2mesh(core,modName,line,im)
                core.dictype[modName][line][im]='one_value'
        else : # OA added 28/10/20
            a = zone2mesh(core,modName,line,im)            
        for il in range(int(lilay[im])): # several layers can exist in each media
            m0[lay]=a
            lay +=1
    return m0
        
def zone2mesh(core,modName,line,media=0,iper=0,loc='elements',val='value'):
    """return a vector of values for one property over a mesh, values are given
    through zones
    can return values on elts or nodes according to loc value
    if val = 'value' returns the value, if val ='nb', return the zone number"""
    vbase = 0
    if val =='value': 
        lval = core.dicval[modName][line]
        if media<len(lval): vbase=float(lval[media])
        else : vbase =float(lval[0])
    mesh = core.addin.mesh
    nbc = mesh.getNumber(loc);
    value = ones(nbc)*vbase
    dicz = core.diczone[modName].dic#;print 'geom 292 shapval',shape(value),dicz
    if line not in list(dicz.keys()):
        return value
    else: # OA 4/8/19
        if len(dicz[line]['name'])==0: return value
    dicz = dicz[line]
    zval,onev = [],True # zval will contain all the values in the zone list
    for iz,v in enumerate(dicz['value']): # OA 20/2/20
        if '$' in v: # here we take only the 3 value (for riv drn, ghb)
            onev = False # OA 20/2/20
            v1 = v.split('$')[2]
            try : zval.append(int(v1))
            except ValueError : pass
        elif '\n' in v: # a transient zone added OA 20/2/20
            zval.append(float(core.ttable[line][iper,iz]))
        else : 
            zval.append(float(v))
    #if onev: zval.insert(0,vbase) # not for multiple value where the domain is a zone
    pindx = zeros(nbc) # this will be the index, 0 for the background, then poly number
    for i in range(len(dicz['name'])):
        if (dicz['name'][i]=='domain'): continue
        idx,zv0 = zmesh(core,dicz,media,i)
        if type(zv0) != type(0) : value[idx]=zv0[idx]
        else : value[idx]=zval[i]
        try : len(idx) 
        except TypeError : continue # the zone media is not the right one
        pindx[idx] = i+1
    if val == 'nb': return pindx
    else : return value
    
def zmesh(core,dicz,media,iz):
    '''returns the index of the cells that are under (line) or in (poly) a zone
    it can also provide the z value if the zone is the variable polygon'''
    mesh = core.addin.mesh
    xc, yc = mesh.elcenters[:,0],mesh.elcenters[:,1]
    poly = dicz['coords'][iz];#print poly
    if len(poly[0])==3 : x,y,z = list(zip(*poly)) # OA 2/5/20
    else : x,y = list(zip(*poly))
    if len(x)>1: d = sqrt((x[1]-x[0])**2+(y[1]-y[0])**2)
    llcoefs = lcoefsFromPoly(poly)
#    zmedia = dicz['media'][i] # a media or a list of media for the zone OA removed 24/7/20
#    if type(zmedia)!=type([5]): zmedia=[zmedia]
#    zmedia = [int(a) for a in zmedia]
#    if media not in zmedia: return None # the zone is not in the correct media
    if len(poly)==1: # one point
        dst = sqrt((poly[0][0]-xc)**2+(poly[0][1]-yc)**2)
        idx = amin(dst)==dst;#where(amin(dst)==dst)[0] # OA 19/4/20
        zval = 0
    elif (abs(x[0]-x[-1])<d/20) & (abs(y[0]-y[-1])<d/20): #a closed polygon
        idx = where(pointsInPoly(xc,yc,poly,llcoefs));
        zval = 0 #OA 2/5/20 add zval
    else : # a line
        id0,zval = cellsUnderPoly(core,dicz,media,iz) # OA 2/5/20 added zval
        idx = id0>0
    return idx,zval
        
def blockRegular(core,modName,line,intp,opt,iper):
    """creates a block for one variable which has a size given by the x,y,z 
    given by addin
    intp : flag for interpolation can be one bool (for all layers) or a list of bool"""
    nx,ny,xvect,yvect = getXYvects(core)
    nmedia = getNmedia(core)
    nlayers = getNlayers(core)
    lilay = getNlayersPerMedia(core)
    g = core.addin.getFullGrid()
    if type(intp)==type(0): intp=[intp] #OA 16/02/20 # EV 05/05/20
    if len(intp)<nmedia: intp=intp*nmedia ## to extend intp at the number of media
    #### interpolation and array option for 2D and 3D model #EV 07/02/20
    if core.addin.getDim() in ['2D','3D']:
        m0 = ones((nlayers,ny,nx))
        lay = 0
        for im in range(nmedia): # 3D case, includes 2D
            #print('lin',line, 'im',im, 'intp', intp[im])
            if intp[im]==1 :
                parms = core.dicinterp[modName][line][im] # EV 19/02/20
                a,mess = zone2interp(core,modName,line,im,parms,iper=iper) # EV 19/02/20
            elif intp[im]==3 :
                a = zone2grid(core,modName,line,im,opt,iper)
            elif intp[im]==4 :
                #try : 
                a = zone2array(core,modName,line,im) # EV 20/02/20
                if a.size == 0 : #EV 01/04/20
                    a = zone2grid(core,modName,line,im,opt,iper)
                    core.dictype[modName][line][im]='one_value'
            for il in range(int(lilay[im])): # several layers can exist in each media
                m0[lay]=a
                lay +=1
    #### interpolation and array option for Xsection & radial model #EV 07/02/20
    else : 
        for im in range(nmedia): 
            if intp[im]==1 :
                parms = core.dicinterp[modName][line][im] # EV 19/02/20
                a,mess = zone2interp(core,modName,line,im,parms,iper=iper) # EV 19/02/20
            elif intp[im]==3 :
                a = zone2grid(core,modName,line,im,opt,iper)
            elif intp[im]==4 :
                a = zone2array()
                if a.size == 0 : #EV 01/04/20
                    a = zone2grid(core,modName,line,im,opt,iper)
                    core.dictype[modName][line][im]='one_value'
        #a = zone2grid(core,modName,line,0,opt,iper)
        m0 = reshape(a,(ny,1,nx))
        m0 = m0[-1::-1]
    linesRadial = ['bcf.4','bcf.6','bcf.7','bcf.8','lpf.8','lpf.10','lpf.11',
        'btn.11','rct.2a','rct.2b','rch.2']
    dx = array(g['dx'])
    if (core.addin.getDim()=='Radial') & (line in linesRadial):
        for i in range(ny): m0[i] = m0[i]*(cumsum(dx)-dx/2.)*6.28
    if line=='lpf.9' and core.getValueFromName('Modflow','CHANI')==0:  # ratio Kh/Kv shall not be scaled
        for i in range(ny): m0[i] = m0[i]*(cumsum(dx)-dx/2.)*6.28
    return m0

def zone2grid(core,modName,line,media,opt=None,iper=0):
    """ Put zone information on the correct grid.
    If the grid is 3D then looks for information of layers
    zone keys 'number','name','coords','media','value','type'.
    Shall work with regular and vairable grids, but does not provide filled zone 
    if the zone contains z values (variable polygon).
    If line in observation (obs.1), opt is a list of zone to consider and 
    the nb of zone instead of value is used.
    """
    lval = core.dicval[modName][line]
    if media<len(lval): vbase=float(lval[media])
    else : vbase =float(lval[0])
    if line != 'obs.1' and opt in['BC','zon']: vbase=0 #EV 26/02/20
    nx,ny,xvect,yvect = getXYvects(core)
    m0=ones((ny,nx))*vbase
    #print 'geom 322',line,core.dictype[modName][line],#core.diczone[modName].dic
    if core.dictype[modName][line][0]=='array': # type of array, dont calculate
        arr = core.dicarray[modName][line];#print 'geom zone2g type arr',line
        if len(shape(arr))==3: arr = arr[media]
        return arr
    elif line in list(core.diczone[modName].dic.keys()):
        diczone = core.diczone[modName].dic[line]
    else : 
        return array(m0) # returns an array of vbase
    lz = len(diczone['name'])
    #print 'geo z2g',line,core.ttable[line]
    
    for i in range(lz):  # loop on zones
        if line != 'obs.1' and diczone['value'][i]=='': continue
        xy = diczone['coords'][i]#;core.gui.onMessage(str(xy))
        zmedia = diczone['media'][i] # a media or a list of media for the zone
        if type(zmedia)!=type([5]): zmedia=[zmedia]
        zmedia = [int(a) for a in zmedia]
        if int(media) not in zmedia: continue # the zone is not in the correct media
        #if line in list(core.ttable.keys()): zv0=float(core.ttable[line][iper,i])
        if line in list(core.ttable.keys()): 
            if line == 'obs.1' : #EV 26/02/20
                if diczone['name'][i] in opt: 
                    zv0= int(diczone['name'].index(diczone['name'][i]))+1  #; print('zv0',diczone['name'][i],zv0)
                else : continue
            else :
                if line in ['drn.1','riv.1','ghb.1']: #EV 26/11/20
                    zv0=float(core.ttable[line][iper,i].split()[opt])
                else : zv0=float(core.ttable[line][iper,i])
        if line != 'obs.1' and opt=='zon': zv0=i+1 # OA added 8/5/19 #EV 26/02/20
        #if type(zv0)!=type(5.) and '$' in diczone['value'][i]: zv0 = float(diczone['value'][i].split('$')[2])# added 17/04/2017
        #else : zv0 = float(diczone['value'][i])
        if len(xy)==1:  # case point
            x,y=list(zip(*xy))
            ix,iy = minDiff(x[0],xvect),minDiff(y[0],yvect)
            m0[iy,ix] = zv0
            
        else: # other shapes
            ndim=len(xy[0]);
            if ndim==3: x,y,z = list(zip(*xy)) # there are z values (variable polygon)
            else :
                x,y = list(zip(*xy))
                z=[zv0]*len(xy)
            nxp,nyp,nzp=zone2index(core,x,y,z);
            put(m0,nyp*nx+nxp,nzp)
            if ndim==3: continue # a zone with z value is not filled!!
            ind = fillZone(nx,ny,nxp,nyp,nzp)
            putmask(m0, ind==1, [zv0])
    #if line =='drn.1':print('m0',m0)
    return array(m0)

def fillZone(nx,ny,nxp,nyp,nzp):
    ind = zeros((ny,nx))
    if type(nxp)!=type(array([5])): #just one index OA Three lines added 31/07
        ind[nyp,nxp]= 1
        return ind
    # fill zone with vertical lines
    js = argsort(nxp);
    nxs, nys = take(nxp,js), take(nyp,js)            
    lls = len(nxs)
    ind1 = zeros((ny,nx))
    mn, mx = int(nys[0]),int(nys[0])
    for j in range(1,lls):
        if nxs[j]!=nxs[j-1]:
            ind1[mn:mx+1,int(nxs[j-1])] = [1]*(mx-mn+1)
            mn,mx = int(nys[j]),int(nys[j])                    
        mn,mx = int(min(mn,nys[j])),int(max(mx,nys[j]))
    ind1[mn:mx+1,int(nxs[lls-1])] = [1]*(mx-mn+1)
    # fill zone with horizontal lines
    js = argsort(nyp)
    nxs, nys = take(nxp,js),take(nyp,js)
    ind2 = zeros((ny,nx))
    mn, mx = int(nxs[0]),int(nxs[0])
    for j in range(1,lls):
        if nys[j]!=nys[j-1]:
            ind2[int(nys[j-1]),mn:mx+1] = [1]*(mx-mn+1)
            mn,mx = int(nxs[j]),int(nxs[j])                    
        mn,mx = int(min(mn,nxs[j])),int(max(mx,nxs[j]))                
    ind2[int(nys[lls-1]),mn:mx+1] = [1]*(mx-mn+1)
    ind = ind1*ind2  ## both indices must equal 1 to fill the zone
    return ind

def zone2index(core,x,y,z,opt=None):
    """ returns indices of cells below a poly zone"""
    nx,ny,xvect,yvect = getXYvects(core)
    nxp,nyp,nzp=array([minDiff(x[0],xvect)]),array([minDiff(y[0],yvect)]),array(z)
    #print 'zoneindex',nxp,nyp,nzp
    if len(x)==1:
        if opt==None: return nxp,nyp,nzp
        elif opt=='angle': return nxp,nyp,nzp,[0],[0]
    nxp=array([],dtype='int');nyp=nxp*1;
    nzp=array([],dtype='float');nsn=nzp*1;ncs=nzp*1
    for j in range(1,len(x)):
        ix1=minDiff(x[j-1],xvect);ix2=minDiff(x[j],xvect);
        iy1=minDiff(y[j-1],yvect);iy2=minDiff(y[j],yvect);#print y[0],yvect
        if ix1==ix2 and iy1==iy2: continue
        sensx,sensy=sign(ix2-ix1),sign(iy2-iy1);
        lx,ly=abs(ix2-ix1),abs(iy2-iy1)
        ll=max(lx,ly)
        dz=z[j]-z[j-1];#print ix1,ix2,iy1,iy2,sensx,xvect
        if lx>=ly:
            ixp=arange(ix1,ix2+sensx,sensx);
            iyp=ixp*0;xv2=xvect[ixp]
            yv2=y[j-1]+(xv2-xv2[0])*(yvect[iy2]-yvect[iy1])/(xv2[-1]-xv2[0]);
            for k in range(ll+1): iyp[k]=minDiff(yv2[k],yvect)
            zp=z[j-1]+dz*(xv2-x[j-1])/(x[j]-x[j-1]);#print ixp,xv2,iyp,yv2,zp
        else:
            iyp=arange(iy1,iy2+sensy,sensy);
            ixp=iyp*0;yv2=yvect[iyp];
            xv2=x[j-1]+(yv2-yv2[0])*(xvect[ix2]-xvect[ix1])/(yv2[-1]-yv2[0]);
            for k in range(ll+1): ixp[k]=minDiff(xv2[k],xvect)
            zp=z[j-1]+dz*(yv2-y[j-1])/(y[j]-y[j-1]);
        nxp=concatenate([nxp,clip(ixp,0,nx-1)],axis=0)           
        nyp=concatenate([nyp,clip(iyp,0,ny-1)],axis=0)
        nzp=concatenate([nzp,zp],axis=0)     
        dx,dy = x[j]-x[j-1],y[j]-y[j-1]
        sn0,cs0 = dy/sqrt(dx**2+dy**2),dx/sqrt(dx**2+dy**2)
        sn,cs = ones((ll+1))*sn0,ones((ll+1))*cs0
        nsn,ncs = concatenate([nsn,sn],axis=0),concatenate([ncs,cs],axis=0)                   
    mix=nxp*1000+nyp;
    a,ind=unique(mix,return_index=True);ind=sort(ind)
    if len(nxp)<1:
        nxp,nyp,nzp=array([minDiff(x[0],xvect)]),array([minDiff(y[0],yvect)]),array(z)
        ind,nsn,ncs = 0,[0],[0] # OA added 11/6/19
    #print 'geom zoneind',nxp,nyp,nzp
    if opt==None:
        if type(ind)==type(5): return [int(nxp[ind])],[int(nyp[ind])],[nzp[ind]] # OA added 7/6/20
        else : return nxp[ind].astype(int),nyp[ind].astype(int),nzp[ind]
    elif opt=='angle':
        return nxp[ind].astype(int),nyp[ind].astype(int),nzp[ind],nsn[ind],ncs[ind]

def minDiff(x,xvect):
    d=x-xvect; d1=d[d>0.]; #print('d',d)
    if len(d1)==0: return 0
    else :
        a=where(d==amin(d1));return int(a[0])
        
def isclosed(core,x,y):
    nx,ny,xv,yv=getXYvects(core) #print('x0',x[0],'y0',y[0])
    ix,iy = minDiff(x[0],xv),minDiff(y[0],yv)
    dmin=min(xv[ix+1]-xv[ix],yv[iy+1]-yv[iy]) # modified 22/3/2017 for variable grid
    d0=sqrt((x[0]-x[-1])**2+(y[0]-y[-1])**2);
    if d0<dmin: return True
    else : return False     
    
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
            if l<=nbdy: dcout['bc'+str(l-1)] = c3[c3[:,4]==l,-2:]-1
            else : 
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
    
################################### mesh utilities ##########################
#############################################################################    

def planeFromPoints(poly):
    '''coordinates of  aplane from n points in a poly
    sur http://www.les-mathematiques.net/phorum/read.php?13,728203,728210
    renvoie les fact a,b,c,d tels que ax+by+cz+d = 0 (ou z=-(d+ax+by)/c) eq plan'''
    x,y,z = list(zip(*poly))
    xg,yg,zg = mean(x),mean(y),mean(z)
    x,y,z = r_[x,x[0]],r_[y,y[0]],r_[z,z[0]]
    a = sum((y[:-1]-y[1:])*(z[:-1]+z[1:]))
    b = sum((z[:-1]-z[1:])*(x[:-1]+x[1:]))
    c = sum((x[:-1]-x[1:])*(y[:-1]+y[1:]))
    d = -a*xg-b*yg-c*zg
    return (a,b,c,d)
    
def planesFromPointMat(xmat,ymat,zmat):
    x,y,z = c_[xmat,xmat[:,0]],c_[ymat,ymat[:,0]],c_[zmat,zmat[:,0]]   
    a = sum((y[:,:-1]-y[:,1:])*(z[:,:-1]+z[:,1:]),axis=1)
    b = sum((z[:,:-1]-z[:,1:])*(x[:,:-1]+x[:,1:]),axis=1)
    c = sum((x[:,:-1]-x[:,1:])*(y[:,:-1]+y[:,1:]),axis=1)
    d=0 # not required
    return a,b,c,d

def isPolyTrigo(lpoints): 
    '''verifies if a poly is trigono direction or not, lpoints is a list of pts
    does NOT work for complex polygons, but here it is only for the domain
    poly must be closed'''
    x,y = list(zip(*lpoints))
    x,y = array(x),array(y)
    trigo = sum((x[1:]-x[:-1])*(y[1:]+y[:-1])) # if that sum is >0 it is clockw
    return trigo<0
    
def lcoefsFromPoly(poly):
    '''finds the list of coefficients for the equation of each line 
    segment in a polygon'''
    if len(poly[0])==3 : x,y,z = list(zip(*poly)) # OA 2/5/20
    else : x,y = list(zip(*poly))
    n = len(x)
    lcoefs,a=zeros((2,n-1)),zeros((2,2))
    for i in range(n-1): # lcoefs are the coefficient of the lines ax+by=1
        if x[i+1]==x[i]  and y[i+1]==y[i]: continue
        elif y[i+1]==y[i]: lcoefs[:,i]= [0,1/(y[i]+1e-18)] # 1e-18 not to have 0
        elif x[i+1]==x[i]: lcoefs[:,i]= [1/(x[i]+1e-18),0] # 1e-18 not to have 0
        else :
            a[:,0]=x[i:i+2];a[:,1]=y[i:i+2];#print i,a
            lcoefs[:,i] = solve(a,ones((2,1))[:,0])
    return lcoefs

def pointsInPoly(ptx,pty,poly,lcoefs): # ray algorithm
    '''finds if a point is in a polygon for a list of points of coords ptx pty
    we use an horizontal line that start from the ptx,pty point
    lcoefs [2,n] for n pts in poly are the eq of ax+by=1 with a,b=lcoefs[:,i]
    polygon is in a trigo direction and must be closed'''
    ptx,pty =array(ptx),array(pty) # OA 20/2/20
    x,y = list(zip(*poly))
    n,count = len(x),ptx*0
    for i in range(n-1):
        ymin,ymax = min(y[i:i+2]),max(y[i:i+2])
        pos = lcoefs[0,i]*ptx+lcoefs[1,i]*pty-1
        #if pos==0 : return True
        count += (pty>=ymin)*(pty<=ymax)*( sign(pos)==sign(lcoefs[0,i]))
    try :
        len(count)
        return list((np.mod(count,2)==1)*1)# an odd number gives point in the poly   
    except TypeError : 
        return [(np.mod(count,2)==1)*1]

def dstPointPoly(lpts,poly):
    '''finds the shortest distance from each point in a list (lpts) to a polygon (poly)'''   
    xpo,ypo = list(zip(*poly));
    dst,i = [],0
    for x,y in lpts[1:]:
        dpoly = sqrt((x-array(xpo))**2+(y-array(ypo))**2);#print i
        idx = where(dpoly==amin(dpoly))[0][0]
        xym = array([[xpo[idx-1],ypo[idx-1]],[xpo[idx],ypo[idx]]])
        if abs(xym[0,0]-xym[1,0])<1e-9: a,b,c = 1,0,xym[0,0] # vert line
        elif abs(xym[0,1]-xym[1,1])<1e-9: a,b,c = 0,1,xym[0,1] # hor line
        else : a,b = solve(xym,ones((2,1))[:,0])
        yp = (b+a*(-b*x+a*y))/(a**2+b**2);xp=(1-b*yp)/a
        if ((xp>min(xym[:,0]))&(xp<max(xym[:,0]))&(yp>min(xym[:,1]))&(yp<max(xym[:,1])))|(idx==len(xpo)-1):
            dst.append(abs(a*x+b*y-1)/sqrt(a**2+b**2))
        else :
            xym = array([[xpo[idx],ypo[idx]],[xpo[idx+1],ypo[idx+1]]])
            a,b = solve(xym,ones((2,1))[:,0])
            dst.append(abs(a*x+b*y-1)/sqrt(a**2+b**2))
        i += 1
    return dst
        
def cellsUnderPoly(core,dicz,media,iz):
    '''finds the cells that are under each line of a poly. 
    It uses an array of coordinates of the cell vertices
    Then calculates the position of the cell points
    relative to a line (can be one line of a  poly) and if there are both
    positive and negative values, the cell is under the line
    the input is two arrays for x and y vertices positions and the line pos'''
    mesh = core.addin.mesh
    xc,yc,idcell = mesh.elcenters[:,0],mesh.elcenters[:,1],mesh.idc
    elxa, elya = array(mesh.elxa),array(mesh.elya)
    poly = dicz['coords'][iz]
    if len(poly)==1: #OA 18/12/20
        x,y = poly[0];dst=(x-xc)**2+(y-yc)**2
        indx=where(dst==amin(dst))[0]
        return indx,dicz['value'][iz]
    lcoefs=lcoefsFromPoly(poly)
    #l0 = zptsIndices(core,dicz)[iz] # OA 24/7/20
    indx,zval = zeros(len(idcell)),zeros(len(idcell))
    for i in range(len(poly)-1):
        if len(poly[0])==3 : x,y,z = list(zip(*poly[i:i+2])) # OA 2/5/20
        else : x,y = list(zip(*poly[i:i+2]));z=0
        pos = sign(lcoefs[0,i]*elxa+lcoefs[1,i]*elya-1) # posit. vs the line
        dpos = [abs(amax(pos[id[0]:id[1]])-amin(pos[id[0]:id[1]])) for id in idcell]
        dpos = array(dpos)
        dst0,dst1 = sqrt((x[0]-xc)**2+(y[0]-yc)**2),sqrt((x[1]-xc)**2+(y[1]-yc)**2) # OA 24/7/20
        idx0 = amin(dst0)==dst0;idx1 = amin(dst1)==dst1; # index of ends of line OA 24/7/20
        x0,x1 = min(xc[idx0],xc[idx1]),max(xc[idx0],xc[idx1])  # OA 28/10/20
        y0,y1 = min(yc[idx0],yc[idx1]),max(yc[idx0],yc[idx1])  # OA 28/10/20
        idx = (dpos>0)*(xc<=x1)*(xc>=x0)*(yc<=y1)*(yc>=y0)*1
        ipts = where(idx==1)[0]  # OA 2/5/20 this and two lines below added
        indx += idx
        dst = sqrt((xc[ipts]-x[0])**2+(yc[ipts]-y[0])**2)/sqrt((x[1]-x[0])**2+(y[1]-y[0])**2)
        if z == 0 : zval = 0
        else : zval[ipts] = z[0]+(z[1]-z[0])*dst
    return clip(indx,0,1),zval

def zptsIndices(core,dicz):
    '''finds the indices of the points in the zones'''
    mgroup = core.dicaddin['Model']['group']
#    if core.getValueFromName('Modflow','MshType')<1: 
#        nx,ny,xv,yv=getXYvects(core)
#        lindx = []
#        for iz in range(len(dicz['name'])):
#            l0 = []
#            for coo in dicz['coords'][iz]:
#                ix = where(argsort(r_[xv,coo[0]])==nx+1)[0][0]-1
#                iy = where(argsort(r_[yv,coo[1]])==ny+1)[0][0]-1
#                l0.append((ix,iy))
#            lindx.append(l0)
#    else : # fully unstruct
    mesh = core.addin.mesh
    xc,yc = mesh.getCenters()
    lindx = []
    for iz in range(len(dicz['name'])):
        l0 = []
        for coo in dicz['coords'][iz]:
            dst=sqrt((coo[0]-xc)**2+(coo[1]-yc)**2)
            iclosept = argsort(dst)[:3]
            for ip in iclosept:
                xp,yp = mesh.elx[ip],mesh.ely[ip] # OA 19/4/20 modif
                poly=list(zip(xp,yp))
                if isclosed(core,xp,yp)==False: poly.append(poly[0]) #OA 19/4/20 added
                lcoefs=lcoefsFromPoly(poly)
                if pointsInPoly(coo[0],coo[1],poly,lcoefs)[0]:
                    l0.append(ip)
        lindx.append(l0)
    return lindx
    
def findSideNodes(core,nodes):
    '''retunrs index of nodes that are on the sides of the domain
    nodes is the array having the x,y coordinates (z are always 0 in nodes)'''
    grd = core.addin.getFullGrid()
    eps = grd['epsilon']
    x,y = nodes[:,0],nodes[:,1]
    indx = x*0
    indx[(x>grd['x0']-eps)&(x<grd['x0']+eps)]=1
    indx[(x>grd['x1']-eps)&(x<grd['x1']+eps)]=1
    indx[(y>grd['y0']-eps)&(y<grd['y0']+eps)]=1
    indx[(y>grd['y1']-eps)&(x<grd['y1']+eps)]=1
    return indx
    
def writeVTKstruct(core,data):
    a,b,xv1,yv1 = getXYvects(core);#print 'geom l 141',plane,layer
    nlay = getNlayers(core)
    ym,zm,xm = meshgrid(yv1,ones(nlay+1),xv1)
    zb = core.Zblock
    zb1 = concatenate([zb[:,:,:1],(zb[:,:,:-1]+zb[:,:,1:])/2, zb[:,:,-1:]],axis=2)
    zm = concatenate([zb1[:,:1,:],(zb1[:,:-1,:]+zb1[:,1:,:])/2, zb1[:,-1:,:]],axis=1)
    #print 'geom 210',shape(xm),shape(zb),shape(zm),shape(data)
    s= '# vtk DataFile Version 3.0 \n scalar \nASCII \n'
    nz,ny,nx = shape(xm)
    x,y,z,data = ravel(xm),ravel(ym),ravel(zm),ravel(data)
    npts,ndata = len(x),len(data);#print 'foem', npts,ndata
    s += 'DATASET STRUCTURED_GRID \nDIMENSIONS '+str(nx)+' '+str(ny)+' '+str(nz)+'\n'
    s += 'POINTS '+str(npts)+' float \n'
    s += ' '.join([str(x[i])+' '+str(y[i])+' '+str(z[i])+'\n' for i in range(npts)])
    s += 'CELL_DATA '+str(ndata)+'\nSCALARS cellData float\nLOOKUP_TABLE default\n'
    s += ' '.join([str(data[i])+'\n' for i in range(ndata)])
    return s

###################### INTERPOLATION ####################

def zone2interp(core,modName,line,media,option,refer=None,iper=None): # EV 19/02/20 #EV 21/02/20
    """ interpolation case from points or zones to matrix
    option[0] : interpolation method 
        0: kriging
        1: Inverse distance
        2: Triangulation
    option[1-8] : option for interpolation 
        kringing: 1.Automatic parameter(True,False), 2.Log Values(True,False), 3.Plot Variogram(True,False),
                  4.nlags(6), 5.Variogram model(power,gaussian,Spherical,Exponential), 
                  6.Sill(Text), 7.Range(Text), 8.Nugget(Text), 
        Inverse distance : Power(1,2,3), Radius(text)
        Triangulation : """
    #### Interpolation method and base
    intpMtd = int(option[0]) ## interpolation method
    vbase=float(core.dicval[modName][line][media])
    #variable=None EV 20/02/20
    logOpt = 0  ## option of log transformation for interpolation initializing to 0 (false) 
    
    #### create the vector of points on which interpolation will be done
    mess = None
    if modName[:5] == 'Opgeo' or core.addin.mesh != None: # OA 1/5/20
        mesh = core.addin.mesh  # OA 1/5/20
        xc, yc = mesh.elcenters[:,0],mesh.elcenters[:,1]
    else :
        nx,ny,xv,yv = getXYvects(core);nz=0
        xv,yv = (xv[1:]+xv[:-1])/2.,(yv[1:]+yv[:-1])/2.;#print 'xv',xv,yv
        xm,ym = meshgrid(xv,yv)
        xc, yc = ravel(xm),ravel(ym)
        if type(refer) == type(zeros(5)) : refer = ravel(refer) # OA 21/8/20 
    m = zeros(len(xc)) + float(vbase)
    if line in list(core.diczone[modName].dic.keys()):
        diczone = core.diczone[modName].dic[line]
        coords = diczone['coords']
        nz = len(coords)
    else : 
        if modName[:5] != 'Opgeo':
            m = reshape(m,(ny,nx))
        return m, mess # strange!!
    
    #### creates the list of point values
    xpt,ypt,zpt=[],[],[];
    for iz in range(nz):  # loop on zones to get coordinates
        lmed = int(diczone['media'][iz])#*1
        if type(lmed) != type([5,6]): lmed = [lmed]# not a list
        if media not in lmed : 
            continue
        xy = coords[iz]
        if len(xy[0])==2: #normal zone
            x,y = list(zip(*xy))
            if line in list(core.ttable.keys()): #EV 21/02/20
                z=float(core.ttable[line][iper,iz]) 
            #if variable != None: z = float(diczone['value'][iz].split('$')[1].split('\n')[variable]) EV 20/02/20
            else : z = float(diczone['value'][iz])
            z = [z]*len(x);#print i,z
        elif len(xy[0])==3: # variable polygon
            x,y,z = list(zip(*xy))
        #print(diczone['name'][iz], z)
        xpt.extend(list(x));ypt.extend(list(y));zpt.extend(list(z))            
    xpt,ypt,zpt = array(xpt),array(ypt),array(zpt);#print 'geom, xpt',xpt,ypt,zpt
    if len(xpt)==0: 
        if modName[:5] != 'Opgeo': m = reshape(m,(ny,nx))
        return m, mess
    
    #### makes the difference if refer is not none to avoid negative thickness (dis.6, dis.7)
    #print shape(xc),shape(yc)
    if type(refer) == type(zeros(5)): # OA 21/8/20
        for i in range(len(xpt)):
            d = sqrt((xpt[i]-xc)**2+(ypt[i]-yc)**2)
            ix = where(d==amin(d))[0];#print ix
            zpt[i] = refer[ix[0]]-zpt[i]
            
    #### interpolation method
    if intpMtd ==1: ## inverse distance
        pw = float(option[1])
        pol=polyfit2d(xpt,ypt,zpt,order=1)
        z0=polyval2d(xpt,ypt,pol);#print 'z0',z0
        dz=zpt-z0;#print 'dz',xpt,ypt,zpt,z0,dz
        m0=polyval2d(xc,yc,pol);#print 'm0',m0
        m1 = invDistance(xpt,ypt,dz,xc,yc,power=pw);#print 'm1',m1
        m2 = m0+m1 

    elif intpMtd ==0 : ## Kriging
        vparm = int(option[1])
        vtype =str(option[5]) ## variogram model
        if vparm == 1 : vparm = None ## automatic parameter calculation
        else :
            sill = float(option[6]) ## sill
            rg = float(option[7]) ## range
            ng = float(option[8]) ## nugget
            vparm =[sill, rg, ng] ## variogram parameter 
        nlags = int(option[4]) 
        logOpt = int(option[2]) ## log option
        if logOpt==1 : zpt=log10(zpt)
        if len(xpt)<5: 
            m2 = m
        #else : m2,mess = krige(xpt,ypt,zpt,rg,xc,yc,vtype)
        else : m2,mess = krige(xpt,ypt,zpt,xc,yc,vtype,vparm,nlags)#,vplot)
        if logOpt==1 : m2 = 10**m2
        #print('Mess',mess)
        
    elif intpMtd ==2: ## Thissen polygons (Triangulation)
        vrange=float(option[1])
        listP = list(zip(xpt,ypt))
        xb,yb = [min(xc),max(xc)],[min(yc),max(yc)]
        polylist = getPolyList(listP,xb,yb);#print 'geom 830 thiess',polylist
        m2 = ones(len(xc))*vbase
        for i,poly in enumerate(polylist):
            if poly==None: continue
            x,y = list(zip(*poly));z=[0]*len(x);#plot(x,y)
            nxp,nyp,nzp=zone2index(core,x,y,z);
            ind = fillZone(nx,ny,nxp,nyp,nzp)
            putmask(m2, ind==1, [zpt[i]]);#print 'geom 611',xpt[i],ypt[i],zpt[i],poly

    if modName[:5] != 'Opgeo' and core.addin.mesh == None: # OA 1/5/20
        m2 = reshape(m2,(ny,nx))
        if intpMtd ==2:   #smoothing of thiessen polys  # OA 17/3/19 added int()   
            for n in range(int(vrange)): m2 = smoo2d(m2)

    if logOpt!=1: m2 = clip(m2,amin(zpt)*0.9,amax(zpt)*1.1);#print amin(amin(m2))

    return m2, mess
    
def smoo2d(m):
    m1=m*1
    m1[1:-1,1:-1]=m1[1:-1,:-2]/5+m1[1:-1,1:-1]/5+m1[1:-1,2:]/5+m1[:-2,1:-1]/5+m1[2:,1:-1]/5
    return m1

###################### ARRAY ####################

'''Import array from one media (im) EV 07/02/20'''
def onMessage1(core,txt):
    cfg = Config(core)
    try : 
        dialogs = cfg.dialogs
        dialogs.onMessage(core.gui,txt) 
    except AttributeError: 
        print(txt)

def zone2array(core,modName,line,im):
    fNameExt = core.dicarray[modName][line][im] #EV 05/05/20
    #fNameExt = file.split('/')[-1] 
    #fDir = file.replace(fNameExt,'')
    fDir = core.fileDir
    ext=fNameExt[-3:]
    arr=array([]) ; zdx,zdy,ysign=None,None,-1 # OA 26/7/20 set ysign default to -1 (modflow)
    txt1 = ('The file '+'"'+core.fileDir+fNameExt+'"'+' does not exist.'+'\n\n'+
                   'Default values or zones values will be used.'+'\n\n'
                   +'Please select an other file to import an array for the parameter '+
                   str(line)+' at media '+str(im)+'.')
    txt2 = ('The format of file '+'"'+core.fileDir+fNameExt+'"'+' is not suitable.'+'\n\n'+
                   'Default values or zones values will be used.'+'\n\n'
                   +'Please select an other file to import an array for the parameter '+
                   str(line)+' at media '+str(im)+'.')
    if ext == 'asc':
        try : 
            arr=core.importAscii(fDir,fNameExt)
            arr=arr.astype(np.float)
        except OSError : onMessage1(core,txt1) 
        except : onMessage1(core,txt2) 
    elif ext == 'var' :
        try : 
            ysign,zdx,zdy,arr=core.importGridVar(fDir,fNameExt) # OA 13/6/20 add ysign
            #zdy = zdy[::-1] # OA added 9/4/20
            arr=arr.astype(np.float)
        except OSError : onMessage1(core,txt1) 
        except : onMessage1(core,txt2) 
    else : #if ext == 'txt' or ext == 'dat' :
        print(fDir+os.sep+fNameExt)
        try : arr=loadtxt(fDir+os.sep+fNameExt)  # OA 6/6/20 changed file did not exist
        except OSError :onMessage(core,txt1) 
        except : onMessage1(core,txt2) 
    #print(type(arr),shape(arr),arr[:1])
    if arr.size != 0:
        grd = core.addin.getFullGrid()
        intp = False # OA 3/4/20 this l an dl. below
        #if line in ['lpf.8']: intp=True #'dis.6','dis.7',
        if core.addin.mesh == None: xx,yy=getXYmeshCenters(core,'Z',0) # OA 24/7/20
        else : m = core.addin.mesh.getCenters();xx,yy = m[0],m[1] # OA 24/7/20
        if ysign==-1: # OA 13/6/20 added this and below
            arr = arr[-1::-1,:]*1
            if zdy != None: zdy = zdy[-1::-1]*1 #OA 26/7/20 added condition
        arr2 = linIntpFromGrid(grd,arr,xx,yy,intp,zdx,zdy) # removed [::-1]
        return arr2
        #else : return arr[::-1]
    else : 
        return arr        
       
    
