#import matplotlib
#matplotlib.use('WX')
#import matplotlib.backends.backend_wxagg
#import matplotlib
#matplotlib.use("agg")
import matplotlib.backends.backend_qt5agg
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as Toolbar
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.artist import Artist
from matplotlib.axes import Subplot
from matplotlib import rcParams,cm
import matplotlib.pylab as pl
from matplotlib.mlab import dist_point_to_segment
from matplotlib.patches import RegularPolygon,Polygon
from matplotlib.lines import Line2D
from matplotlib.collections import PolyCollection,LineCollection

#pour l'affichage d'une carte de fond
import matplotlib.image as Im
import matplotlib as mpl
import matplotlib.tri as mptri
import matplotlib.pyplot as plt
import time
#try: import aColor
#except ImportError: aColor = wx.Colour

from scipy.interpolate import interp1d
import numpy.ma as ma

from .geometry import *
from .qtDialogs import *
#from mayavi.mlab import *

def smoo(v):
    """ fonction qui permet de lisser une matrice par des moyennes mobiles"""
    [l,c] = shape(v)
    v1 = concatenate([v[:1,:],(v[:-2,:]+v[1:-1,:]+v[2:,:])/3,v[l-1:l,:]],axis=0)
    v2 = concatenate([v1[:,:1],(v1[:,:-2]+v[:,1:-1]+v1[:,2:])/3,v1[:,c-1:c]],axis=1)
    return v2

def aColor(b,g,r):
    '''to replace wx.Color'''
    return (b,g,r)

"""
classe representant un objet graphique qui sera ajoute aux axes (a l'objet Axes
du canvas en fait)
type : type de l'objet graphique, pour une zone le type est dans la liste de
    constante listeNomZone de Visualisation, sinon c'est une string, 'imshow'
    pour une image par exemple ou 'grid', 'contour'
object : l'objet graphique en lui-meme cree par python
visible : booleen qui indique si l'objet est visible (True) ou non (False)
color : couleur ou ensemble de couleurs de l'objet graphique
"""        
#/////////////////////////////////////////////////////////////////////
class qtVisualisation(FigureCanvasQTAgg):

    def __init__(self,gui):
    
        #super(wxVisualisation, self).__init__(parent=None, id=-1) 
        self.fig = pl.figure()
        #self.figurepanel = FigureCanvasWxAgg(self.panel, -1, self.fig)        
        FigureCanvasQTAgg.__init__(self, self.fig) #gui, -1, self.fig)
        policy = QSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.setSizePolicy(policy)
        self.xlim, self.ylim, self.dataon = (), (),False
        
        # c'est la GUI e tle modele
        self.gui,self.core = gui,gui.core
        # polygone d'interaction sur une zone (pour pouvoir la modifier)
        self.polyInteract = None
        
        # liste de variables, sert pour la GUI et le combobox sur les variables
        self.listeNomVar= []
        for mod in self.core.modelList:
            for k in list(gui.linesDic[mod].keys()):
                self.listeNomVar.extend(gui.linesDic[mod][k])
        self.curVar,self.curContour = None,'Charge' # variable courante selectionne       
        self.curMedia,self.curOri, self.curPlan,self.curGroupe=0,'Z',0,None
        # variable pour savoir si on est en cours de tracage d'une zone
        self.typeZone = -1
        
        # coordonnees et zone de la zone que l'on est en train de creer
        self.curZone = None  # objet graphique (ligne, point, rect..)
        self.x1, self.y1 = [], []
        self.tempZoneVal = []  # liste de values pour polyV
        self.calcE=0;self.calcT=0;self.calcR=0; # dit si calcule effectue ou non
        
        # dictionnaire qui est compose des variables de l'Aquifere
        # a chaque variable est associe une liste de zones
        self.listeZone, self.listeZoneText, self.listeZmedia = {},{},{}
        for i in range(len(self.listeNomVar)):
            #print self.listeNomVar[i]
            self.listeZone[self.listeNomVar[i]] = []
            self.listeZoneText[self.listeNomVar[i]] = []
            self.listeZmedia[self.listeNomVar[i]] = []

        # toolbar de la visu, de type NavigationToolbar2Wx
        self.toolbar = Toolbar(self,gui)
        self.toolbar.setFixedWidth(350)
        # ajout du subplot a la figure
        self.cnv = self.fig.add_axes([.05,.05,.9,.9]) #left,bottom, wide,height     
        #self.toolbar.update()    
        self.pos = self.mpl_connect('motion_notify_event', self.onPosition)
        
        # create the major objects:
        self.Contour, self.ContourF, self.ContourLabel, self.Vector = None,None,None,None
        self.Grid, self.Image, self.Map = None, None, None
        self.Particles = {'line':[],'txt':[],'data':[],'color':(255,0,0)}
        self.curOri,self.curLayer = 'Z',0
        self.mesh,self.triangles = None,None
        
    #####################################################################
    #                     Divers accesseur/mutateurs
    #####################################################################
    
    def GetToolBar(self):return self.toolbar
    def getcurVisu(self): return [self.curGroupe,self.curNom,self.curObj]
    def onPosition(self,evt):
        sx,sy = str(evt.xdata)[:6],str(evt.ydata)[:6]
        sz = self.gui.guiShow.getPointValue(evt.xdata,evt.ydata)
        self.gui.onPosition(' x: '+sx +' y: '+sy+'  z: '+sz)
        
    def delAllObjects(self):
        for v in self.listeZone:
            self.listeZone[v]=[]
            self.listeZoneText[v]=[]                
        self.cnv.lines=[]
        self.cnv.collections = [None,None]
        self.cnv.artists = []
        self.cnv.images = []
        self.cnv.cla()
        self.draw()
        
    def setVisu(self, core):
        """creer les objets graphiques a partir des caracteristiques d'un modele
        importe.
        creer les zones avec setAllzones, puis le contour et les vectors ecoult,
        les lignes etle contour pour tracer, contour pour reaction
        depend de l'etat du systeme de la liste graphique
        comme ca tout ca pourra etre visualise sans faire de nouveau calcul
        """
        self.delAllObjects()
        for mod in self.core.modelList:
            self.setAllZones(core.diczone[mod].dic)
        self.initDomain()
        self.draw()
        
    def setDataOn(self,bool):
        """definit l'affichage ou non des donnees qaund contour"""
        self.dataon=bool
        
    def redraw(self):
        #self.cnv.set_xlim(self.xlim)
        #self.cnv.set_ylim(self.ylim)
        self.draw()

#    def changeTitre(self,titre):
#        s='';ori=self.curOri
#        if ori in ['X','Y','Z']:
#            plan=self.curPlan;
#            x1,y1 = self.model.Aquifere.getXYticks()
#            zl = self.model.Aquifere.getParm('zList')
#        if ori=='Z': s=' Z = '+ str(zl[plan])
#        if ori=='X': s=' X = '+ str(x1[plan])
#        if ori=='Y': s=' Y = '+ str(y1[plan])
#        pl.title(self.traduit(str(titre))+s[:9],fontsize=20)

    def createAndShowObject(self,dataM,dataV,opt,value=None,color=None):
        """create the Contour, Vector, dataM are data for matrix, V for vector,
        opt is contour or vector
        """
        #print('qtVis',self.gui.guiShow.swiImg,dataM)
        if dataM == None : self.drawContour(False)
        else : # modifs below 22/8/19
            if self.gui.guiShow.swiImg == 'Contour':
                self.createContour(dataM,value,color)
            elif self.gui.guiShow.swiImg == 'Image':
                self.createImage(dataM)
                self.drawContour(False) # OA 22/8/19
            elif self.gui.guiShow.swiImg == 'ImageContour':
                self.createImage(dataM)
                self.drawContour(True) # OA 22/8/19
        if dataV == None : self.drawVector(False)
        else : self.createVector(dataV)
        
    def drawObject(self,typObj,bool):
        if typObj == 'Map' and self.Map == None and bool == False : return
        exec('self.draw'+typObj+'('+str(bool)+')')
        
    def changeObject(self,groupe,name,value,color):
        if name=='Grid':self.changeGrid(color)
        elif name=='Particles': self.changeParticles(value=value,color=color)
        elif name=='Veloc-vect': self.changeVector(value,color)
       #elif name=='Visible': self.changeData(value,color)
        else :self.changeContour(value,color)
        
    #####################################################################
    #             Gestion de l'affichage de la grid/map
    #####################################################################
    # methode qui change la taille du domaine d'etude (les values de l'axe 
    # de la figure matplotlib en fait) et la taille des cellules d'etude
    def initDomain(self):
        # change value of the axes
        grd = self.core.addin.getFullGrid()
        self.xlim = (grd['x0'],grd['x1']);
        self.ylim = (grd['y0'],grd['y1']);
        p,= pl.plot([0,1],'b')
        p.set_visible(False)
        self.transform=p.get_transform()
        self.cnv.set_xlim(self.xlim)
        self.cnv.set_ylim(self.ylim)
        if type(self.core.Zblock) != type(array(5)): 
            self.core.Zblock = makeZblock(self.core)
        grp = self.core.addin.getModelGroup()
        self.mesh=None;#print('in qtvisu unstr:',self.core.mfUnstruct)
        if grp[:2]=='Op' and self.core.dicval[grp+'Flow']['domn.1'][0]>0:
            self.mesh = self.core.addin.opgeo
        if grp[:2]=='Mi' and self.core.dicval[grp+'Flow']['spat.4'][0]>0:
            self.mesh = self.core.addin.min3p
        if self.core.mfUnstruct: 
            self.mesh = self.core.addin.mfU
        self.createGrid()
        # add basic vector as a linecollection
        dep=rand(2,2)*0.;arr=dep*1.
        self.Vector = LineCollection(list(zip(dep,arr)))
        self.Vector.set_transform(self.transform)
        self.Vector.set_visible(False)
        #pl.setp(lc,linewidth=.5);
        self.cnv.collections.append(self.Vector)
        self.Vector.data = [0,0,None,None]

#    def changeDomain(self):
#        self.changeAxesOri('Z',0)
#        
    def changeAxesOri(self,ori):
        # change orientation de la visu
        zb = self.core.Zblock
        self.curOri = ori
        zlim = (amin(zb),amax(zb))
        if ori=='Z': 
            self.cnv.set_xlim(self.xlim);self.cnv.set_ylim(self.ylim)
        elif ori=='X':
            self.cnv.set_xlim(self.ylim);self.cnv.set_ylim(zlim)
        elif ori=='Y':
            self.cnv.set_xlim(self.xlim);self.cnv.set_ylim(zlim)
        self.draw()
        
    def createGrid(self,col=None,ori='Z',layer=0): 
        if self.Grid == None:
            col=(.6,.6,.6);self.Grid=[0,0,col];
        else :
            for i in range(2): self.Grid[i].set_visible(False)
        if col == None: col=self.Grid[2]
        else : self.Grid[2]=col
        if len(self.cnv.collections)<2: self.cnv.collections=[0,0]

        if self.mesh != None and self.core.getValueFromName('Modflow','MshType')>0:# case irregular mesh (from matplotlib2dviewer)
            xcoo = self.mesh.elx
            ycoo = self.mesh.ely
            lines = [list(zip(x,y)) for x,y in list(zip(xcoo,ycoo))]
            self.Grid[0] = LineCollection(lines)
            self.Grid[1] = LineCollection(lines[:2])
            self.Triangles = self.mesh.trg

        else: #rectangular cases
            xl,yl = getXYmeshSides(self.core,self.curOri,self.curLayer)
            dep=list(zip(ravel(xl[:-1,:]),ravel(yl[:-1,:])))
            arr=list(zip(ravel(xl[1:,:]),ravel(yl[1:,:])))
            self.Grid[0] = LineCollection(list(zip(dep,arr)))
            dep=list(zip(ravel(xl[:,:-1]),ravel(yl[:,:-1])))
            arr=list(zip(ravel(xl[:,1:]),ravel(yl[:,1:])))
            self.Grid[1] = LineCollection(list(zip(dep,arr)))
          
        self.cnv.collections[0]= self.Grid[0]
        self.cnv.collections[1]= self.Grid[1]
        for i in [0,1]: 
            self.Grid[i].set_transform(self.transform)
            self.Grid[i].set_color(col);
        self.redraw()
        
    def drawGrid(self, bool):
        col = self.Grid[2]
        if self.curOri != 'Z':
            self.createGrid(col=col,ori=self.curOri,layer=self.curLayer)
        for i in [0,1]: 
            self.Grid[i].set_visible(bool)
            self.Grid[i].set_color(col)
        self.redraw()
        
    def changeGrid(self,color):
        a=color.Get();col=(a[0]/255,a[1]/255,a[2]/255)
        for i in [0,1]: self.Grid[i].set_color(col)
        self.Grid[2]=col
        self.redraw()
    #####################################################################
    #             Draw a map or a variable as an image
    #####################################################################
    # l'image se met en position 1 dans la liste des images
    def createMap(self):
        file = self.gui.map
        mat=Im.imread(file)
        org='upper';ext=(self.xlim[0],self.xlim[1],self.ylim[0],self.ylim[1])
        self.Map=pl.imshow(mat,origin=org,extent=ext,aspect='auto',interpolation='nearest');
        self.cnv.images=[self.Map] #
        self.cnv.images[0].set_visible(True)
        self.redraw()
        
    def drawMap(self, bool):
        if self.Map == None: self.createMap()
#        self.Map.set_visible(bool)
        self.cnv.images=[self.Map] #
        self.cnv.images[0].set_visible(bool)
        self.redraw()

    def createImage(self,data):
        #print 'vis img',len(xt),len(yt),shape(mat)
        X,Y,Z = data;Z1 = array(Z);#print shape(Z1)
        modgroup = self.core.addin.getModelGroup();#print 'visu l 301',modgroup
        if self.mesh==None: # classical square grid
            image=pl.pcolormesh(X,Y,Z,cmap='jet') #,norm='Normalize') #EV 27/08/19
            self.cnv.images=[image]
        elif modgroup[:5] == 'Modfl' :
            data = data[-1][0]# first values are coordinates, and its 3D
            data_sc = (data-amin(data))/(amax(data)-amin(data));#print data
            vts=self.mesh.myv.cverts # for mfU
            xyv=self.mesh.myv.xyverts
            nnod,a = shape(self.mesh.myv.nodes)
            pc=[list(zip(xyv[vts[i],0],xyv[vts[i],1])) for i in range(nnod)]
            self.polyC = PolyCollection(pc,facecolors=cm.jet(data_sc))
            if len(self.cnv.collections)==3:
                self.cnv.collections.append(self.polyC) #self.polyC
            else :
                self.cnv.collections[3] = self.polyC
            self.polyC.set_visible(True)
            self.polyC.set_transform(self.transform)
        elif modgroup[:5] == 'Opgeo' :
            opgeo = self.core.addin.opgeo
            image=plt.tripcolor(self.Triangles,Z1[0],shading='gouraud')#PolyCollection(polys) #plt.tripcolor(triang, z, shading='gouraud', cmap=plt.cm.rainbow)
            #image.set_transform(self.transform)
            #image.set_edgecolor('none')
        #print 'qtVis 336',self.cnv.collections,self.fig.dpi             
        self.redraw()
        
    def drawImage(self,bool):
        if len(self.cnv.images)>0:
            self.cnv.images[0].set_visible(bool)
            self.redraw()

    #####################################################################
    #             Gestion de l'affichage des contours
    #####################################################################
        
    def createContour(self,data,value=None,col = None):
        """ calculation of a contour deom value : value[0] : min
        [1] : max, [2] nb contours, [3] decimals, [4] : 'lin', log' or 'fix',
        if [4]:fix, then [5] is the series of contours"""
        X,Y,Z = data; #print 'visu controu',value,col
        self.cnv.collections=self.cnv.collections[:3]
        self.cnv.artists = []
        V = 11;Zmin=amin(amin(Z));Zmax=amax(amax(Z*(Z<1e5)));
        if Zmax==Zmin : # test min=max -> pas de contour
            onMessage(self.gui,' values all equal to '+str(Zmin))
            return
        if value == None or len(value)<4: value = [Zmin,Zmax,(Zmax-Zmin)/10.,2,'auto',[]]
        # adapt the number and values of the contours
        val2=[float(a) for a in value[:3]]
        if value[4]=='log':  # cas echelle log
            n = int((log10(val2[1])-log10(max(val2[0],1e-4)))/val2[2])+1
            V = logspace(log10(max(val2[0],1e-4)),log10(val2[1]),n)
        elif (value[4]=='fix') and (value[5]!=None) : # fixes par l'utilisateur
            V = value[5]*1;V.append(V[-1]*2.);n=len(V)#;print 'visu 366 ',V
        elif value[4]=='lin' :  # cas echelle lineaire
            n = int((val2[1]-val2[0])/val2[2])+1
            V = linspace(val2[0],val2[1],n)
        else : # cas automatique
            n=11;V = linspace(Zmin,Zmax,n)
        # ONE DIMENSIONAL
        if self.mesh==None:
            r,c=shape(X);
            if r==1: 
                X=concatenate([X,X]);Y=concatenate([Y-Y*.45,Y+Y*.45]);Z=concatenate([Z,Z])
        Z2=ma.masked_where(Z.copy()>1e5,Z.copy());#print value,n,V
        # definir les couleurs des contours
        if col==None or type(col)<type(5): # or (col==[(0,0,0),(0,0,0),(0,0,0),10]):
            cmap = mpl.cm.jet
            col=[(0,0,255),(0,255,0),(255,0,0),10]
        else :
            r,g,b=[],[],[];
            lim=((0.,1.,0.,0.),(.1,1.,0.,0.),(.25,.8,0.,0.),(.35,0.,.8,0.),(.45,0.,1.,0.),\
                 (.55,0.,1.,0.),(.65,0.,.8,0.),(.75,0.,0.,.8),(.9,0.,0.,1.),(1.,0.,0.,1.))
            for i in range(len(lim)):
                c1=lim[i][1]*col[0][0]/255.+lim[i][2]*col[1][0]/255.+lim[i][3]*col[2][0]/255.
                r.append((lim[i][0],c1,c1))
                c2=lim[i][1]*col[0][1]/255.+lim[i][2]*col[1][1]/255.+lim[i][3]*col[2][1]/255.
                g.append((lim[i][0],c2,c2))
                c3=lim[i][1]*col[0][2]/255.+lim[i][2]*col[1][2]/255.+lim[i][3]*col[2][2]/255.
                b.append((lim[i][0],c3,c3))
            cdict={'red':r,'green':g,'blue':b}
            cmap=mpl.colors.LinearSegmentedColormap('my_colormap', cdict, 256)
        if self.mesh==None or self.core.getValueFromName('Modflow','MshType')<1: # OA 29/2/20
            cf = pl.contourf(pl.array(X),pl.array(Y),Z2,V, cmap=cmap)
            c = pl.contour(pl.array(X),pl.array(Y),Z2,V, cmap=cmap)
        else :
            #print 'qtVi 385',shape(Z2),shape(V),shape(self.Triangles.x)
            cf = pl.tricontourf(self.Triangles,Z2,V, cmap=cmap)
            c = pl.tricontour(self.Triangles,Z2,V, cmap=cmap)            
        #print col[3]
        for c0 in cf.collections:
            c0.set_alpha(int(col[3])/100.);#print cl
        if value==None: fmt = '%1.3f' 
        else : fmt = '%1.'+ str(int(value[3]))+'f'
        #print fmt
        cl = pl.clabel(c,colors='black',fontsize=9,fmt=fmt)
        self.Contour = c
        self.ContourF = cf
        self.ContourLabel = cl
        self.Contour.data = data
        self.draw()
            
    def changeContour(self,value,col):
        """ modifies the values of an existing contour"""
        self.drawContour(False)
        self.createContour(self.Contour.data,value,col)

    def drawContour(self,bool):
        # OA 22/8/19 modified to have true and false
        if self.Contour == None : return # added 2/10/19
        for c in self.Contour.collections :c.set_visible(bool)
        for c in self.ContourF.collections :c.set_visible(bool)
        for c in self.ContourLabel: c.set_visible(bool)
        self.draw()

    #####################################################################
    #             Gestion de l'affichage de vectors
    #####################################################################
    """vector has been created as the first item of lincollection list
    during domain intialization"""
    def createVector(self,data):
        X,Y,U,V = data
        if self.mesh !=None :
            X,Y = self.mesh.nodes[:,1],self.mesh.nodes[:,2]; #print X,Y,U,V
        """ modifie les values de vectors existants"""
        if self.Vector.data[3] == None: #first vector no color
            a=ravel(X);ech =( max(a)-min(a))/40;col=(0,0,1)
        else : 
            a,b,ech,col = self.Vector.data
            self.drawVector(False)
        l=len(ravel(X))
        dep=concatenate([X.reshape((l,1)),Y.reshape((l,1))],axis=1)
        b=X+U*ech;c=Y+V*ech;
        arr=concatenate([b.reshape((l,1)),c.reshape((l,1))],axis=1)
        self.Vector = LineCollection(list(zip(dep,arr)))
        self.Vector.set_transform(self.transform)
        self.Vector.set_color(col);
        if len(self.cnv.collections)>2:self.cnv.collections[2]=self.Vector
        else : self.cnv.collections.append(self.Vector)
        self.Vector.set_visible(True)
        self.Vector.data = [dep,arr,ech,col];#print self.Vector.data
        self.redraw()

    def drawVector(self, bool):
        """ dessine les vectors vitesse a partir de x,y,u,v et du
        booleen qui dit s'il faut dessiner ou non """
        self.Vector.set_visible(bool);#print self.Vector
        self.redraw()
        
    def changeVector(self,ech,col=(0,0,255)):
        """ modifie les values de vectors existants"""
        #self.drawVector(False)
        ech = float(ech)
        #change coordinates
        dep,arr_old,ech_old,col_old = self.Vector.data;#print shape(dep),shape(arr_old),ech,ech_old
        arr=dep+(arr_old-dep)*ech/ech_old
        # new object
        #self.Vector = LineCollection(zip(dep,arr))
        self.Vector.set_segments(list(zip(dep,arr)))
        #self.Vector.set_transform(self.transform)
        #a=col.Get();col=(a[0]/255,a[1]/255,a[2]/255)
        self.Vector.set_color(col);
        self.Vector.set_visible(True)
        #self.cnv.collections[2]=self.Vector
        self.Vector.data = [dep,arr,ech,col]
        self.redraw()
    
    #####################################################################
    #             Gestion de l'affichage de particules
    #####################################################################
    def startParticles(self):
        if self.Particles != None:
            self.partVisible(False)
        self.Particles = {'line':[],'txt':[],'data':[],'color':(255,0,0)}
        self.mpl_disconnect(self.toolbar._idPress)
        self.mpl_disconnect(self.toolbar._idRelease)
        self.mpl_disconnect(self.toolbar._idDrag)
        # on capte le clic gauche de la souris
        self.m3 = self.mpl_connect('button_press_event', self.mouseParticles)
        self.stop = False

    def mouseParticles(self, evt):
        #test pour savoir si le curseur est bien dans les axes de la figure
        if self.stop : return
        if evt.inaxes is None: return
        if evt.button==1:
            self.core.addin.calcParticle(evt.xdata,evt.ydata);#print xp,yp,tp
            [xp,yp,tp] = self.core.addin.particle['data']
            self.updateParticles(xp,yp,tp)
        elif evt.button==3:
            self.mpl_disconnect(self.m3);self.stop =True
            self.gui.actions('zoneEnd')

    def updateParticles(self,X,Y,T,freq=10):
        """ adds a line in the gropu of particles"""
        ligne, = pl.plot(pl.array(X),pl.array(Y),'r');
        if freq>0 and type(T) == type(ones(5)): # OA 18/3/19
            tx,ty,tt = X[0::freq],Y[0::freq],T[0::freq]
            txt = []
            for i in range(len(tx)):
                a=str(tt[i]);b=a.split('.');ln=max(4,len(b[0]))
                txt.append(pl.text(tx[i],ty[i],a[:ln],fontsize='8'))
            self.Particles['txt'].append(txt)
        self.Particles['line'].append(ligne)
        self.Particles['data'].append((X,Y,T))
        #self.gui_repaint() # bug matplotlib v2.6 for direct draw!!!
        self.draw()
        
    def drawParticles(self,bool,value=None):
        if self.Particles==None: return
        self.partVisible(bool)
        #self.gui_repaint()
        self.draw()
                
    def changeParticles(self,value=None,color=(255,0,0)):
        self.partVisible(False)
        self.Particles['color'],self.Particles['txt'] = color,[]
        for i,data in enumerate(self.Particles['data']):
            X,Y,T = data
            self.Particles['line'][i].set_data(X,Y)
            try:
                tx,ty,tt = self.ptsPartic(X,Y,T,float(value))
                txt = []
                for i in range(len(tx)):
                    a=str(tt[i]);b=a.split('.');ln=max(4,len(b[0]))
                    txt.append(pl.text(tx[i],ty[i],a[:ln],fontsize='8'))
                self.Particles['txt'].append(txt)
            except: 
                a=1
        self.partVisible(True)
        #self.gui_repaint()  #OA removed 20/11/19
        self.draw()
        
    def partVisible(self,bool):
        a=(255,0,0) #self.Particles['color'].Get()
        color=(a[0]/255,a[1]/255,a[2]/255)
        for line in self.Particles['line']:
            line.set_visible(bool);line.set_color(color)
        for points in self.Particles['txt']: 
            for txt in points : txt.set_visible(bool)
        
    def ptsPartic(self,X,Y,T,dt):
        #tx,ty,tt,i1=iphtC1.ptsLigne(X,Y,T,dt);
        tmin=amin(T);tmax=amax(T);
        t1=linspace(tmin,tmax,int((tmax-tmin)/dt))
        f=interp1d(T,X);xn=f(t1)
        f=interp1d(T,Y);yn=f(t1)
        return xn,yn,t1
    #####################################################################
    #                   Gestion des zones de la visu
    #####################################################################
    # affichage de toutes les zones d'une variable
    def showVar(self, var, media):
        self.setUnvisibleZones()
        self.curVar, self.curMedia =var, media
        for i in range(len(self.listeZone[var])):
            #print 'vis showvar',len(self.listeZone[var]),len(self.listeZmedia[var][i])
            if (media in self.listeZmedia[var][i]) or (media==-1):
                self.listeZone[var][i].set_visible(True)
                self.visibleText(self.listeZoneText[var][i],True)
        #self.changeTitre(var)
        self.redraw()
        
    def showData(self,liForage,liData):
        self.setUnvisibleZones();self.curVar='data'
        self.listeZoneText['data']=[]
        for zone in self.listeZone['Forages']: zone.set_visible(True)
        lZone=self.model.Aquifere.getZoneList('Forages');txt=[]
        for z in lZone:
            x,y=list(zip(*z.getXy()));name=z.getNom()
            if name in liForage:
                ind=liForage.index(name)
                txt.append(pl.text(mean(x),mean(y),name+'\n'+str(liData[ind])))
        obj = GraphicObject('zoneText', txt, True, None)        
        self.addGraphicObject(obj)
        self.redraw()
        
    def changeData(self,taille,col):
        obj=self.listeZoneText['data'][0].getObject()
        for txt in obj:
            txt.set_size(taille);txt.set_color(col)
                
    def getcurZone(self) : return self.curZone
    def setcurZone(self,zone) : self.curZone = zone
    # methode qui efface toutes les zones de toutes les variables   
    def setUnvisibleZones(self):
        for v in self.listeZone:
            for zone in self.listeZone[v]: zone.set_visible(False)
            for txt in self.listeZoneText[v]: 
                if type(txt)==type([5,6]):
                    for t in txt : t.set_visible(False)
                else : txt.set_visible(False)

    # method called by the GUI to create a new zone
    def setZoneReady(self,typeZone, curVar):
        self.typeZone = typeZone
        self.curVar = curVar
        self.tempZoneVal = []
        # on deconnecte la toolbar pour activer la formaiton de zones
        self.mpl_disconnect(self.toolbar._idPress)
        self.mpl_disconnect(self.toolbar._idRelease)
        self.mpl_disconnect(self.toolbar._idDrag)
        # on capte le clic gauche de la souris
        self.m1 = self.mpl_connect('button_press_event', self.mouse_clic)
        
    def setZoneEnd(self,evt):
        # the GUI is informed about end of drawing
        xv, yv = self.getcurZone().get_xdata(),self.getcurZone().get_ydata()
        #dx,dy = (amax(xv)-amin(xv))/25, (amax(yv)-amin(yv))/25 # closing the poly if needed
        #if (abs(xv[0]-xv[1])<dx) & (abs(yv[0]-yv[1])<dy): 
        #    xv[-1], yv[-1] = xv[0], yv[0]
        if len(self.tempZoneVal)>1: xy = list(zip(xv,yv,self.tempZoneVal))
        else : xy = list(zip(xv,yv))
        # remove zone if cancel, reorder
        self.curZone.set_visible(False)
        self.curZone = None
        self.x1,self.y1 = [],[]
        self.gui.addBox.onZoneCreate(self.typeZone, xy)
        
    def addZone(self,media,name,val,coords,visible=True):
        """ add a zone in visu 
        """
        #print 'visu',coords
        a = list(zip(*coords)); txt = [];#print name,a
        if len(a)==0:return
        if len(a)==2: x,y = a
        elif len(a)==3: x,y,z = a
        if len(x)==1: 
            zone = Line2D(x, y,marker='+',markersize=10,markeredgecolor='r')
        else : 
            zone = Line2D(x, y)
        zone.verts = coords
        zone.set_visible(visible)
        if type(media)!=type([2]): media = [int(media)]
        #self.curMedia = media
        self.cnv.add_line(zone)
        if self.typeZone == "POLYV" or len(coords[0])==3:
            txt = [pl.text(mean(x)*.1+x[0]*.9,mean(y)*.1+y[0]*.9,name+'\n'+str(val)[:16])]
            for i in range(len(x)):
                t=pl.text(x[i],y[i],str(z[i]))
                t.set_visible(visible)
                txt.append(t)
        else :
            txt = pl.text(mean(x)*.1+x[0]*.9,mean(y)*.1+y[0]*.9,name+'\n'+str(val)[:16])
            txt.set_visible(visible)
        self.listeZone[self.curVar].append(zone)
        self.listeZmedia[self.curVar].append(media)
        self.listeZoneText[self.curVar].append(txt)
        if visible: self.redraw()
        
    def delZone(self, Variable, ind):
        """methode de suppression de la zone d'indice ind de Variable
        """
        if (Variable in self.listeZone)==False: return
        if len(self.listeZone[Variable])>ind:
            self.listeZone[Variable][ind].set_visible(False)
            self.visibleText(self.listeZoneText[Variable][ind],False)
            self.listeZone[Variable][ind:ind+1] = []
            self.listeZoneText[Variable][ind:ind+1] = []
            self.listeZmedia[Variable].pop(ind)
            self.redraw()
            
    def visibleText(self,text,bool):
        if type(text)==type([5,6]): 
            for t in text : t.set_visible(bool)
        else :
            text.set_visible(bool)
            
    def delAllZones(self,Variable):
        lz=self.listeZone[Variable]
        for i in range(len(lz)):
            lz[i].setVisible(False);
            self.listeZoneText[Variable][i].set_visible(False)
        self.listeZone[Variable] = []
        self.listeZmedia[Variable] = []
        self.listeZoneText[Variable] = []
        self.redraw()
        
    def modifValZone(self, nameVar, ind, val,xy):
        """modify the value (or list of value) for the zone nameVar 
        the text contains name et value"""

        
    def modifZoneAttr(self,nameVar,ind,val,media,xy):
        # modify xy
        zone = self.listeZone[nameVar][ind]
        if len(xy[0])==3: x,y,z =list(zip(*xy))
        else : x,y =list(zip(*xy))
        zone.set_data(x,y)
        # modify media
        if type(media)!=type([2]): media = [int(media)]
        self.listeZmedia[nameVar][ind] = media
        # modify text
        textObj = self.listeZoneText[nameVar][ind];
        if type(textObj)==type([2,3]):
            name = pl.getp(textObj[0],'text').split('\n')[0]
            pl.setp(textObj[0],text=name+'\n'+str(val)[:16])
            for i in range(len(z)):
                pl.setp(textObj[i+1],text=str(z[i]))
        else:
            name = pl.getp(textObj,'text').split('\n')[0]
            pl.setp(textObj,text=name+'\n'+str(val)[:16])
        self.redraw()

    def modifZone(self, nameVar, ind):
        """ modification interactive des points de la zone d'indice ind de name nameVar
        """
        zone = self.listeZone[nameVar][ind]
        self.polyInteract = PolygonInteractor(self,zone, nameVar, ind)
        zone.set_visible(False)
        self.cnv.add_line(self.polyInteract.line)
        self.draw()

    def finModifZone(self):
        """fonction qui met fin a la modification de la zone courante"""
        if self.polyInteract != None:
            self.polyInteract.set_visible(False)
            self.polyInteract.disable()
            # on informe la GUI des nouvelles coordonnees
            var, ind = self.polyInteract.typeVariable,self.polyInteract.zind #OA 6/2/20 ind-zind
            x,y=self.polyInteract.lx,self.polyInteract.ly;#print x,y
            xy=list(zip(x,y));self.gui.modifBox.onModifZoneCoord(var, ind, xy)
            zone =self.listeZone[var][ind]
            zone.set_data(x,y);zone.set_visible(True)
            # on modifie la position du texte
            txt = self.listeZoneText[var][ind]
            if type(txt)==type([5,6]):
                txt[0].set_position((x[0],y[0]))
                for i in range(1,len(txt)):
                    txt[i].set_position((x[i-1],y[i-1]))
            else:
                txt.set_position((mean(x)*.1+x[0]*.9,mean(y)*.1+y[0]*.9))                
            self.draw()            
    
    def setAllZones(self,dicZone):
        """updates all zones when a file is imported
        """
        for var in list(dicZone.keys()):
            self.listeZone[var]=[]
            self.curVar = var
            lz = dicZone[var];
            nbz = len(lz['name'])
            for i in range(nbz):
                if lz['name'][i]=='': continue
                coords = lz['coords'][i]
                self.addZone(lz['media'][i],lz['name'][i],lz['value'][i],coords,False)
        self.setUnvisibleZones()
        self.redraw()

    #####################################################################
    #             Gestion de l'interaction de la souris
    #             pour la creation des zones
    #####################################################################    

    #methode executee lors d'un clic de souris dans le canvas
    def mouse_clic(self, evt):
        if evt.inaxes is None:
            return      
        if self.curZone == None:  # au depart
            self.x1 = [float(evt.xdata)] 
            self.y1 = [float(evt.ydata)]
            self.setcurZone(Line2D(self.x1,self.y1))
            self.cnv.add_line(self.curZone)
            self.m2 = self.mpl_connect('motion_notify_event', self.mouse_motion)
            if self.typeZone=="POLYV":
                self.polyVdialog()
            if self.typeZone=="POINT":
                self.deconnecte()
                self.setZoneEnd(evt)

        else :  # points suivants
            if self.typeZone=="POLYV": # and evt.button ==1:
                rep = self.polyVdialog()  # dialog for the current value of z
                if rep==False : return
            self.x1.append(float(evt.xdata))# OA 14/2 removed [:6] for long coordinates
            self.y1.append(float(evt.ydata))
            if self.typeZone=="LINE" or self.typeZone=="RECT" :
                self.deconnecte()  #fin des le 2eme point
                self.setZoneEnd(evt)
            if self.typeZone in ["POLY","POLYV"] and evt.button==3: # fin du polygone
                self.deconnecte()
                self.setZoneEnd(evt)

    #done after a mouse clic durng mouse movt in the canvas
    def mouse_motion(self, evt):
        time.sleep(0.1)
        if evt.inaxes is None: return
        lx, ly = self.x1*1, self.y1*1    
        if self.typeZone == "RECT":
            xr,yr = self.creeRectangle(self.x1[0],self.y1[0],evt.xdata,evt.ydata)
            self.curZone.set_data(xr,yr)      
        else : # autres cas
            lx.append(evt.xdata);ly.append(evt.ydata)
            self.curZone.set_data(lx,ly)
        self.draw()

    def polyVdialog(self):
        lst0=[('Elevation','Text',0)] # OA 25/4/19 changed value to elevation
        dialg = genericDialog(self.gui,'value',lst0)
        values = dialg.getValues()
        if values != None:
            val = float(values[0]);#print val*2
            self.tempZoneVal.append(val)
            return True
        else: 
            return False

    def creeRectangle(self, x1, y1, x2, y2):
            xr=[x1,x2,x2,x1,x1]
            yr=[y1,y1,y2,y2,y1]
            return [xr,yr]
        
    def deconnecte(self):
        # deconnecter la souris
        self.mpl_disconnect(self.m1)
        self.mpl_disconnect(self.m2)
    ###################################################################
    #   deplacer une zone ##############################

    def startMoveZone(self, nameVar, ind):
        """ methode qui demarre les interactions avec la souris"""
        # reperer la zone et rajouter un point de couleur
        self.nameVar, self.ind = nameVar, ind;
        zone = self.listeZone[nameVar][ind]
        self.curZone = zone
        self.lx, self.ly = zone.get_xdata(),zone.get_ydata()
        self.xstart=self.lx[0]*1.;self.ystart=self.ly[0]*1.;
        self.ptstart=Line2D([self.xstart],[self.ystart],marker='o',
                markersize=7,markerfacecolor='r')
        self.cnv.add_line(self.ptstart)
        self.m1 = self.mpl_connect('button_press_event', self.zoneM_clic)
        self.draw()
        
    def zoneM_clic(self,evt):
        """ action au premier clic"""
        if evt.inaxes is None: return
        #if evt.button==3: self.finMoveZone(evt) # removed OA 6/2/13
        d=sqrt((evt.xdata-self.xstart)**2+(evt.ydata-self.ystart)**2)
        xmn,xmx=self.xlim;ymn,ymx=self.ylim
        dmax=sqrt((xmx-xmn)**2+(ymx-ymn)**2)/100;
        if d>dmax: return
        self.m2 = self.mpl_connect('motion_notify_event', self.zone_motion)
        self.m3 = self.mpl_connect('button_release_event', self.finMoveZone)
        self.mpl_disconnect(self.m1)

    def zone_motion(self,evt):
        """ methode pour deplacer la zone quand on deplace la souris"""
        # reperer le curseur proche du point de couleur
        time.sleep(0.1);
        if evt.inaxes is None: return
        # changer els coord du polygone lorsque l'on deplace la souris
        lx=[a+evt.xdata-self.xstart for a in self.lx]
        ly=[a+evt.ydata-self.ystart for a in self.ly]
        self.ptstart.set_data(lx[0],ly[0])
        self.curZone.set_data(lx,ly);
        self.draw()
        
    def finMoveZone(self,evt):
        """ methode pour arret de deplacement de la zone"""
        # lorsque l'on relache la souris arreter les mpl connect
        self.mpl_disconnect(self.m2)
        self.mpl_disconnect(self.m3)
        # renvoyer les nouvelels coordonnes au modele
        lx, ly = self.curZone.get_xdata(),self.curZone.get_ydata()
        self.listeZone[self.nameVar][self.ind].set_data(lx,ly)
        xy = list(zip(lx,ly))
        self.gui.modifBox.onModifZoneCoord(self.nameVar, self.ind, xy)
        # on modifie la position du texte
        txt = self.listeZoneText[self.nameVar][self.ind]
        if type(txt)==type([5,6]):
            txt[0].set_position((lx[0],ly[0]))
            for i in range(1,len(txt)): 
                txt[i].set_position((lx[i-1],ly[i-1])) #-1 because 1st position zone names
        else:
            txt.set_position((mean(lx)*.1+lx[0]*.9,mean(ly)*.1+ly[0]*.9))                
        self.ptstart.set_visible(False);self.ptstart = None
        self.curZone=None
        self.draw()

class PolygonInteractor:
    """
    A polygon editor.

    Key-bindings
      't' toggle vertex markers on and off.  When vertex markers are on,
          you can move them, delete them
      'd' delete the vertex under point      
      'i' insert a vertex at point.  You must be within epsilon of the
          line connecting two existing vertices          
    """

    showverts = True

    def __init__(self, gui,poly, typeVariable, zind):
        if poly.figure is None:
            raise RuntimeError('You must first add the polygon to a figure or canvas before defining the interactor')
        self.canvas,self.gui,self.zind = poly.figure.canvas,gui,zind #OA 6/2/20
        self.poly = poly
        self.epsilon=(gui.xlim[1]-gui.xlim[0])/100
        x, y = self.poly.get_xdata(),self.poly.get_ydata()
        self.lx=list(x);self.ly=list(y)
        self.line = Line2D(x,y,marker='o', markerfacecolor='r')
        self.typeVariable = typeVariable
        
        cid = self.poly.add_callback(self.poly_changed)
        self._ind = None # the active vert
        self.canvas.setFocusPolicy( Qt.ClickFocus) # OA added for Qt to focus on keyboard 19/11/19
        self.canvas.setFocus()

        self.c1 = self.canvas.mpl_connect('button_press_event', self.button_press_callback)
        self.c2 = self.canvas.mpl_connect('key_press_event', self.key_press_callback)        
        self.c3 = self.canvas.mpl_connect('button_release_event', self.button_release_callback)
        self.c4 = self.canvas.mpl_connect('motion_notify_event', self.motion_notify_callback)                

    def poly_changed(self, poly):
        'this method is called whenever the polygon object is called'
        # only copy the artist props to the line (except visibility)
        vis = self.line.get_visible()
        Artist.update_from(self.line, poly)
        self.line.set_visible(vis)  # don't use the poly visibility state        

    def get_ind_under_point(self, event):
        'get the index of the vertex under point if within epsilon tolerance'
        x, y = self.lx,self.ly
        d = sqrt((pl.array(x)-event.xdata)**2 + (pl.array(y)-event.ydata)**2)
        indseq = nonzero(equal(d, amin(d)));ind = indseq[0];
        if len(ind)>1: ind=ind[0]
        if d[ind]>=self.epsilon: ind = None
        return ind
        
    def button_press_callback(self, event):
        'whenever a mouse button is pressed'
        if not self.showverts: return 
        if event.inaxes==None: return
        if event.button == 1: 
            self._ind = int(self.get_ind_under_point(event))
        if event.button==3:
            self.disable();self.gui.finModifZone()
        #print(self._ind)

    def button_release_callback(self, event):
        'whenever a mouse button is released'
        if not self.showverts: return
        if event.button != 1: return
        self._ind = None      
        # modification du poly passe en parametre en fonction 
        # des nouveaux sommets crees par le PolygonInteractor
        x = self.line.get_xdata()
        y = self.line.get_ydata()  
        v = self.poly.verts
        for i in range(len(x)):
            v[i] = (x[i],y[i])  
        
    def key_press_callback(self, event):
        'whenever a key is pressed'
        if not event.inaxes: return
        # if event.key=='t':
        #     self.showverts = not self.showverts
        #     self.line.set_visible(self.showverts)
        #     if not self.showverts: self._ind = None
        if event.key=='d':
            ind = self.get_ind_under_point(event);print(ind)
            if ind is not None:
                self.lx.pop(ind[0]);self.ly.pop(ind[0]) # added 0 OA 19/11/19
                self.line.set_data(self.lx,self.ly)
        elif event.key=='i':            
            #xys = self.poly.get_transform().seq_xy_tups(self.poly.verts)
            p = (event.xdata, event.ydata) # display coords
            for i in range(len(self.lx)-1):                
                s0 = (self.lx[i],self.ly[i]); #xys[i]
                s1 = (self.lx[i+1],self.ly[i+1]); #xys[i+1]
                d = dist_point_to_segment(p, s0, s1)
                if d<=self.epsilon:
                    self.lx.insert(i+1,event.xdata)
                    self.ly.insert(i+1,event.ydata)
                    self.line.set_data(self.lx,self.ly)
                    break                          
        self.canvas.draw_idle()

    def motion_notify_callback(self, event):
        'on mouse movement'
        time.sleep(0.1);
        if not self.showverts: return 
        if self._ind is None: return
        if event.inaxes is None: return
        if event.button != 1: return
        x,y = event.xdata, event.ydata;
        self.lx[self._ind],self.ly[self._ind] = x,y
        self.line.set_data(self.lx,self.ly)
        self.canvas.draw_idle()
        
    def disable(self):
        self.canvas.mpl_disconnect(self.c1)
        self.canvas.mpl_disconnect(self.c2)
        self.canvas.mpl_disconnect(self.c3)
        self.canvas.mpl_disconnect(self.c4)

    def set_visible(self,bool):
        self.line.set_visible(bool)
