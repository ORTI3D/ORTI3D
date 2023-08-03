from .config import *
from .geometry import *

class guiShow:
    def __init__(self,gui,core):

        self.core = core
        self.Glist = {}
        #self.SetBackgroundColour('#EDFAFF')
        self.groups={
            'Model':[0,['Plane',['Z','X','Y']],['Layer',['___']], #OA 20/11/19  order
                ['Tstep',['_____']],'Grid','Map'], #,'Variable'],
            'Flow':[1,'Head','Wcontent','Veloc-vect','Veloc-magn','Particles'],
            'Transport':[2,'Tracer','Temperature'],
            'Chemistry':[3,['Species',['_______']],
                     ['Units',['mol/L','mmol/L','umol/L','nmol/L']],
                     ['User',['_______']]],
            #'Observation':[4,['Type',['Profile','Breakthrough','XYplot']],
                        #   ['Zone',['_______']]]}
            'Observation':[4,['Type',
                              ['Time-series Graphs','Horizontal profile Graphs',
                               'Vertical profile Graphs','Calibration Graphs', #EV 20/04/20
                               'Mass balance Graphs','Zone budget Graphs']], # EV 02/03/20 
                           ['Result',['Flow','Transport','Chemistry']]]} # OA 21/2/2019 # EV 02/03/20 'W content',
        self.dicVisu = {'Model':{'Plane':'Z','Layer':0,'Tstep':0,
                                 'Grid':False,'Map':False,'Variable':False},
                'Flow':{'Head':False,'Wcontent':False,'Veloc-vect':False,
                        'Veloc-magn':False,'Particles':False},
                'Transport':{'Tracer':False,'Temperature':False},
                'Chemistry':{'Species':False,'Units':'mmol/L','User':False},
                #'Observation':{'Type':'Profile','Zone':' '}}
                'Observation':{'Type':'Time-series Graphs','Result':'Head'}}
        self.change={'Grid':None,'Veloc-vect':['scale',1.],
            'Particles':['time',10.],
            'Visible':['size',10]} # change that are not contours
        self.Vtypes = {'Modif':['Plane','Layer','Tstep','Units'],
            'Array' : ['Head','Wcontent','Veloc-magn','Tracer','Species','User'],
            'Image': ['Variable','Map'],'Grid':['Grid'],'Particles':['Particles'],
            'Vector' : ['Veloc-vect']}
        self.gui,self.visu,self.Tstep = gui,gui.visu,0  # OA 20/11/19 added tstep
        cfg = Config(self.core)
        self.gtyp = cfg.gtyp
        if self.gtyp=='wx': 
            self.dlgShow = cfg.show.Show(self,self.gui,self.core)
        elif self.gtyp in ['qt','qgis']: 
            self.dlgShow = gui.dlgShow
        self.dialogs = cfg.dialogs
        self.dicplots={'X_head':None,'X_tracer':None}
        self.curName, self.arr3= None, None #EV 11/12/19
        self.init()
        
    def init(self):     
        self.swiImg = 'Contour'
        self.userSpecies,self.data = {},None
        self.visu.Glist = self.Glist
        self.curVar, self.curVarView = {},None
        modgroup = self.core.addin.getModelGroup();self.modgroup = modgroup
        self.MshType = 0
        if 'USG' in modgroup: self.MshType = self.core.getValueFromName('Modflow','MshType')
        if modgroup == 'Openfoam': self.MshType = self.core.getValueFromName('OpenFlow','MshType')
        
    def openModel(self):
        self.init()
        if self.core.addin.getDim() == '3D': self.setNames('Model_Plane_L',['Z','X','Y'])#OA 20/11/19  order
        else : self.setNames('Model_Plane_L',['Z'])
        '''
        if self.core.dicaddin['Model']['group'] == 'Min3p':
            self.dlgShow.getBoxNames('Flow_Wcontent_B',False)
        else : 
            mm,mval = self.core.dicaddin['usedM_Modflow']
            v1, v2 = mval[mm.index('UPW')],mval[mm.index('UZF')]
            mod=self.core.dicaddin['Model']['group']
            if (v1==2 or v2==2) and (mod =='Modflow series'):
                self.gui.guiShow.dlgShow.getBoxNames('Flow_Wcontent_B',False)
            else : self.gui.guiShow.dlgShow.getBoxNames('Flow_Wcontent_B',True)
        '''
        if self.core.dicaddin['Model']['group'] == 'Modflow USG':
            self.gui.guiShow.dlgShow.getBoxNames('Flow_Particles_B',True)
            self.gui.onParticle(False) 
        else : 
            self.gui.guiShow.dlgShow.getBoxNames('Flow_Particles_B',False)
            self.gui.onParticle(True)
            

    def getCurrentTime(self): return self.dlgShow.getCurrentTime()
    def getNames(self,nameBox): return self.dlgShow.getNames(nameBox)
    def setNames(self,nameBox,names,opt='strings'): 
        self.dlgShow.setNames(nameBox,names,opt)
    def setChemSpecies(self,specList):
        self.setNames('Chemistry_Species_L',specList)
    def setUserSpecies(self,dicU): self.userSpecies = dicU.copy();#print dicU

    def getGlist(self,group,name):
        #print('guish getglist',group,name,self.Glist)
        if group in self.Glist:
            if name in self.Glist[group]: return self.Glist[group][name] 
            else :
                self.Glist[group][name] = {'value':None,'color':None}
        else :
            self.Glist[group] = {}
            self.Glist[group][name] = {'value':None,'color':None}
        return self.Glist[group][name] 
    
    def setGlistParm(self,group,name,parm,value):
        #print 'guish 67', group,name,parm, value
        self.Glist[group][name][parm]= value;#print('guishow setglist',group,name,self.Glist[group][name])
        
    # def resetGlist(self):
    #     """reset all Glist tags 'calc' to False"""
    #     for group in list(self.Glist.keys()):
    #         for name in list(self.Glist[group].keys()): 
    #             self.Glist[group][name]['calc'] = False

    def onClick2(self,group,name,retour):
        """after the click the type of object is defined and data are retrieved 
        according to what has been typed
        group is the box, name is the item and retour is wether a boolean
        for a tickbox or the name of the item for a list
        the time and species has been stored previously (in onclick below)
        there are four types of action : 
          -  if a False arrives in retour just make the object unvisible. 
          - for observation go elsewhere to show the observations
          - if plane/time has changed show the same object but for different plane/time
          - for the other case show the object
        objects : grid (true/false), vectors(true/false), map(true/false)
        variable (?), contour(group,name):True/false, types in Vtypes
        """
        #current visu data are stored (by wx dialog) in dicVisu
        self.currentGroup = group
        listSpecies = self.getNames('Chemistry_Species_L')
        opt,bool = 'contour',False
        # set the first steps GRID, MAP, VARIABLE, PARTICLE
        m = self.dicVisu['Model'] 
        plane, layer, tstep = m['Plane'],m['Layer'],m['Tstep'];
        self.Tstep = tstep
        #self.visu.drawObject('Grid',self.dicVisu['Model']['Grid'])
        if m['Variable'] : 
            self.visu.createImage(self.getCurrentVariable(plane,layer))
            m['Map'] = False
        elif m['Map'] : 
            self.visu.drawObject('Map',True)
        elif name == 'Grid' :#EV 15/02/2021
            self.visu.drawObject('Grid',retour)
        else :
            self.visu.drawObject('Image',False)
        self.visu.drawObject('Particles',self.dicVisu['Flow']['Particles'])
        # find the current CONTOUR and if needs to be dranw
        Cgroup,Cname,species = self.getCurrentContour();print('quishow 147',Cname,species)
        self.dlgShow.uncheckContours(Cgroup,Cname,species) # OA 9/6/19
        self.curGroup,self.curName,self.curSpecies = Cgroup,Cname,species;# OA 10/5/17
        if group=='Observation': # observation for the group that is currently drawn
            #if name=='Zone': self.dlgShow.onObservation(Cgroup,tstep)
            if name=='Result': self.dlgShow.onPlot()
            return
        # get the data for contours
        dataM = None
        #print 'guish 116', name,species,self.userSpecies
        if Cgroup != None : 
            self.arr3 = self.getArray3D(Cgroup,Cname,tstep,species)
            if self.arr3 is None : # EV 8/12/21
                mess=onMessage(self.gui,'No result')
                self.dlgShow.onTickBox(group,name,'B',False)
                self.dicVisu[group][name]=False
                return
            if species in list(self.userSpecies.keys()):
                dataM = self.getUserSpecies(species,plane,layer)
            else :
                dataM = self.getArray2D(Cgroup,self.arr3,plane,layer)
        self.data = dataM; #print('guish 133',shape(dataM))
        # get VECTORS
        dataV = None
        if self.dicVisu['Flow']['Veloc-vect']: 
            dataV = self.getVectors(plane,layer,tstep)
            opt = 'vector'
        #print self.dicVisu
        toshow = species
        if type(species)==type(bool): toshow = Cname # 28/3/17 oa to keep contour values for 
        glist = self.getGlist(Cgroup,toshow)
        value,color = glist['value'],glist['color'];#print('guishow 171',Cgroup,Cname,value,color)
        if layer !=0 : self.visu.changeAxesOri(plane)
        self.visu.curLayer = layer
        self.visu.createAndShowObject(dataM,dataV,opt,value,color)

    def redraw(self):
        group,name,species = self.getCurrentContour()
        self.onClick2(group,name,species)
        
    def resetDicContour(self):
        # put all cntour values to none
        for k in list(self.dicVisu.keys()):
            for k1 in list(self.dicVisu[k].keys()):
                if k1 in self.Vtypes['Array']:
                    self.dicVisu[k][k1] = False
 
    def getCurrentContour(self):   
        #returns the contour that is currently visible
        for k in list(self.dicVisu.keys()):
            for k1 in list(self.dicVisu[k].keys()):
                if k1 in self.Vtypes['Array']:
                    if self.dicVisu[k][k1] != False :  
                        return k,k1,self.dicVisu[k][k1]
        return None,None,None
        
    def getArray3D(self,group,name,tstep,spec):
        """get an array of data in the files written by a model, using the reader"""
        #print group,name,tstep
        arr = None
        if group=='Flow':
            if name=='Head':
                arr = self.core.flowReader.readHeadFile(self.core,tstep)
            elif name=='Wcontent':
                arr = self.core.flowReader.readWcontent(self.core,tstep)
            elif name=='Veloc-magn':
                vx,vy,vz = self.core.flowReader.readFloFile(self.core,tstep)
                if vx is not None : #EV 8/12/21
                    if vz is None: #EV 01/02/19 & 19/07/19
                        arr=sqrt((vx[:,:,1:]/2+vx[:,:,:-1]/2)**2+(vy[:,1:,:]/2+vy[:,:-1,:]/2)**2)
                    else :
                        arr=sqrt((vx[:,:,1:]/2+vx[:,:,:-1]/2)**2+(vy[:,1:,:]/2+vy[:,:-1,:]/2)**2+(vz[1:,:,:]/2+vz[:-1,:,:]/2)**2)
                    #arr = arr[:,-1::-1,:] already done in flofile
                    #arr=arr[0] #EV 14/06/21
                else : return None
        if group=='Transport':
            if name=='Temperature':
                arr = self.core.transReader.readUCN(self.core,'T',tstep,-1,'Tracer');#print shape(arr),arr                
            else:
                arr = self.core.transReader.readUCN(self.core,'Mt3dms',tstep,-1,'Tracer');#print shape(arr),arr
        if group=='Chemistry':
            if name=='Species':
                iesp = self.getNames('Chemistry_Species_L').index(spec)
                arr = self.core.transReader.readUCN(self.core,'Pht3d',tstep,iesp,spec) #iesp=0
            elif name=='User':
                #print(self.userSpecies)
                arr = self.userSpecies[spec][tstep]
        return arr
        
    def getArray2D(self,group,arr3,plane,section):
        '''uses the 3D array to and plane to return a 2D array'''
        X,Y = getXYmeshCenters(self.core,plane,section)# X,Y can be in Z for presention purpose
        self.Umult, units = 1.,self.dicVisu['Chemistry']['Units']
        #print('guish l 206',shape(arr3))
        if group =='Chemistry':
            if self.dicVisu['Chemistry']['Species'] not in ['pH','ph','pe']:
                if units == 'mmol/L': self.Umult = 1000.
                elif units == 'umol/L': self.Umult = 1e6
                elif units == 'nmol/L': self.Umult = 1e9 
        #print('guish 196',shape(arr3),self.mesh,self.core.getValueFromName('Modflow','MshType'))
        if self.MshType>0: # OA 24/10/20 mfUnstruct
            return None,None,arr3[section,:]*self.Umult # OA 7/12/21
        else:
            if self.core.addin.getDim() in ['Radial','Xsection']:
                if self.modgroup[:4]=='Modf': 
                    if self.core.addin.getModelType()=='free' and self.curName=='Head':
                        Y = self.getXyHeadFree(arr3,Y)
                    data = (X,Y,arr3[::-1,0,:]*self.Umult)
                else : 
                    data = (X,Y,arr3[:,0,:]*self.Umult)
            else : # 2 or 3D
                if plane=='Z': data = (X,Y,arr3[section,:,:]*self.Umult); #-1 for different orientation in modflow and real world
                elif plane=='Y': data = (X,Y,arr3[:,section,:]*self.Umult)
                elif plane=='X': data = (X,Y,arr3[:,:,section]*self.Umult)
            #☻print('getA2',shape(X),shape(Y),shape(data[2]))
            return data
                
    def getPointValue(self,x,y):
        """using a coordinate get the value of the current variable at
        this coordinates, using self.data which is a 2D array
        or """
        if x==None or y==None: #OA 24/10/20 removed mesh
            return ' '
        vbox=self.gui.varBox
        var = vbox.choiceV.currentIndex()
        if vbox.chkView.isChecked(): 
            X,Y,mat = vbox.base.getCurVariable(var)
        if self.MshType>0:
            c = self.core.addin.mesh.elcenters
            xc,yc = c[:,0],c[:,1]
            d = sqrt((x-xc)**2+(y-yc)**2);
            inod = where(d==amin(d))[0][0];# print(self.curVarView )
            if self.data == None: zval = 0
            else : zval = self.data[2][inod]
            if vbox.chkView.isChecked(): 
                return '%g, node nb : %i, parm value : %g'%(zval,inod,mat[inod])
            else :
                return '%g, node nb : %i'%(zval,inod)
        else :
            if self.data == None: zval = 0
            else : 
                x0,y0 = self.data[0][0,:],self.data[1][:,0]
                d=x-x0; d1=d[d>0.]
                if len(d1)==0 : return ' '
                ix=where(d==amin(d1))[0][0]
                d=y-y0; d1=d[d>0.]
                if len(d1)==0 : return ' '
                iy=where(d==amin(d1))[0][0] 
                zval = self.data[2][iy,ix]
            if vbox.chkView.isChecked():
                return '%g, parm value %g'%(zval,mat[iy,ix])
            else :
                return '%g'%zval # OA 6/11/18
            
    def getXyHeadFree(self,arr3,Y):
        """modifies the Y coord to show the elevation of the water table
        for cases or radial and xsection. !Y not in the same dir as head"""
        nl,ny,nc = shape(arr3)
        dy=Y[-1,0]-Y[-2,0]
        for ic in range(nc):
            nbl0 = len(where(arr3[:,0,ic]==0)[0]);#print ic,nbl0
            if nbl0>0: 
                Y[nl-nbl0-1,ic] = arr3[nbl0,0,ic]
                Y[nl-nbl0,ic] = Y[nl-nbl0-1,ic]+dy/5
        return Y
        
            
    def getUserSpecies(self,name,plane,layer):
        if 'postfix.phrq' in os.listdir(self.core.fileDir):
            # user species from postfix (pht3d only) adde OA 8/6
            arr3 = self.userSpecies[name][self.Tstep]
            #print shape(arr3),amax(arr3)
        else :
            # classical user species
            formula0 = self.userSpecies[name]
            formula1 = formula0*1
            species = self.getNames('Chemistry_Species_L')
            tsp = str(self.Tstep)
            for e in species:
                if e in formula0: 
                    rp ='self.getArray3D(\'Chemistry\',\'Species\','+tsp+',\''+e+'\')'
                    formula1 = formula1.replace(e,rp)
            arr3 =eval(formula1); # OA  3/10/18
        return self.getArray2D('Chemistry',arr3,plane,layer)        
        
    def get3Dvectors(self,tstep):
        mod = self.core.dicaddin['Model']['group']
        mesh = self.core.addin.mesh
        if mod=="Modflow":
            qx,qy,qz = self.core.flowReader.readFloFile(self.core,tstep);
            qx = qx[:,:,1:]/2+qx[:,:,:-1]/2
            qy = qy[:,1:,:]/2+qy[:,:-1,:]/2
            if qz is None: pass # EV 19/07/19
            else : qz = qz[1:,:,:]/2+qz[:-1,:,:]/2
        elif "USG"in mod:
            qx,qy,qz=self.core.flowReader.readFlowMesh(self.core,mesh,tstep)
        return qx,qy,qz
        
    def getVectors(self,plane,layer,tstep):
        X,Y = getXYmeshCenters(self.core,plane,layer);#print shape(X),shape(Y)
        qx,qy,qz = self.get3Dvectors(tstep)
        if self.core.addin.getDim() =='3D':
            if plane=='Z':
                U=qx[layer]
                V=qy[layer]
        elif self.core.addin.getDim() =='2D':
            U=qx[0]
            V=qy[0]
        elif self.core.addin.getDim() in ['Radial','Xsection']:
            U=qx[:,0,:]
            V=qz[:,0,:]            
        return X,Y,U,V
            
    def getCurrentVariable(self,plane,section):
        mod = self.gui.varBox.parent.currentModel
        line = self.gui.varBox.parent.currentLine
        media, opt, iper = 0,None,0;print('guisho 330',line,self.curVar)
        if line in self.curVar: 
            mat = self.curVar[line]*1
        else:
            mat = self.core.getValueLong(mod,line,0)*1;#print 'show 260',shape(mat)
            self.curVar[line] = mat*1
        X,Y = getXYmeshSides(self.core,plane,section)
        if self.core.addin.getModelGroup()=='Opgeo':
            return None,None,mat[0][0]
        if self.core.addin.getDim() not in ['Radial','Xsection']:
            if plane=='Z': m2 = mat[section,:,:] #-1 for different orientation in modflow and real world
            elif plane=='Y': m2 = mat[:,section,:]
            elif plane=='X': m2 = mat[:,:,section]
        else :
            m2 = mat[-1::-1,0,:]
        self.curVarView = m2
        return X,Y, m2
        
