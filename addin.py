import os, sys, inspect,types
from .config import *
from .Pht3d import *
from .Min3p import *
from .modflowUsg import *
from .geometry import *
#from matplotlib import pylab # for having grafs in batch
from .Pest import *
from .Opgeo import *
from .instantFit import *
from .modflowWriter import *
from .multiPlot import * # oa 28/11/18

#from matplotlib.tri import CubicTriInterpolator
from PyQt5.QtGui import *
from PyQt5 import Qt
from PyQt5.QtWidgets import *
def onMessage(gui,text):  QMessageBox.information(gui,"Info",text)

class addin:
    """the addin class is used to add buttons or menus inmesh the interface to
        manage things that are around the model disctionnaries
        the structure dict stores the name and location of addins
        the data are stored in dicaddin (a copy of core.dicaddin)"""
    def __init__(self,core,gui=None):
        self.gui,self.core = gui,core
        self.dickword = core.dickword
        self.grd = None
        self.initAddin()
        
    def setGui(self,gui):
        self.gui = gui
        cfg = Config(self.core)
        if gui != None:
            self.dialogs = cfg.dialogs
            self.gtyp = cfg.gtyp
        
    def initAddin(self):
        self.structure={'button':{},'menu':{}}
        self.pht3d = PHT3D(self.core)
        self.min3p = Min3p(self.core)
        self.mfU = modflowUsg(self.core)
        self.opgeo = Opgeo(self.core)
        self.pest = Pest(self.core)
        # creating usedModules addin
        self.lastBatch = ''
        name='usedM'
        self.core.dicaddin[name] = {}
        for mod in self.core.modelList:
            lmodules = self.dickword[mod].grpList;#print(lmodules)
            if mod[:4] in ['Min3','Opge','Sutr','Pest']: val = [True]*len(lmodules) #,'Fipy'
            else : val = [False]*len(lmodules)
            for i in range(len(lmodules)):
                if (mod=='Modflow') & (lmodules[i] in ['DIS','BAS6','LPF','WEL']): val[i]=True # OA 10/3/19 removed pcg
                if (mod=='Mt3dms') & (lmodules[i] in ['BTN','ADV','DSP','GCG']): val[i]=True
                if (mod=='MfUsgTrans') & (lmodules[i] in ['BCT','PCB','CRCH','CWELL']): val[i]=True # OA 27/7/19
                if (mod=='Pht3d') & (lmodules[i] in ['PH']): val[i]=True
                if (mod=='Observation') & (lmodules[i] in ['OBS']): val[i]=True
            if mod=='Modflow': val = self.addSolver(lmodules,val) # OA 10/3/19
            self.structure['menu'][mod]={'name':name,'position': 0,
                'function': 'onUsedModules','short':'M'}
            self.core.dicaddin[name+'_'+mod] = (lmodules,val) # only data are stored in dicaddin  
        # creating dic for observation data
        name='obsHead'
        self.core.dicaddin[name] = {}
        self.structure['menu'][mod]={'name':name,'position': 0,
                'function': 'OnImportHead','short':'oH'}
        name='obsTracer'
        self.core.dicaddin[name] = {}
        self.structure['menu'][mod]={'name':name,'position': 0,
                'function': 'OnImportTracer','short':'oH'}
        name='obsChemistry'
        self.core.dicaddin[name] = {}
        self.structure['menu'][mod]={'name':name,'position': 0,
                'function': 'OnImportChemistry','short':'oH'}
        # creating the structure for the buttons
        name = 'Model'
        model = {'dimension':'2D horizontal','type':'Confined','group':'Modflow series'} #EV 22/07/2018 2D -> 2D horizontal
        self.core.dicaddin[name] = model
        self.structure['button']['1.Model'] = [{'name':name,'pos':0,'short':'M'}]
        self.checkFlowMod=False # EV 27/04/20
        name = 'Grid'
        grid = {'x0':'0','y0':'0','x1':'100','y1':'50','dx':'5','dy':'2'} # default grid
        self.core.dicaddin[name] = grid
        self.structure['button']['1.Model'].append({'name':name,'pos':0,'short':'G'})
        self.checkDomain = False # do not consider the diczone in dis.1 as domain limit OA 14/11/18
        name = '3D'
        dim = {'zmin':0.,'topMedia':[10],'listLayers':[1]}
        self.core.dicaddin[name] = dim
        self.structure['button']['1.Model'].append({'name':name,'pos':0,'short':'3D'})
        name = 'Time' # select times for the model
        self.core.dicaddin[name] = {'final':'10.','steps':'1.'} # EV 18/02/19 remove 'mode':'linear'} # default time
        self.structure['button']['1.Model'].append({'name':name,'pos':0,'short':'T'})

        name = 'Particle' # create particles
        self.structure['button']['2.Flow'] = [{'name':name,'pos':1,'short':'P'}]
        
        name = 'MtSpecies' # llow to select species for mt3d
        self.core.dicaddin[name] = {'flag':'Mt3dms','species':[]}
        self.structure['button']['3.Transport'] = [{'name':name,'pos':0,'short':'Spc'}]
        name = 'MtReact' # values of parameters for each species
        self.core.dicaddin[name] = {}
        self.structure['button']['3.Transport'].append({'name':name,'pos':0,'short':'Rct'})
        
        name = 'ImpDb' # import the phreeqc database (no data needed)
        self.structure['button']['4.Chemistry']=[{'name':name,'pos':0,'short':'I'}]
        name = 'Chemistry' # opens the chemistry dialog
        self.core.dicaddin[name] = {}
        self.core.dicaddin['MChemistry'] = {} # for Min3p to store a different chemistry
        self.structure['button']['4.Chemistry'].append({'name':name,'pos':0,'short':'C'})        

        name = 'Pback' # dict for the pest zoens parameters
        self.core.dicaddin['Pback1'] = {} # to choose the parameters
        self.core.dicaddin['Pback2'] = {} # to provide values
        self.structure['button']['5.Pest']=[{'name':name,'pos':0,'short':'Pbk'}]

        name = 'Pzones' # dict for the pest zoens parameters
        self.core.dicaddin['Pzones1'] = {}
        self.core.dicaddin['Pzones2'] = {}
        self.structure['button']['5.Pest'].append({'name':name,'pos':0,'short':'Pz'})
        
        name = 'Pestchek' # run Pestchek EV 07/11
        self.structure['button']['5.Pest'].append({'name':name,'pos':0,'short':'Pz'})

        name = 'InitialChemistry' # to set specific initial chemistry
        self.core.dicaddin[name] = {'name':'','formula':'','tstep':''} #ev05/02/19 removed value =
        self.grd = makeGrid(self.core,self.core.dicaddin['Grid'])
        self.mesh = None  # OA added 25.9.18
        self.setChemType()
        self.setMfUnstruct()
        self.fit = instant(self.gui,self.core)
        self.particle = {}
        
    def addSolver(self,lmodules,val): # OA adde 10/3/19 to select one solver, at least pcg
        a = False
        for i in range(len(val)): 
            if (lmodules[i] in ['PCG','NWT','SIP','SOR','DE4']) and (val[i]) : a=True
        if a==False : val[lmodules.index('PCG')] = True
        return val
        
    def setMfUnstruct(self): #OA added 17/9/17 for modflow UNS
        lmodules,val = self.core.dicaddin['usedM_Modflow']
        bool = self.core.mfUnstruct;#print('in mfuns',bool,lmodules)
        val[lmodules.index('DIS')] = not bool
        val[lmodules.index('PCG')] = not bool
        val[lmodules.index('DISU')] = bool
        val[lmodules.index('SMS')] = bool
        self.core.dicaddin['UsedM_Modflow'] = (lmodules,val)
        
    def setChemType(self):
        self.chem = self.pht3d
        if self.core.dicaddin['Model']['group'] == 'Min3p': 
            self.chem = self.min3p
        self.chem.resetBase()
        
    def update1(self,dict1):
        cdict = self.core.dicaddin.copy()
        for k in list(dict1.keys()):
            if k not in list(cdict.keys()): continue
            if type(dict1[k])==type({'d':1}):
                for k1 in list(dict1[k].keys()):
                    cdict[k][k1] = dict1[k][k1]
            # elif k[:5]=='usedM': # some moduls can have been added in new version # OA removed 10/3/19
            #     lmod,val = cdict[k]
            #     lmod1,val1 = dict1[k]
            #     for i,md in enumerate(lmod): 
            #         if md in lmod1: 
            #             if val1[lmod1.index(md)]: val[i]=True
            #     cdict[k] = (lmod,val)
            else :
                cdict[k] = dict1[k]
        self.core.dicaddin = cdict.copy()
                
    def initMenus(self):
        '''add the menus in the gui interface''' 
        self.gui.addMenu(502,'Modflow_modules',self.onUsedModules)   
        self.gui.addMenu(503,'Mt3dms_modules',self.onUsedModules)   
        self.gui.addMenu(504,'MfUsgTrans_modules',self.onUsedModules)  #OA 27/5/21 
        #self.gui.addMenu(504,'Interactive fitting',self.onInstantFit)   
        #self.gui.addMenu(505,'MultiPlot',self.onMultiPlot)   #OA 28/11/18
        self.gui.addMenu(506,'Batch',self.onBatchDialog)   
        self.gui.addMenu(507,'Initial chemistry',self.onInitialChemistryDialog)   
        
    def addButton(self, location, panelName):
        """this method is called by an external panel (for instance parameters)
        and provides a button that will open a dialog able to modify the content
        of the addin dictionnary
        the action to be done is stored in the parametersGui dict of actions """
        sButt = self.structure['button']
        but = None
        if panelName in sButt:
            but = []
            for bdic in sButt[panelName]:
                name = 'Ad_'+bdic['name']
                but.append((bdic['short'],name,bdic['pos']))
        return but
        
    def doaction(self,actionName):
        """the action to be done when an addin button is clicked"""
        if actionName == 'Ad_Model':
            m = self.core.dicaddin['Model']
            if m['dimension']=='2D': # EV 22/07/2019 to transform 2D in 2D horizontal & free in Unconfined & confined in Confined in the existing model
                m['dimension']='2D horizontal'
            if m['type']=='free':
                m['type']='Unconfined' 
            if m['type']=='confined':
                m['type']='Confined'
            data = [('Dimension','Choice',(m['dimension'],['2D horizontal','3D','Radial','Xsection'])), #EV 22/07/2018 2D -> 2D horizontal
                    ('Type','Choice',(m['type'],['Confined','Unconfined','Mix (for 3D model)'])),  #EV 22/07/2018 free -> Unconfined
                    ('Group','Choice',(m['group'],['Modflow series','Modflow USG','Min3p','Opgeo'])),#EV 18.10.18 removed 'Sutra'
                    ('Use other flow model?','Check',self.checkFlowMod)] # EV 27/04/20
            dialg = self.dialogs.genericDialog(self.gui,'Model',data)
            retour = dialg.getValues()
            nmed=getNmedia(self.core) #EV 25/09/19
            if retour != None:
                self.gui.onGridMesh('Grid') # default grid button
                m['dimension'],m['type'],m['group'],self.checkFlowMod = retour # EV 27/04/20
                self.core.mfUnstruct = False;self.setMfUnstruct() # OA 1/3/20 added cond
                if m['group'] == 'Modflow USG':
                    self.core.mfUnstruct = True
                    if self.core.getValueFromName('Modflow','MshType')>0: self.gui.onGridMesh('Mesh') # to chang the button
                    self.mesh = self.mfU # OA 22/8/19
                    self.setMfUnstruct()
                self.gui.varBox.chooseCategory(m['group'])
                if m['type']=='Confined': #EV 25/09/19
                    self.core.setValueFromName('Modflow','LAYTYP',[0]*nmed) # 0 for confined, 1 for unconfined
                    self.core.setValueFromName('Mt3dms','TLAYCON',[0]*nmed)
                if m['type']=='Unconfined': #EV 22/07/2018 free -> Unconfined
                    self.core.setValueFromName('Modflow','LAYTYP',[1]*nmed) # 0 for confined, 1 for unconfined
                    self.core.setValueFromName('Mt3dms','TLAYCON',[1]*nmed)
                if m['dimension'] in ['Xsection','Radial']:
                    self.core.setValueFromName('Modflow','TOP',1.)
                    self.core.setValueFromName('Modflow','BOTM',0.)
                    self.core.setValueFromName('Modflow','DELC',1.)
                    self.core.setValueFromName('Modflow','NROW',1)
                    self.gui.on3D(False) #EV 26/09/19
                self.gui.onWriteModflow(self.checkFlowMod) # EV 27/04/20
                self.set3D()
                self.setChemType()
                self.gui.on3D(m['dimension']=='3D')  # OA 14/3/21
                self.gui.onSetMediaNb(nmed,getNlayers(self.core))  # OA 14/3/21
                
        if actionName =='Ad_Grid':
            g = self.core.dicaddin['Grid']
            x0,x1,y0,y1 = g['x0'],g['x1'],g['y0'],g['y1'] # OA 6/11/18
            dicz = self.core.diczone['Modflow'].dic # OA 6/11/18
            if ('dis.1' in dicz.keys()) and (self.gtyp=='qgis'): # OA 6/11/18, modif 20/11/19
                data = [('Use the dis.1 domain zone?','Check',self.checkDomain)]
                dlg1 = self.dialogs.genericDialog(self.gui,'domain zone',data) # OA 14/11
                retour = dlg1.getValues() # OA 14/11
                if retour[0] : # dilg returns true OA 14/11
                    self.checkDomain = True
                    x,y = zip(*dicz['dis.1']['coords'][0]) # OA 6/11/18
                    x0,x1,y0,y1 = min(x),max(x),min(y),max(y) # OA 6/11/18
            dvert = 'dy';vert='Y'
            if self.getDim() in ['Radial','Xsection']:
                dvert = 'dz';vert = 'Z'
            data = [('Xmin','Text',x0),(vert+'min','Text',y0),
                    ('Xmax','Text',x1),(vert+'max','Text',y1),
                    ('dx','Textlong',g['dx']),(dvert,'Textlong',g['dy'])]
            if self.mesh == None or self.core.getValueFromName('Modflow','MshType')==0: #OA 21/3/20
                dialg = self.dialogs.genericDialog(self.gui,'Grid',data)
                retour = dialg.getValues()
            else : #mesh case no dialog
                self.setGridInModel('new')
                self.gui.visu.initDomain() 
                self.gui.onGridMesh('Mesh') # to chang the button # EV 27/04/20
                retour = None               
            if retour != None:
                g['x0'],g['y0'],g['x1'],g['y1'],g['dx'],g['dy'] = retour;#print g
                if self.checkDomain:
                    x0,x1,y0,y1 = float(g['x0']),float(g['x1']),float(g['y0']),float(g['y1']) # OA 14/11/18
                    dicz['dis.1']['coords'][0] = list(zip([x0,x1,x1,x0,x0],[y1,y1,y0,y0,y1])) # OA 14/11/18 set dis.1 domain new bdy
                self.setGridInModel('new')
                self.gui.visu.initDomain()

        if actionName == 'Ad_3D':
            m = self.core.dicaddin['3D']
            data = [('Top of Media','Textlong',m['topMedia']),
                    ('Bottom','Text',m['zmin']),
                    ('Nb of layers','Textlong',m['listLayers'])
                    ]
            dialg = self.dialogs.genericDialog(self.gui,'3D',data)
            retour = dialg.getValues()
            if retour != None:
                m['topMedia'],m['zmin'],m['listLayers'] = retour
                self.set3D()
                dic = self.make3DTable()
                dialg = self.dialogs.myNoteBook(self.gui,'3D',dic)
                retour = dialg.getValues()
            self.gui.onSetMediaNb(getNmedia(self.core),getNlayers(self.core))  # OA 14/3/21
                
        if actionName == 'Ad_Time':
            t = self.core.dicaddin['Time']
            data = [('Total simulation time','Textlong',t['final']),
                    ('Step size','Textlong',t['steps'])] # EV 18/02/19 remove 'mode'
                    #('Step mode','Choice',(t['mode'],['linear','log']))]
            dialg = self.dialogs.genericDialog(self.gui,'Time',data)
            retour = dialg.getValues()
            if retour != None:
                t['final'],t['steps']= retour # EV 18/02/19 removed t['mode'] 
                self.setTime()
                
        if actionName == 'Ad_Particle':
            self.particle={'direction':1,'type':'transient'} # forward
            data = [('direction','Choice',('forward',['forward','backward'])),
                    ('type','Choice',('transient',['transient','steady']))]
            dialg = self.dialogs.genericDialog(self.gui,'Particle',data)
            retour = dialg.getValues()
            if retour != None: 
                if retour[0] == 'backward' : self.particle['direction'] = -1
                self.particle['type'] = retour[1]
                #self.gui.actions('zoneStart')
                self.gui.visu.startParticles()
                self.gui.guiShow.dlgShow.onTickBox('Flow','Particles','B',True)
            else : return
                
        if actionName == 'Ad_MtSpecies':
            m = self.core.dicaddin['MtSpecies']
            data = [('Type','Choice',(m['flag'],['Mt3dms','Pht3d'])),
                    ('Species','Textlong',m['species'])
                    ]
            dialg = self.dialogs.genericDialog(self.gui,'Select Species',data)
            retour = dialg.getValues()
            if retour != None:
                m['flag'],m['species'] = retour
                #if m['flag']=='Pht3d': m['species'] = self.pht3d.getListSpecies()
            self.setMtSpecies(m['flag'],m['species'])
                
        if actionName == 'Ad_MtReact':
            dic = {'Parameters':self.core.dicaddin['MtReact'].copy()}
            dialg = self.dialogs.myNoteBook(self.gui,"Rct parameters",dic)
            dic2 = dialg.getValues()
            if dic2 != None:
                self.core.dicaddin['MtReact'] = dic2['Parameters']
            #dialg.Destroy()  

        if actionName == 'Ad_ImpDb':
            if self.core.dicaddin['Model']['group'] == 'Min3p':
                self.chem = self.min3p
                self.chem.importDB(self.min3p.Base['MChemistry'])
                self.callCheckDialog()
                self.chem.updateChemistry()
                bs = self.chem.Base['MChemistry'];#print bs
                self.chem.readKinetics(bs['redox']['rows'])
                self.chem.readMinerals(bs['mineral']['rows'])
#                for name in ['comp','complex','gases']:
#                    self.chem.readOtherDb(name,bs[name]['rows'])
#                print bs
                self.dialogs.onMessage(self.gui,'Database imported')
            else :
                self.chem = self.pht3d
                if 'pht3d_datab.dat' not in os.listdir(self.core.fileDir):
                    self.dialogs.onMessage(self.gui,'pht3d_datab.dat file missing')
                else : #EV 26/08/19
                    fname = str(self.core.fileDir+os.sep+'pht3d_datab.dat')
                    self.pht3d.tempDbase,self.pht3d.npk= self.pht3d.importDB(fname);
                    self.chem.updateChemistry()
                    self.dialogs.onMessage(self.gui,'Database imported')
            
        if actionName == 'Ad_Chemistry':
            if self.core.dicaddin['Model']['group'] == 'Min3p': 
                nameB='MChemistry'
                dic = self.chem.getBase()
            else :
                nameB = 'Chemistry'
                dic = self.chem.Base[nameB].copy()
            dialg = self.dialogs.myNoteBook(self.gui,"Chemistry",dic)
            dic2 = dialg.getValues()
            if dic2 != None:
                for k in list(dic2.keys()): self.chem.Base[nameB][k] = dic2[k]
                if nameB == 'MChemistry': # OA 1/3/19 for exchange species
                    self.chem.Base[nameB]['exchange']['text'] = dic['exchange']['text'] 
            self.core.dicaddin[nameB] = self.chem.Base
        
        if actionName == 'Ad_Pback':
            dic = self.pest.getDicBack1() # choose the line to modify
            dialg = self.dialogs.myNoteBookCheck(self.gui,"Pest background",dic)
            dic2 = dialg.getValues();#print 'addin 304',dic2
            if dic2 != None:
                self.core.dicaddin['Pback1'] = dic2
            else : return
            dic = self.pest.getDicBack2() # provide the values for the selected lines
            dialg = self.dialogs.myNoteBook(self.gui,"Pest background",dic)
            dic2 = dialg.getValues()
            if dic2 != None:
                self.core.dicaddin['Pback2'] = dic2
                self.setPestParm()

        if actionName == 'Ad_Pzones':
            dic = self.pest.getDicZones1() # choose the line to modify
            dialg = self.dialogs.myNoteBookCheck(self.gui,"Pest zones",dic)
            dic2 = dialg.getValues();#print 'addin 304',dic2
            if dic2 != None:
                self.core.dicaddin['Pzones1'] = dic2
            else : return
            dic = self.pest.getDicZones2()
            dialg = self.dialogs.myNoteBook(self.gui,"Pest zones",dic)
            dic2 = dialg.getValues()
            if dic2 != None:
                self.core.dicaddin['Pzones2'] = dic2
                self.setPestParm()
        
        if actionName == 'Ad_Pestchek':
            self.onPestchek()

    def getUsedModulesList(self,modName):
        """returns only the modules that are used as a list"""
        #print self.core.dicaddin['usedM_'+modName]
        modules,val=self.core.dicaddin['usedM_'+modName]
        l0=[]
        for i in range(len(modules)): 
            if val[i]: l0.append(modules[i])
        return l0
        
    def setUsedModulesList(self,modName,mlist):
        modules,val=self.core.dicaddin['usedM_'+modName]
        for i in range(len(modules)):
            if modules[i] in mlist: val[i]=True
            else : val[i]=False
    
    def onUsedModules(self):
        if self.gtyp=='qt': txt = self.gui.menuBar().sender().text()
        else : txt = self.gui.file.File.sender().currentText()
        modName = str(txt).split('_')[0]
        data = self.core.dicaddin['usedM_'+modName]
        # below not really elegant but to have mnwt, rct and vdf for older version
        for n in ['MNWT']:
            if (modName=='Modflow') & (n not in data[0]): 
                data[0].append(n);data[1].append(False)
        for n in ['SSMs','RCT','VDF','UZT']:
            if (modName=='Mt3dms') & (n not in data[0]): 
                data[0].append(n);data[1].append(False)
        typ = ['Check']*len(data[0])
        dialg = self.dialogs.genericDialog(self.gui,'Select Modules',list(zip(data[0],typ,data[1])))
        chkStates = dialg.getValues()
        if chkStates != None:
            self.core.dicaddin['usedM_'+modName] = [data[0],chkStates]
            item = self.gui.varBox.choiceG
            l0 = []
            for i in range(len(data[0])):
                if data[1][i]: l0.append(data[0][i])
            self.gui.varBox.setChoiceList(item,l0)
        mtm,mtval = self.core.dicaddin['usedM_Mt3dms']
        if 'DIS' in mtm: bDis = mtval[mtm.index('DIS')] # OA 02/20
        bool = False
        if 'RCT' in mtm : bool = mtval[mtm.index('RCT')]
        self.gui.onRCT(bool)
        
    def callCheckDialog(self):
        dicIn = self.chem.temp['Dbase'].copy()
        cpl = self.chem.temp['Dbase']['complex']*1
        dicIn.pop('complex')
        #print 'incallcheck',dicIn
        dialg = self.dialogs.myNoteBookCheck(self.core.gui,'Choose species',dicIn,'sort') # OA 3/4
        retour = dialg.getValues()
        if retour!= None:
            self.chem.temp['Dbase'] = retour
        self.chem.temp['Dbase']['complex']=cpl
        
    def onInstantFit(self,evt):
        '''creates the box for instant fitting and starts the observer of changes'''
        #create the object that observe the cnage in topbar
        self.fit.setObserver(self.gui,self.gui.modifBox.obs)
        self.fit.startDialog()

    def onBatchDialog(self,evt):
        #from matplotlib import rcParams  # OA 25/1/21
        #rcParams['interactive']=True  # OA 25/1/21
        head='insert python commands below'
        dialg = self.dialogs.textDialog(self.gui,'Batch program',(500,300),self.lastBatch)
        retour = dialg.getText();#print retour
        if retour != None:
            txt = retour #dialg.GetTextAsList()
            self.lastBatch=txt
            txt1=txt.replace('core','self.core');#print(type(txt1))
            self.formExec(txt1) # OA 25/10/18
        else : return
        
    def formExec(self,s): # added OA 25/10/18
        s1 = s.replace('\n','\n\t')
        s1 = 'def yoyo(self):\n\t'+s1
        dct={}
        exec(s1,globals(),dct)
        b = types.MethodType(dct['yoyo'], self)
        b()

    def onInitialChemistryDialog(self,evt):
        self.head='calculate initial chemistry'
        inC = self.core.dicaddin['InitialChemistry']
        self.txt = 'name: '+inC['name'] +'\n'+'formula: '+inC['formula'] +'\n'+'tstep: '+inC['tstep'] # EV 05/02/19
        #self.txt+= '\n#name:Ca, formula: value=ones((25,30))*5e-4' # EV 05/02/19
        self.txt+= '\n#name:All, formula: value=importUCN, tstep:0'
        #self.txt+= '\n#name:Hfo_w, formula: value=core.importLayerValues(\'Hfo_layers.txt\',\'Hfo_w\')'
        #print(self.txt)
        dialg = self.dialogs.textDialog(self.gui,'Initial chemistry',(500,300),self.txt)
        retour = dialg.getText();print ('ok',retour)
        if retour != '':
            f0 = retour.split('#')[0] # EV 05/02/19
            name=f0.split('formula: ')[0].split(':',1)[1].strip()
            formula=f0.split('formula: ')[1].split('tstep')[0].strip()
            tstep=f0.split('formula: ')[1].split('tstep:')[1].strip()
            self.core.dicaddin['InitialChemistry']={'name':name,'formula':formula,'tstep':tstep}
        else : 
            self.core.dicaddin['InitialChemistry']={'name':'','formula':'','tstep':''}
            return
        #dialg.Destroy()
        
    def setGrd(self,grd): self.grd = grd
    
    def getFullGrid(self): return self.grd      
        
    def setGridInModel(self,opt='old'):
        '''when the grid dialog has been filled transmit info to the models'''
        g0 = self.core.dicaddin['Grid']
        g = makeGrid(self.core,g0)
        self.grd = g
        self.mesh, self.xsect = None, False
        mgroup = self.core.dicaddin['Model']['group'];#print 'addin line 316 mgroup', mgroup
        if self.getDim() in ['Radial','Xsection']: self.xsect = True # oa 26/5
        if mgroup in ['Modflow series']: # oa 9/2/20 added rect, 29/2 removed
            self.core.dicval['Modflow']['dis.4']=list(g['dx']) # pb with setting list
            self.core.dicval['Modflow']['dis.5']=list(g['dy'])
            self.core.setValueFromName('Modflow','NCOL',g['nx'])
            self.core.setValueFromName('Modflow','NROW',g['ny'])
        elif mgroup =='Modflow USG':
            self.core.lcellInterp = []  # OA 16/1/21
            self.mfU.buildMesh(opt)
            ncell = self.mfU.ncell
            self.core.setValueFromName('Modflow','NCELL',int(ncell))
            self.core.setValueFromName('Modflow','NODELAY',self.mfU.ncell_lay)
            self.core.setValueFromName('Modflow','NJAG',self.mfU.nconnect+ncell)
        elif mgroup =='Min3p' :
            l2,l3 = 'spat.2','spat.3'
            if self.xsect: l2,l3 = 'spat.3','spat.2'
            self.core.dicval['Min3pFlow']['spat.1']=[1,g['nx'],g['x0'],g['x1']] 
            self.core.dicval['Min3pFlow'][l2]=[1,g['ny'],g['y0'],g['y1']] 
            self.core.dicval['Min3pFlow'][l3]=[1,1,0.0,1.0] 
            self.min3p.buildMesh()
        elif mgroup =='Sutra' :
            nx,ny = int(g['nx']),int(g['ny'])
            if self.getDim() == '2D': nbL,nbL1 = 1,1
            else : 
                nbL = getNlayers(self.core)
                nbL1 = nbL + 1
            nno,nel = (nx+1)*(ny+1)*nbL1, nx*ny*nbL
            self.core.dicval['Sutra']['glob.2b'][2:] =[nx+1,ny+1,nbL1] 
            self.core.dicval['Sutra']['glob.3'][:2] = [nno,nel]
        elif mgroup =='Opgeo' : # OA 26/3/19 Added 5 next lines
            nz = getNlayers(self.core)
            nx,ny = int(g['nx']),int(g['ny'])
            if self.xsect : 
                nz = ny*1; ny = 1
            self.core.dicval['OpgeoFlow']['domn.1'][1:4] = [nz,nx,ny]
            self.core.dicval['OpgeoFlow']['domn.2']=list(g['dx']) 
            self.core.dicval['OpgeoFlow']['domn.3']=list(g['dy'])   
            self.opgeo.buildMesh()
        
    def getModelType(self): return self.core.dicaddin['Model']['type']
    def getModelGroup(self): return self.core.dicaddin['Model']['group']
        
    def getDim(self): 
        dm = self.core.dicaddin['Model']['dimension']
        if dm == '2D horizontal' : dm = '2D' #EV 22/07/2018 2D -> 2D horizontal
        #if dm=='Cross_Section': dm= 'Xsection' # OA 9/6
        return dm
        
    def get3D(self): return self.core.dicaddin['3D']
    
    def set3D(self):
        dm = self.getDim()
        if dm not in ['2D','3D'] : return
        med = self.core.dicaddin['3D']
        nbL = getNlayers(self.core)
        z0 = med['zmin']
        botM = med['topMedia'][1:];botM.append(z0)
        if self.getModelGroup()=='Modflow series': 
            self.core.dicval['Modflow']['dis.6'] = med['topMedia']
            self.core.dicval['Modflow']['dis.7'] = botM
        if self.getModelGroup()=='Modflow USG': 
            self.core.dicval['Modflow']['disu.7'] = med['topMedia']
            self.core.dicval['Modflow']['disu.8'] = botM
        if self.getModelGroup()=='Opgeo': 
            self.core.setValueFromName('OpgeoFlow','O_TOP',med['topMedia'])
            self.core.setValueFromName('OpgeoFlow','O_BOTM',botM)
        self.core.setValueFromName('Modflow','NLAY',nbL)            
        self.core.setValueFromName('Modflow','UNLAY',nbL)            
            
    def make3DTable(self):
        nbL,lilay,dzL = makeLayers(self.core)
        toplist = [float(a) for a in self.core.addin.get3D()['topMedia']]
        nbM = len(lilay)
        bot = float(self.core.addin.get3D()['zmin'])
        toplist.append(bot);#print toplist
        dic={'cols':['Media','Layer','z'],'rows':[str(a) for a in range(nbL)],'data': []}
        top,nl = toplist[0],0
        for im in range(nbM):
            ep = toplist[im]-toplist[im+1]
            for il in range(lilay[im]):
                dz = dzL[im][il]*ep
                dic['data'].append([im,nl,nice(top)+' to '+nice(top-dz)]) # OA 26/8/19 passed to nice format
                top -= dz
                nl += 1
        return {'3Dlayers':dic}
        
    def setTime(self):
        a = self.core.makeTtable()
        tlist = self.core.getTlist2();nper = len(tlist)
        self.gui.guiShow.setNames('Model_Tstep_L',tlist,'numbers')
        mgroup = self.getModelGroup()
        if mgroup =='Modflow series': # OA 3/10 f was missing
            self.core.setValueFromName('Modflow','NPER',nper)
        elif mgroup =='Modflow USG':
            self.core.setValueFromName('Modflow','UNPER',nper)
        elif mgroup [:3]=='Min':
            self.core.setValueFromName('Min3pFlow','Tfinal',tlist[-1])
            self.core.setValueFromName('Min3pFlow','Tmaxstep',tlist[-1]/100.)
        elif mgroup [:3]=='Sut':        
            tl2 = ' '.join([str(t) for t in tlist])
            self.core.setValueFromName('Sutra','TLIST',tl2)
            self.core.setValueFromName('Sutra','NTLIST',nper)
            stepdiv = self.core.getValueFromName('Sutra','STEPDIV')
            self.core.setValueFromName('Sutra','NCOLPR',stepdiv)
            self.core.setValueFromName('Sutra','LCOLPR',stepdiv)

    def setMtSpecies(self,flag,species):
        rows = species
        #nrows = len(rows)
        cols = ['SP1','SP2','RC1','RC2']
        data,rowIn,dataIn = [],[],[]
        mtreact = self.core.dicaddin['MtReact']
        if 'rows' in mtreact: 
            rowIn,dataIn = mtreact['rows'],mtreact['data']
        for sp in rows :
            if sp in rowIn: 
                data.append(dataIn[rowIn.index(sp)])
            else : 
                data.append([0]*len(cols))
        self.core.dicaddin['MtReact']={'rows':rows,'cols':cols,'data':data}
        
    def getMtSpecies(self): return self.core.dicaddin['MtSpecies']
    
    def getMtReact(self): return self.core.dicaddin['MtReact']
    
    def setPestParm(self): #EV 06/11
        self.pgrp = [] #list of parameter group
        self.dicPback= self.core.dicaddin['Pback2']
        self.dicPzones=self.core.dicaddin['Pzones2']
        for md in list(self.dicPback.keys()):
            for i,line in enumerate(self.dicPback[md]['rows']):
                dp = self.dicPback[md]['data'][i]
                if dp[0]:self.pgrp.append(dp[6]) #EV 25/07/19
        for line in list(self.dicPzones.keys()):
            cols = self.dicPzones[line]['cols']
            i1 = cols.index('Use')
            for i,zname in enumerate(self.dicPzones[line]['rows']):
                dp = self.dicPzones[line]['data'][i]
                if dp[i1]:self.pgrp.append(dp[6])
        self.nparmgrp = len(unique(self.pgrp))
        self.core.dicval['Pest']['pgr.1'][0]= self.nparmgrp 
        for i in range(self.nparmgrp):
            self.core.dicval['Pest']['pgr.'+str(i+2)][0]=str(unique(self.pgrp)[i])
            
    def onPestchek(self): #EV 07/11
        pest=Pest(self.core)
        message=pest.Pestchek()
        self.dialogs.onMessage(self.gui,message)
            
###################### CALCULATION OF PARTICLES #################
    def calcParticle(self,xp0,yp0,zoneMatrix=None):
        """ represent the particle in the visu from x0,y0 coordinates"""
        grp = self.getModelGroup()
        iper=self.gui.guiShow.Tstep
        ilay = self.gui.visu.curLayer
        try : startTime = float(self.gui.guiShow.getCurrentTime())
        except ValueError : startTime = 0.001
        if grp[:2]=='Op' and self.core.dicval[grp+'Flow']['domn.1'][0]>0:
            xp,yp,tp = self.calcPartMesh(xp0,yp0)
        else :
            grd  = self.getFullGrid()
            dx,dy = array(grd['dx']), array(grd['dy']);
            data = array([grd['x0'],grd['y0'],xp0,yp0])
            #[xp,yp,tp,cu1,nt] = iphtC1.calcPart1(data,dx,dy,vx[0],vy[0],clim)
            [xp,yp,tp,cu1,nt] = self.calcPartGrid(data,dx,dy,startTime,iper,ilay,zoneMatrix) #,clim)
            xp,yp,tp = xp[:nt],yp[:nt],tp[:nt]
            xp,yp,tp = xp[tp>=0],yp[tp>=0],tp[tp>=0] # OA added 14/10/19
        self.particle['data'] = [xp,yp,tp]
        
    def calcPartGrid(self,datain,dxi,dyi,startTime,iper,ilay,zoneMatrix=None):
        """transient particle tracking for 2D, xsection and semi-3D (along layer)
        fully validated in 2D, except that with pollock there are non continuity
        at the cell boundaries (velocity field is not fully interpolated )"""
        eps=1.e-6;cst1=5.e-5;#print self.particle
        infos = self.core.flowReader.getPart()
        poro = self.core.getValueFromName('Mt3dms','PRSTY')
        tlist = self.core.getTlist2()
        #print startTime, iper, ilay
        thick = self.core.flowReader.getThickness(self.core,[iper])[0];# can be 3D
        x0,y0,xp0,yp0 = datain;
        #ny,a=shape(vxi);nx=a-1;
        nx,ny = len(dxi),len(dyi)
        nt = int(max(nx,ny)*3.);#print 'addin 564',datain,nt # mod nx
        xg, yg = r_[x0,x0+cumsum(dxi)],r_[y0,y0+cumsum(dyi)];#print 'addin 568',xg,yg
        #lx,ly = arange(nx+1),arange(ny+1)

        #* start calculations */
        xp=zeros((nt+1));yp=xp*0.;tp=xp*0.;cu=xp*0.;
        xp1,yp1,cu1,tp1 = [],[],[],[]
        it = 0;xp[0]=xp0;yp[0]=yp0;tp[0]=startTime;ptin=1;
        if zoneMatrix==None: zoneMatrix = ones((ny,nx))
        #jp = lx[xp0>xg][-1]; 
        #ip = ly[yp0>yg][-1]; 
        while (it<nt) and (ptin==1): # and (tp[it]>0):
            it+=1
            dxc = 0.; dyc = 0.;dt = 0.; #dx=dxi[jp]; dy=dyi[ip]; print ip,jp
            if tp[it-1]>tlist[-1] : # OA 20/11/19 changed iper+1
                if self.particle['type'] == 'steady': iper = 0 #0 transietn, 1 steady
                elif tp[it-1]<tlist[-1]:  
                    a=tlist-tp[it-1];iper=min(where(a>0)[0])-2
                    if iper<0: break
                else : break
            #print 'addin 577',it,xp[it-1],yp[it-1],tp[it-1],iper
            #if (xp[it-1]<xg[-2]) & (xp[it-1]>xg[1]): jp = lx[xp[it-1]>xg][-1]
            #else : ptin = 0
            #if (yp[it-1]<yg[-2]) & (yp[it-1]>yg[1]): ip = ly[yp[it-1]>yg][-1]
            #else : ptin = 0
            #dx=dxi[jp];dy=dyi[ip];#print it,ip,jp,xp[it-1],yp[it-1]
            dist=xg[nx-1]-x0;jp=0;
            for j in range(nx):
                a=xp[it-1]-xg[j]
                if ((a<dist)and(a>=0.)):dist=a;jp=j
            jp=clip(jp,0,nx-1);dx=dxi[jp];

            dist=yg[ny-1]-x0;ip=0;
            for i in range(ny):
                a=yp[it-1]-yg[i];
                if ((a<dist)and(a>=0.)):dist=a;ip=i
            ip=clip(ip,0,ny-1);dy=dyi[ip];
            
            cl = zoneMatrix[ip,jp] #clim[ip,jp];
            if ((jp<nx-1)and(jp>1)and(ip<ny-1)and(ip>1)and(cl>0))or(it<5): ptin=1
            else : ptin=0
            vx,vy=self.core.flowReader.getLocalTransientV(self.core,infos,thick,ilay,ip,jp,iper)
            vx,vy=self.particle['direction']*vx,self.particle['direction']*vy;#forward or backward
            vx1,vx2=vx/poro;vy1,vy2=vy/poro;
            Ax = (vx2-vx1)/dx;vxm = 2*vx1-vx2;
            Ay = (vy2-vy1)/dy;vym = 2*vy1-vy2;

            x0m = (xp[it-1]-xg[jp]+dx);
            y0m = (yp[it-1]-yg[ip]+dy);
            vxp0 = vxm+Ax*x0m; vyp0 = vym+Ay*y0m;
            sensx = sign(vxp0);sensy = sign(vyp0)
            #* on differencie les cas */
            if ( (abs(vy1)+abs(vy2)) < ((abs(vx1)+abs(vx2))*cst1) ): # sens x
                dt = ((1.5+0.5*sensx)*dx-x0m)/vxp0*(1.+eps);dt1 = linspace(dt/4,dt,4)
                dxc = vxp0*dt; dyc = 0; jp += sensx
                dxc1 = vxp0*dt1; dyc1 = zeros(4)
            elif ( (abs(vx1)+abs(vx2)) < ((abs(vy1)+abs(vy2))*cst1) ): # sens y
                dt = ((1.5+0.5*sensy)*dy-y0m)/vyp0*(1.+eps);dt1 = linspace(dt/4,dt,4)
                dyc = vyp0*dt; dxc = 0; ip += sensy
                dyc1 = vyp0*dt1; dxc1 = zeros(4);
            else :
                lb1 = (vxp0-vxm)/x0m;
                lb2 = (vyp0-vym)/y0m;
                ax1 = max((lb1*dx+vxm)/vxp0,eps);
                ax2 = max((lb1*dx*2+vxm)/vxp0,eps);
                ay1 = max((lb2*dy+vym)/vyp0,eps);
                ay2 = max((lb2*dy*2+vym)/vyp0,eps);
                dtx1 = log(ax1)/lb1;dtx2 = log(ax2)/lb1;dtx = max(dtx1,dtx2);
                dty1 = log(ay1)/lb2;dty2 = log(ay2)/lb2;dty = max(dty1,dty2);
                if (dtx<=0) : dtx=1.e5; 
                if (dty<=0) : dty=1.e5;
                dt = min(dtx,dty)*(1+eps); #print it,dt,tp[it-1]
                dt1 = linspace(dt/4,dt,4)
                dxc1 = ( vxp0*exp(lb1*dt1)-vxm )/lb1-x0m; 
                dyc1 = ( vyp0*exp(lb2*dt1)-vym )/lb2-y0m;
                #if (dxc==0) & (dyc==0): break
            #* mis a jour des matrices */
            #print xp[it-1],dxc,xp[it-1]+dxc,yp[it-1],dyc,yp[it-1]+dyc
            cu[it] = cu[it-1] + sqrt(dxc1[-1]*dxc1[-1]+dyc1[-1]*dyc1[-1]);
            xp[it] = xp[it-1] + dxc1[-1];
            yp[it] = yp[it-1] + dyc1[-1];
            tp[it] = tp[it-1] + dt*self.particle['direction'];
            cu1.extend(cu[it-1]+sqrt(dxc1*dxc1+dyc1*dyc1))
            xp1.extend(xp[it-1]+dxc1)
            yp1.extend(yp[it-1]+dyc1)
            tp1.extend(tp[it-1]+dt1)
        #print 'addin 632',it,c_[xp[:it],yp[:it],tp[:it]]
        return array(xp1),array(yp1),array(tp1),array(cu1),it*4
        
    def calcPartInterp(self,xp0,yp0,iper):
        """make particle tracks with interpolated velocities in a permanent
        flow field. Works only for a regular grid"""
        poro = self.core.getValueFromName('Mt3dms','PRSTY')
        qx,qy,qz = self.core.flowReader.readFloFile(self.core,iper)
        qx,qy = qx[0]/poro,qy[0]/poro
        xp,yp,tp = [xp0],[yp0],[0]
        grd  = self.getFullGrid()
        x0,y0,nx,ny,dx,dy = grd['x0'],grd['y0'],grd['nx'],grd['ny'],grd['dx'][0], array(grd['dy'][0]);
        xmx,ymx = x0+nx*dx,y0+ny*dy
        it = 0
        
        while (it<1000)&(xp[it]>x0)&(xp[it]<xmx-dx)&(yp[it]>y0)&(yp[it]<ymx-dy): # 
            ip,jp = floor((yp[it]-y0)/dy),floor((xp[it]-x0)/dx)
            delty,deltx = (yp[it]-ip*dy)/dy,(xp[it]-jp*dx)/dx
            # velocity bilinear interpolation 
            u0 = qx[ip,jp]*(1-deltx)+qx[ip,jp+1]*deltx
            u1 = qx[ip+1,jp]*(1-deltx)+qx[ip+1,jp+1]*deltx
            u2 = u0*(1-delty)+u1*delty
            v0 = qy[ip,jp]*(1-delty)+qy[ip+1,jp]*delty
            v1 = qy[ip,jp+1]*(1-delty)+qy[ip+1,jp+1]*delty
            v2 = v0*(1-deltx)+v1*deltx
            dt = min(dx,dy)/sqrt(u2**2+v2**2)/4
            # move last point, to be close to the last line
            xp.append(xp[it]+u2*dt)
            yp.append(yp[it]+v2*dt)
            tp.append(tp[it]+dt)
            it += 1
        xp,yp,tp = array(xp),array(yp),array(tp)
        cu = r_[0,cumsum(sqrt((xp[1:]-xp[:-1])**2+(yp[1:]-yp[:-1])**2))]
        return xp,yp,tp,cu
        

    def calcPartMesh(self,vxin,vyin,xp0,yp0,ilay):
        """first attempt to calculate particles tracks in a mesh
        be careful this version computes only 80 triangles, the time is wrong
        and does not stop at boundaries..."""
        #iper=self.gui.guiShow.Tstep;#print 'addin 628 part',iper
        mesh = self.mesh
        lx,ly = [xp0],[yp0]
        xmn,xmx = amin(mesh.elx,axis=1),amax(mesh.elx,axis=1)
        ymn,ymx = amin(mesh.ely,axis=1),amax(mesh.ely,axis=1)
        lelt = where((xp0>xmn)&(xp0<xmx)&(yp0>ymn)&(yp0<ymx))[0] # thre can be several
        for iel in lelt:
            inod = mesh.elements[iel,-3:]
            inod1 = r_[inod, inod[0]] # to close the polygon
            poly = list(zip(mesh.nodes[inod1,1],mesh.nodes[inod1,2]))
            cf = array(mesh.lcoefs[iel]).T
            if pointsInPoly(array([xp0]),array([yp0]),poly,cf)[0]: break
            
        x,y = xp0,yp0
        for it in range(80): # loop on 80 triangles
            print(iel)
            vx,vy = vxin[iel],vyin[iel]
            inod=mesh.elements[iel,-3:]
            dref = sqrt(mesh.carea[iel])
            xno,yno=mesh.nodes[inod,1],mesh.nodes[inod,2]
            cf = array(mesh.lcoefs[iel]).T
            a,b = cf[0,:],cf[1,:]
            #equation y=(yo.vx-x0.vy+vy/a)/(vx+bvy/a)
            ym = (y*vx-x*vy+vy/a)/(vx+b*vy/a)
            xm = (1-b*ym)/a;#plot(x,y,'+')
            idx = where((sign(xm-x)==sign(vx))&(sign(ym-y)==sign(vy))&(abs(x-xm)>dref/1e3))[0]# the correct points
            if len(idx)==0:break#print x,y,xm,ym,idx
            if len(idx)>1: 
                dst = sqrt((x-xm)**2+(y-ym)**2)[idx]
                idx=idx[argsort(dst)][0]
            else :
                idx = idx[0]
            x,y = xm[idx],ym[idx]
            lx.append(x);ly.append(y)
            # find next triangle
            pts = inod[mod([idx,idx+1],3)]
            for n in mesh.cneighb[iel][:3]: # only first 3 others are for 3d
                if pts[0] in mesh.elements[n,-3:] and pts[1] in mesh.elements[n,-3:]:
                    iel = n
        return lx,ly,ones(len(lx))
        
class instant(object):
    def __init__(self, gui,core):
        self.gui, self.core = gui, core
        self.fitter = instantFit(self.core)
        
    def setObserver(self,gui,observer):
        self.gui, self.observer = gui, observer
        self.observer.bind_to(self.update)
        
    def startDialog(self):
        self.react = True # to specify if the fitter reacts or not
        self.first = True # to specify the first use of the isntant
        self.fitter.H_bc = (self.core.getValueLong('Modflow','bas.5',0)*(self.core.getValueLong('Modflow','bas.3',0)<0))[0]
        self.core.flowReader.readHeadFile = self.readHeadFile
        self.core.flowReader.getPtObs = self.getPtObsH
        self.core.transReader.readHeadFile = self.readUCN
        self.core.transReader.getPtObs = self.getPtObsC
        self.dic_options = {'type':'Head','zvalue':0,'aL':10.,'aT':1.}
        cfg = Config(self.core)
        dlg = cfg.dialogs.instantFitDialog(self.gui,self,self.dic_options)
        dlg.raise_()
        retour = dlg.show()
        
    def update(self, obs_value):
        #print 'instant gui',self,self.gui
        if self.react :
            self.fitter.calculate(self.dic_options,obs_value);print('addin 723, calculate')
            if self.dic_options['type']=='Tracer':
                xp1, yp1 = self.fitter.xp1,self.fitter.yp1
                # particle visu
#                nt,np = shape(xp1)
#                if self.first:
#                    for i in range(np):
#                        self.gui.visu.updateParticles(xp1[:,i],yp1[:,i],None) # adds a new particle in the list
#                    self.first = False
#                else :
#                    ldata = []
#                    for i in range(np):
#                        ldata.append((xp1[:,i],yp1[:,i],None))
#                    self.gui.visu.Particles['data'] = ldata
#                    self.gui.visu.changeParticles()
#                self.gui.visu.drawObject('Particles',True)
                data = [xp1,yp1,self.fitter.C]
                self.gui.visu.createAndShowObject(data,None,'contour')
            else :
                self.gui.guiShow.onClick2('Flow','Head',True)
            self.updateXyplot()
        
    def updateXyplot(self):
        '''updates the XY plots that were previously done with current values'''
        #print self.gui.guiShow.dicplots
        plt_h = self.gui.guiShow.dicplots['X_head'];#print plt_h
        if plt_h != None: # there is a graf XY for heads
            dist,val,labl = self.core.onPtObs('X0',0,'Flow','',['Head'],0);#zonenm '',layers 0
            plt_h.lignes[0]._points[:,1] = val[:,0]
            plt_h.cnv.Redraw()
        plt_t = self.gui.guiShow.dicplots['X_tracer'];#print plt_h
        if plt_t != None: # there is a graf XY for transport
            dist,val,labl = self.core.onPtObs('X0',0,'Transport','',['Transport'],0);#zonenm '',layers 0
            plt_t.lignes[0]._points[:,1] = val[:,0]
            plt_t.cnv.Redraw()
    
    def readHeadFile(self,core,tstep):
        return array(self.fitter.H,ndmin=3) # reader requires 3D data and self.H is 2D
    def readUCN(self,core,tstep):
        return array(self.fitter.C,ndmin=3) # reader requires 3D data and self.H is 2D

    def getPtObsH(self,core,iy,ix,iz,iper,esp):
        '''shall return a list of values for a list of points'''
        return self.fitter.H[iy,ix]
    def getPtObsC(self,core,iy,ix,iz,iper,opt):
        return self.fitter.getCfromReglist(ix,iy)
        
    def end(self): # does not work!!!!
        self.react = False
        self.core.flowReader = modflowReader(self.core.fileDir,self.core.fileName)
        
class oldstuff:
    def __init__(self):
        a=0
    def calcPartMesh_old(self,velo,xp0,yp0,ilay):
        """first attempt to calculate particles tracks in a mesh
        be careful this version computes only 80 triangles, the time is wrong
        and does not stop at boundaries..."""
        #iper=self.gui.guiShow.Tstep;#print 'addin 628 part',iper
        m1 = self.mesh
        vx, vy = velo #tci.gradient(mesh.trg_nodes.x, mesh.trg_nodes.y)
        lx,ly = [xp0],[yp0]
        xmn,xmx = amin(m1.elx,axis=1),amax(m1.elx,axis=1)
        ymn,ymx = amin(m1.ely,axis=1),amax(m1.ely,axis=1)
        lelt = where((xp0>xmn)&(xp0<xmx)&(yp0>ymn)&(yp0<ymx))[0] # thre can be several
        for iel in lelt:
            inod = m1.elements[iel,-3:]
            inod1 = r_[inod, inod[0]] # to close the polygon
            poly = list(zip(m1.nodes[inod1,1],m1.nodes[inod1,2]))
            cf = array(m1.lcoefs[iel]).T
            if pointsInPoly(array([xp0]),array([yp0]),poly,cf)[0]: break
                
        x,y = xp0,yp0
        for it in range(80): # loop on 80 triangles
            print(iel)
            inod=m1.elements[iel,-3:]
            xn,yn=m1.nodes[inod,1],m1.nodes[inod,2]
            #un,vn = velo[0][0,nodes],velo[1][0,nodes]
            un,vn = vx[inod],vy[inod]
            # coords of u ,v planes in the x,y dimension
            M = ones((3,3));M[:,0]=xn;M[:,1]=yn
            ucoef,vcoef = solve(M,un),solve(M,vn) 
            # velocity and position at the start point
            u0 = x*ucoef[0]+y*ucoef[1]+ucoef[2]
            v0 = x*vcoef[0]+y*vcoef[1]+vcoef[2]
            #find dt
            dtx,dty = (max(xn)-min(xn))/abs(u0),(max(yn)-min(yn))/abs(v0)
            dt=min(dtx,dty)/6
            cf = array(m1.lcoefs[iel]).T
            pos0 = sign(dot(array([x,y]),cf)-1)
            for i in range(50):
                pos = sign(dot(array([x,y]),cf)-1);#print i,pos
                if any(pos!=pos0): break
                u1 = x*ucoef[0]+y*ucoef[1]+ucoef[2]; 
                v1 = x*vcoef[0]+y*vcoef[1]+vcoef[2];
                x += u1*dt*self.particle['direction']
                y += v1*dt*self.particle['direction']
                lx.append(x*1);ly.append(y*1); #print 'pos', pos, pos0
            if all(pos==pos0): return lx,ly,ones(len(lx)) #no movement
            idff = where(abs(pos-pos0))[0][0] # the index of the line that has been crossed
            # move last point, to be on the last line
            lcx,lcy = cf[0,idff],cf[1,idff]
            #equation lcx(x-a.u1.dt)+lcy(y-a.v1.dt)=1
            a = (x*lcx+y*lcy-1)/(u1*dt*lcx+v1*dt*lcy)*.99
            x,y = x-a*u1*dt, y-a*v1*dt;#plot(x,y,'+')
            # find next triangle
            pts = inod[mod([idff,idff+1],3)]
            neighbs = m1.cneighb[iel][:3] ; #print neighbs# still in 2D, take only the first three neighbors (same layer)
            for n in neighbs:
                if pts[0] in m1.elements[n,-3:] and pts[1] in m1.elements[n,-3:]:
                    iel = n
        return lx,ly,ones(len(lx))
        
    def calcPartMesh0(self,K,H,xp0,yp0,ilay):
        """first attempt to calculate particles tracks in a mesh
        be careful this version computes only 80 triangles, the time is wrong
        and does not stop at boundaries..."""
        #iper=self.gui.guiShow.Tstep;#print 'addin 628 part',iper
        mesh = self.mesh
        lx,ly = [xp0],[yp0]
        xmn,xmx = amin(mesh.elx,axis=1),amax(mesh.elx,axis=1)
        ymn,ymx = amin(mesh.ely,axis=1),amax(mesh.ely,axis=1)
        lelt = where((xp0>xmn)&(xp0<xmx)&(yp0>ymn)&(yp0<ymx))[0] # thre can be several
        for iel in lelt:
            inod = mesh.elements[iel,-3:]
            inod1 = r_[inod, inod[0]] # to close the polygon
            poly = list(zip(mesh.nodes[inod1,1],mesh.nodes[inod1,2]))
            cf = array(mesh.lcoefs[iel]).T
            if pointsInPoly(array([xp0]),array([yp0]),poly,cf)[0]: break
            
        x,y = xp0,yp0
        for it in range(80): # loop on 80 triangles
            print(iel)
            inod=mesh.elements[iel,-3:]
            #xno,yno=mesh.nodes[inod,1],mesh.nodes[inod,2]
            # find gradH
            lelt = [];
            for io in inod: lelt.extend(mesh.nd_elt[io])
            lelt=unique(lelt)
            xyc = mesh.elcenters[lelt,:]
            Hc = H[lelt] # head at cell centers
            cf = array(mesh.lcoefs[iel]).T
            pos0 = sign(dot(array([x,y]),cf)-1)
            pos,dref = pos0,sqrt(mesh.carea[iel])
            for i in range(50):
                pos = sign(dot(array([x,y]),cf)-1);#print i,x,y,pos
                if any(pos!=pos0): break
                d0 = sqrt((x-xyc[:,0])**2+(y-xyc[:,1])**2) # dist from present pt to cell centers
                h0 = sum(Hc/d0)/sum(1/d0)
                dx,dy = dref/10,dref/10
                ddx,ddy = sqrt((x+dx-xyc[:,0])**2+(y-xyc[:,1])**2),sqrt((x-xyc[:,0])**2+(y+dy-xyc[:,1])**2)
                h0dx,h0dy = sum(Hc/ddx)/sum(1/ddx),sum(Hc/ddy)/sum(1/ddy)
                vx,vy = -K[iel]*(h0dx-h0)/dx, -K[iel]*(h0dy-h0)/dy
                dt = min(abs(dref/vx),abs(dref/vy))/10
                x += vx*dt*self.particle['direction']
                y += vy*dt*self.particle['direction'];#plot(x,y,'+')
                lx.append(x*1);ly.append(y*1); #print 'pos', pos, pos0
            if all(pos==pos0): return lx,ly,ones(len(lx)) #no movement
            idff = where(abs(pos-pos0))[0][0] # the index of the line that has been crossed
            # move last point, to be on the last line if we crossed two lines
            lcx,lcy = cf[0,idff],cf[1,idff]
            #equation lcx(x-a.u1.dt)+lcy(y-a.v1.dt)=1
            a = (x*lcx+y*lcy-1)/(vx*dt*lcx+vy*dt*lcy)*.99
            x,y = x-a*vx*dt, y-a*vy*dt;#plot(x,y,'+')
            lx[-1],ly[-1] = x,y
            # find next triangle
            pts = inod[mod([idff,idff+1],3)]
            for n in mesh.cneighb[iel][:3]: # only first 3 others are for 3d
                if pts[0] in mesh.elements[n,-3:] and pts[1] in mesh.elements[n,-3:]:
                    iel = n
        return lx,ly,ones(len(lx))
