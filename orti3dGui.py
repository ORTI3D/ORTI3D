"""this is the main gui for ORTi3D. It creates a wx window where parameter
panel, visualisation panel and visuChoose panel are included.
it uses core for the major job of storing and retrieving the data
writer/readers are available for different models"""
import os, sys,traceback # traceback added OA 25/9/18
from PyQt5.QtCore import *
from PyQt5.QtGui import *
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
#from PyQt5.QtWidgets import *
from qtVisualisation import *
from qtShow import *
from qtParameters import *
from qtTopBar import *
from qtDialogs import *
from menus import *
from core import *
from addin import *
#
from defSPH import *
from configSPH import *
#from pyqtconsole.console import PythonConsole # does not work, pyqt console is nto present
#from subprocess import CREATE_NEW_CONSOLE

class sendMessage: # OA added 13/12/19
    def __init__(self,gui): self.gui=gui
    def write(self,txt):
        onMessage(self.gui,txt)

class orti3dGui(QMainWindow):
    
    def __init__(self,title):
        super(orti3dGui, self).__init__()
        #self.Maximize(True)
        self.gtyp = 'qt'
        self.title = title
        self.screenShape = QDesktopWidget().screenGeometry()
        #self.icons = self.makeIcons()
        lfi=os.listdir(os.getcwd()) # oa 25/5
        if 'menus.py' in lfi: 
            self.mainDir,self.u_dir = os.path.dirname(os.getcwd()),'..'+os.sep+'utils' # python case
        else : 
            self.mainDir,self.u_dir = os.getcwd(),'utils'
        self.core = Core(self)
        self.linesDic,self.linesCommDic = {},{}
        for mod in self.core.modelList:
            self.linesDic[mod] = self.core.diczone[mod].getLinesDic()
            self.linesCommDic[mod] = self.core.diczone[mod].getLinesCommDic()
        #print self.linesDic
        self.makePanelMatplotlib()
        self.makeTopBar()
        self.makePanelParameters()
        self.makePanelShow()
        self.makeMenus()
        self.core.addin.setGui(self)
        self.core.addin.initMenus()
        #self.BackgroundColour = (230,230,230)
                
        frameSizer = QHBoxLayout()
        frameSizer.addItem(self.paramSizer)
        frameSizer.addItem(self.matplot)
        frameSizer.addItem(self.showSizer)

        globalSizer = QVBoxLayout()
        globalSizer.addItem(self.topSizer)
        globalSizer.addItem(frameSizer)
        globalSizer.addItem(self.basSizer)
        
        self.widget = QWidget()
        self.widget.setLayout(globalSizer)

        self.setCentralWidget(self.widget)
        self.setWindowTitle(title)
        sys.excepthook = self.myExceptionHandler
        #sys.stdout=sendMessage(self) # OA added 13/12/19

    def myExceptionHandler(self, type, value, trace_b): # OA added 25/9/18
        """Catch exceptions and show error dialog"""
        f = traceback.format_tb(trace_b, limit=10) # OA modif 1/10
        onMessage(self,'\n'.join(f[-10:])+'\n'+str(value))
        
    def onInitGui(self,core):
        self.on3D(core.addin.getDim()=='3D')
        self.onSetMediaNb(getNmedia(core),getNlayers(core))
        self.onRCT(False)
        self.guiShow.init()
        self.guiShow.openModel()
        self.core.addin.resetAddin()

    def on3D(self,bool):
        boutonVisible(self,'Ad_3D',bool) #EV 05/08/19
        self.varBox.choice3D.setEnabled(bool) #EV 05/08/19
        
    def onParticle(self,bool): # EV 3/12/21
        boutonVisible(self,'Ad_Particle',bool) 
        
    def onSetMediaNb(self,nbM,nbL):
        self.varBox.choice3D.clear()
        for i in range(nbM): self.varBox.choice3D.addItem(str(i))
        if self.core.addin.getDim()=='3D':
            self.guiShow.setNames('Model_Layer_L',list(range(nbL)))
        
    def onRCT(self,bool):
        boutonVisible(self,'Ad_MtSpecies',bool) #EV 05/08/19
        boutonVisible(self,'Ad_MtReact',bool) #EV 05/08/19
        
    def onGridMesh(self,gridtyp):  # OA 22/8/19
        boutonIcon(self,'Ad_Grid','Ad_'+gridtyp+'.png')
        
    def onWriteModflow(self,check): #EV 24/04/20
        if check == 2 : check = False
        else : check = True
        boutonVisible(self,'Fl_Write',check) 
        
    def onMessage(self,text): onMessage(self,text)  # OA 18/12/21

    ####################################################
    #                   make menus
    ####################################################
    def makeMenus(self):
        self.menus = Menus(self, self.core)
        #file menu
        menuFile = self.menuBar().addMenu("&File")
        newAction = QAction("&New", self, shortcut="CTRL+n",
                statusTip="create new file",triggered=self.menus.OnNew)
        openAction = QAction("&Open", self, shortcut="CTRL+o",
                statusTip="open file",triggered=self.menus.OnOpen)
        saveAction = QAction("&Save", self, shortcut="CTRL+s",
                statusTip="save this file",triggered=self.menus.OnSave)
        saveAsAction = QAction("&Save as", self, shortcut="save as",
                statusTip="save as another file",triggered=self.menus.OnSaveAs)
        exitAction = QAction("&Exit", self, shortcut="Ctrl+X",
                statusTip="Quit", triggered=self.closeEvent)
        menuFile.addAction(newAction)
        menuFile.addAction(openAction)
        menuFile.addAction(saveAction)
        menuFile.addAction(saveAsAction)
        menuFile.addSeparator()
        menuFile.addAction(exitAction)
        # import menu
        menuImport = self.menuBar().addMenu("&Import")
        i3dgeoAction = QAction("&3Dgeom", self, 
                statusTip="import disrw file",
                triggered=self.menus.OnImport3DgeomDis)
        iuspecAction = QAction("&User species", self, #EV 11/12/19
                statusTip="import user species text file",
                triggered=self.menus.OnImportUserSpecies)
        ipostfAction = QAction("&Postfix species", self, # OA 30/6/19
                statusTip="import user species from postfix",
                triggered=self.menus.OnImportPostfixSpecies)
        menuImport.addAction(i3dgeoAction)
        menuImport.addAction(iuspecAction)
        menuImport.addAction(ipostfAction)        # OA 30/6/19
        # export menu
        menuExport = self.menuBar().addMenu("&Export")
        eparmAction = QAction("&Current Parameter", self, 
                statusTip="export current parameter as text file",
                triggered=self.menus.OnExportParm)
        eresuAction = QAction("&Current Results", self, 
                statusTip="export current results as text file",
                triggered=self.menus.OnExportResu)
        eparmVAction = QAction("&Current Parameter Vtk", self, 
                statusTip="export current parameter as text file",
                triggered=self.menus.OnExportParmVtk)
        eresuVAction = QAction("&Current Result Vtk", self, 
                statusTip="export current results as text file",
                triggered=self.menus.OnExportResuVtk)
        eresuVaAction = QAction("&Current Result Vtk, all times", self, 
                statusTip="export current results as text file",
                triggered=self.menus.OnExportResuVtkAll)
        evectVAction = QAction("&Current Vector Vtk", self, 
                statusTip="export current vectors as text file",
                triggered=self.menus.OnExportVectorVtk)
        menuExport.addAction(eparmAction)
        menuExport.addAction(eresuAction)
        menuExport.addAction(eparmVAction)
        menuExport.addAction(eresuVAction)
        menuExport.addAction(eresuVaAction)
        menuExport.addAction(evectVAction)
       
       #Tools
        menuTools = self.menuBar().addMenu("&Tools")
        tnodesAction = QAction("&Node search", self, 
                statusTip="Node number",
                triggered=self.menus.OnNodeSearch)
        menuTools.addAction(tnodesAction)
        
       #Add-ins
        self.menuAddin = self.menuBar().addMenu("&Addin")
        
       #Plugins
        menuPlugins = self.menuBar().addMenu("&Plugins")
        self.core.plugins.setGui(self)
        pl_list= self.core.plugins.pl_list
        lactions = []
        for pl_n in pl_list:
            act = QAction("&"+pl_n,self,statusTip=n,
                triggered=self.core.plugins.action)
            menuPlugins.addAction(act)

        #Help
        menuHelp = self.menuBar().addMenu("&Help")
        hActionI = QAction("&Help interface", self,
                statusTip="Help file",
                triggered=self.menus.OnHelpI)
        hActionM = QAction("&Help models", self,
                statusTip="Help file",
                triggered=self.menus.OnHelpM)
        dwnsAction = QAction("&Download stable", self, 
                statusTip="download stable version from github",
                triggered=self.menus.OnDownloadLast)
        #dwndAction = QAction("&Download develop", self, 
                #statusTip="download development version from github",
                #triggered=self.menus.OnDownloadDev)
        #dwnlAction = QAction("&Download local", self, #EV 11/12/19
                #statusTip="download any version from local file",
                #triggered=self.menus.OnDownloadLocal)
        hActionPy = QAction("&Python console", self,
                statusTip="py console",
                triggered=self.menus.OnHelpPy)
        menuHelp.addAction(hActionI)
        menuHelp.addAction(hActionM)
        menuHelp.addAction(dwnsAction)
        menuHelp.addAction(hActionPy)
        #menuHelp.addAction(dwndAction) #EV 11/12/19
        #menuHelp.addAction(dwnlAction)

    def enableMenu(self,nomM,bool):
        pass
        '''id=self.menuBar.FindMenu(nomM)
        if id!=-1:self.menuBar.EnableTop(id,bool)  # pour les griser'''
        
    def addMenu(self,num,menuName,methd):
        action = QAction(menuName, self, triggered=methd)
        self.menuAddin.addAction(action)
        
    def updateTitle(self):
        self.setWindowTitle(self.title + " - " + self.core.fileDir+'/'+ self.core.fileName)
                               
    def closeEvent(self,evt=None):
        self.menus.askSave(evt)
        sys.exit()
    #####################################################
    #                   Panel Matplotlib
    ######################################################
    def makePanelMatplotlib(self):
        self.matplot = QVBoxLayout()
        self.visu = qtVisualisation(self)
        #self.visu.setGeometry(QRect(0,0,760,500))
        self.matplot.addWidget(self.visu)
        
        self.basSizer = QHBoxLayout()
        self.basSizer.addWidget(self.visu.GetToolBar())
        self.pos = QLabel("Media")
        #self.pos.SetOwnBackgroundColour('WHITE')
        self.basSizer.addWidget(self.pos)        
        self.notify = QLabel("")
        #font = wx.Font(16, wx.SWISS, wx.NORMAL, wx.NORMAL)
        #self.notify.SetFont(font)
        self.basSizer.addWidget(self.notify)
        self.onNotify("run time :" )
        self.visu.initDomain() #,self.model.getGlist())
        
    #def getVisu(self): return self.visu
    def onNotify(self,text): self.notify.setText(text)
    def onPosition(self,text): self.pos.setText(text)
    
    #####################################################
    #                   Panel Top et Parametres
    #####################################################
    def makeTopBar(self):
        self.currentModel,self.currentLine,self.currentMedia = 'Modflow',None,0
        self.topSizer = QHBoxLayout()
        self.topSizer.setGeometry(QRect(0,0,int(self.screenShape.width()*0.8),45)) #900
        width = self.screenShape.width()*0.8
        self.topSizer.setSpacing(2)
        policy = QSizePolicy(QSizePolicy.Preferred, QSizePolicy.Preferred)
        qwt1 = QGroupBox()
        qwt1.setSizePolicy(policy)
        qwt1.setMinimumHeight(45)
        qwt1.setMinimumWidth(int(width*0.55))
        self.varBox = Ui_Var()
        self.varBox.setupUi(qwt1,self,self.core)
        self.topSizer.addWidget(qwt1) #,0,0,1,4) 
        qwt2 = QGroupBox()
        qwt2.setSizePolicy(policy)
        qwt2.setMinimumHeight(45)
        qwt2.setMaximumWidth(140)
        self.addBox = Ui_AddZone()
        self.addBox.setupUi(qwt2,self,self.core)
        self.topSizer.addWidget(qwt2) #,0,3,1,1)
        qwt3 = QGroupBox()
        qwt3.setSizePolicy(policy)
        qwt3.setMinimumHeight(45)
        self.modifBox = Ui_ModifZone()
        self.modifBox.setupUi(qwt3,self,self.core)
        self.topSizer.addWidget(qwt3)#,0,5,1,3)

    def makePanelParameters(self):
        self.paramSizer = QVBoxLayout()
        qwp = QWidget()
        qwp.setMinimumWidth(int(self.screenShape.width()*0.085))
        self.dlgParameters = Ui_Parameters()
        self.dlgParameters.setupUi(qwp,self,self.core,self.mainDir)
        self.paramSizer.addWidget(qwp)
        it1=qwp.findChild(QPushButton,'Ad_MtSpecies') #EV 05/08/19
        it2=qwp.findChild(QPushButton,'Ad_MtReact')
        it3=qwp.findChild(QPushButton,'Ad_3D')
        it1.setEnabled(False) #EV 05/08/19
        it2.setEnabled(False)
        it3.setEnabled(False)
        #self.onRCT(False)
        
    #####################################################
    #                   Panel Vue
    #####################################################
    def makePanelShow(self):
        #self.guiShow = guiShow(self,self.core)
        self.showSizer = QVBoxLayout()
        qws = QWidget()
        qws.setMinimumWidth(int(self.screenShape.width()*0.112))
        self.dlgShow = Ui_Show()
        self.dlgShow.setupUi(qws,self,self.core)     
        self.guiShow = self.dlgShow.guiShow
        self.showSizer.addWidget(qws)
        self.guiShow.dlgShow.onTickBox('Model','Grid','B',True) #EV 15/01/21
        self.guiShow.dlgShow.setLine('Chemistry_Units_L',1) #EV 16/12/21
        
    ######################## actions ############################
    def actions(self,action):
        if action == 'zoneStart': self.panelsEnable(False)
        if action == 'zoneEnd': self.panelsEnable(True)

    def panelsEnable(self,bool):
        pass
        '''self.qtParameters.Enable(bool);self.guiShow.dlgShow.Enable(bool)
        self.varBox.Enable(bool);#self.addBox.Enable(bool);self.modifBox.Enable(bool)
        '''
        
    ######################## ICONS ############################
    import sys
    os.path.join(os.path.dirname(sys.executable), 'utils')

    ######################## Common variables for the interface ###################""""
    def common(self):
        self.currentModel,self.currentGroup,self.currentLine,self.currentMedia = None, None, None,None

    def resetCurVar(self):
        global curVar
        curVar = {}        

class Observer(object):
    '''this tool serve to observe what happens somewhere and when some object
    is in the observers list it is instantaneously modified'''
    def __init__(self):
        self._obs_value = 1
        self._observers = []
    def get_obs(self): return self._obs_value
    def bind_to(self, callback): self._observers.append(callback)
    def set_obs(self, value):
        self._obs_value = value
        for callback in self._observers:
            callback(self._obs_value)
    obs_value = property(get_obs, set_obs)
