# -*- coding: utf-8 -*-

from PyQt5.QtGui import *
import os
from .qtDialogs import *
from .geometry import *
from .core import *
from .config import *
from .topBar import BaseTop

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

class Ui_Var(object):
    def setupUi(self,Var,gui,core):
        self.gui,self.core = gui,core
        self.base = BaseTop(gui,core)
        self.cfg = Config(core)
        self.screenShape = QDesktopWidget().screenGeometry() # get the pc screen shape
        width = self.screenShape.width()*0.8
        Var.setObjectName("Var")
        Var.setTitle('Spatial attributes')
        Var.setGeometry(QRect(0, 0,(width*0.55), 45))
        self.hlWidget = QWidget(Var) #(self.group)
        self.hlWidget.setGeometry(QRect(5, 20,(width*0.55), 25))
        self.gridLayout = QHBoxLayout(self.hlWidget)
        self.gridLayout.setContentsMargins(1,1,10,1)
        self.gridLayout.setSpacing(2)

        #policy = QSizePolicy(QSizePolicy.Preferred, QSizePolicy.Preferred)
        self.choiceM = QComboBox(self.hlWidget)
        #self.choiceM.setSizePolicy(policy)
        self.choiceM.setMaximumWidth(100);#setGeometry(QRect(0, 0, 45, 18))
        view =  QListView() ; view.setMinimumWidth(150) # creat and set the width ListView
        self.choiceM.setView(view)
        self.choiceM.activated['QString'].connect(self.onChoiceModel)
        self.gridLayout.addWidget(self.choiceM, 0)

        self.choiceG = QComboBox(self.hlWidget)
        #self.choiceG.setSizePolicy(policy)
        self.choiceG.setMaximumWidth(65);
        self.choiceG.activated['QString'].connect(self.onChoiceGroup)
        self.gridLayout.addWidget(self.choiceG, 1)

        self.choiceL = QComboBox(self.hlWidget)
        self.choiceL.setMaximumWidth(100);
        view2 =  QListView(self.choiceL) ; view2.setMinimumWidth(250) # creat and set the width ListView
        self.choiceL.setView(view2)
        #self.choiceL.setSizePolicy(policy)
        self.choiceL.activated['QString'].connect(self.onChoiceLine)
        self.gridLayout.addWidget(self.choiceL,2)

        label = QLabel(self.hlWidget)
        label.setText("     Media")#;label.setSizePolicy(policy);label.setMaximumWidth(30)
        #label.setFixedWidth(55)
        self.gridLayout.addWidget(label)
        self.choice3D = QComboBox(self.hlWidget)
        #self.choice3D.setSizePolicy(policy)
        self.choice3D.setMaximumWidth(55);
        self.choice3D.activated['QString'].connect(self.onChoiceMedia)
        self.choice3D.activated['QString'].connect(self.onViewVariable)
        self.gridLayout.addWidget(self.choice3D, 3)

        label = QLabel(self.hlWidget)
        label.setText("     Backg.")#;label.setSizePolicy(policy);label.setMaximumWidth(30)
        #label.setFixedWidth(55)
        self.gridLayout.addWidget(label)
        self.backg = QLineEdit(self.hlWidget)
        #self.backg.setSizePolicy(policy)
        self.backg.setMaximumWidth(55);
        #self.backg.returnPressed.connect(self.onBackOk) # OA removed 25/9/18
        self.gridLayout.addWidget(self.backg)
        
        self.butOK = QPushButton(self.hlWidget) # OA added 25/9/18
        self.butOK.setMaximumWidth(30)    # OA added 25/9/18         
        self.butOK.setText('Ok')    # OA added 25/9/18         
        self.butOK.clicked.connect(self.onBackOk)  # OA added 25/9/18
        self.gridLayout.addWidget(self.butOK)  # OA added 25/9/18

        label = QLabel(self.hlWidget)
        label.setText("     Type")#;label.setSizePolicy(policy);label.setMaximumWidth(30)
        #label.setFixedWidth(45)
        self.gridLayout.addWidget(label)
        self.choiceT = QComboBox(self.hlWidget)
        view3 =  QListView() ; view3.setMinimumWidth(100) # creat and set the width ListView
        self.choiceT.setView(view3)
        #self.choiceT.setSizePolicy(policy)
        #self.choiceT.setMaximumWidth(120);
        self.typeList = ['one_value','zone','formula','interpolate','importArray','importZones'] # EV 04/02/20
        self.choiceT.addItems(self.typeList)
        self.choiceT.activated['QString'].connect(self.onChoiceType)
        self.gridLayout.addWidget(self.choiceT)#, 0, 8, 1, 1)

        label = QLabel(self.hlWidget)
        label.setText("     View")#;label.setSizePolicy(policy);label.setMaximumWidth(25)
        #label.setMaximumWidth(45)
        self.gridLayout.addWidget(label)
        self.choiceV = QComboBox(self.hlWidget) #EV 26.11.20
        self.choiceV.setMaximumWidth(55);
        self.chkView = QCheckBox(self.hlWidget)  
        self.choiceV.activated['QString'].connect(self.onViewVariable)
        self.chkView.stateChanged.connect(self.onViewVariable)
        self.gridLayout.addWidget(self.chkView)
        self.gridLayout.addWidget(self.choiceV)
        
        self.retranslateUi(Var)
        QMetaObject.connectSlotsByName(Var)
        
        self.blind = {}
        for k in self.core.modelList: self.blind[k] = []
        self.blind['Mt3dms']=['btn.9','btn.10','uzt.3','uzt.4']
        self.currentMedia = 0

    def retranslateUi(self, Var): pass
    
    def setChoiceList(self,obj,l):
        obj.clear()
        obj.addItems(l)

    def chooseCategory(self,c0):
        self.currentCategory = c0;#print c0,self.core.modelList
        lmodels = self.base.modlistFromGroup(c0)
        self.setChoiceList(self.choiceM,lmodels)
        
    def onChoiceModel(self,evt):
        """contains the models : for modflow series : modflow, mt3d, pht3d..."""
        model = str(self.choiceM.currentText()) # OA 28/7/19 removed 4 lines for upper/lower
        self.gui.currentModel = model
        lmodules = self.core.addin.getUsedModulesList(model); #print 'topbar,l84',lmodules# to set the groups
        lmodules = self.selectGroups(model,lmodules)
        self.setChoiceList(self.choiceG,lmodules)
        
    def selectGroups(self,model,lmodules):
        '''to show only the groups that contain arrays'''
        lmodOut = []
        for grp in lmodules:
            if grp in list(self.gui.linesDic[model].keys()): 
                lmodOut.append(grp)
        return lmodOut
        
    def onChoiceGroup(self,evt):
        '''the group is the modules : DIS...'''
        curGroup = str(self.choiceG.currentText())
        self.gui.currentGroup = curGroup
        self.choiceL.clear()
        if curGroup not in list(self.gui.linesDic[self.gui.currentModel].keys()) : return
        lines = self.gui.linesDic[self.gui.currentModel][curGroup]
        indx = self.base.testConditions(self.gui.currentModel,lines)
        lcomm = self.gui.linesCommDic[self.gui.currentModel][curGroup]
        llines = []
        for il in indx: 
            if lines[il] in self.blind[self.gui.currentModel]: continue
            self.choiceL.addItem("")
            llines.append(lines[il]+' '+lcomm[il])
        self.setChoiceList(self.choiceL,llines)
            #self.choiceL.setItemText(il) #, QApplication.translate("Type", n, None, QApplication.UnicodeUTF8))
        self.gui.currentLine = lines[0]
        
    def onChoiceLine(self,evt):
        '''the line is dis.1 or lpf.8..'''
        line = str(self.choiceL.currentText()).split()[0]
        self.gui.currentLine = line
        media = self.gui.currentMedia
        vallist = self.core.dicval[self.gui.currentModel][line];#print 'topbar valist',vallist
        nval = len(vallist)
        nmedia = getNmedia(self.core)
        if nval<nmedia : vallist.extend([vallist[0]]*(nmedia-nval))
        self.backg.setText(str(vallist[media])) # set to the new value
        typ = self.core.dictype[self.gui.currentModel][line]#[0] # EV 3/2/20
        ntyp = len(typ) # EV 3/2/20
        if ntyp<nmedia : typ.extend([typ[0]]*(nmedia-ntyp)) # EV 3/2/20
        i = self.choiceT.findText(typ[media]) # EV 3/2/20
        self.choiceT.setCurrentIndex(i)
        self.gui.modifBox.updateChoiceZone(line)
        self.base.changeVisu()
        self.onSetVariable() #EV 26.11.20
            
    def onChoiceMedia(self,evt):
        """changes the media in visualization and stores the current media"""
        media = self.choice3D.currentIndex();#print media
        self.gui.currentMedia = media
        line = self.gui.currentLine
        vallist = self.core.dicval[self.gui.currentModel][line]
        self.backg.setText(str(vallist[media]))
        typ = self.core.dictype[self.gui.currentModel][line] # EV 3/2/20
        i = self.choiceT.findText(typ[media]) # EV 3/2/20
        self.choiceT.setCurrentIndex(i) # EV 3/2/20
        self.base.changeVisu()
        
    def onBackOk(self):
        line = self.gui.currentLine
        media = self.gui.currentMedia
        self.core.setValue(self.gui.currentModel,line,media,float(self.backg.text()))
        
    def onChoiceType(self):
        line = self.gui.currentLine
        media = self.gui.currentMedia # EV 3/2/20
        choice = str(self.choiceT.currentText()); #print choice
        self.core.dictype[self.gui.currentModel][line][media] = choice # EV 3/2/20
        if choice =='zone': self.base.changeVisu()
        elif choice == 'formula': self.base.onFormula()
        elif choice == 'interpolate': self.base.onInterpolate()
        elif choice == 'importArray' : self.base.onImportArray() # EV 3/2/20
        elif choice == 'importZones': 
            self.base.onImportZones()
            self.core.dictype[self.gui.currentModel][line][media] = 'zone' # EV 3/2/20
        if line in list(self.cfg.curVar.keys()): 
            self.cfg.curVar.pop(line);#print line,self.curVar.keys() # when a type is selected, it removes the stored  view

    def onEdit(self):
        pass
    
    def onSetVariable(self): #EV 26.11.20
        line = str(self.choiceL.currentText()).split()[0]
        if line in ['drn.1','ghb.1']: lvar =['1','2']
        elif line=='riv.1': lvar = ['1','2','3']
        else :lvar = ['1']
        self.setChoiceList(self.choiceV,lvar)
        #self.onViewVariable()
    
    def onViewVariable(self):
        """used to see the current variable for current medium"""
        var = self.choiceV.currentIndex() #EV 26.11.20
        if self.chkView.isChecked(): 
            X,Y,mat = self.base.getCurVariable(var)
            self.gui.visu.createImage([X,Y,mat])
            self.gui.visu.drawImage(True)
        else :
            self.gui.visu.drawImage(False)
            
            
class Ui_AddZone(object):
    def setupUi(self,AddZone,gui,core):
        AddZone.setObjectName("AddZone")
        self.gui,self.core = gui,core
        self.AddZone,self.visu = AddZone,self.gui.visu
        self.base = BaseTop(gui,core)
        self.currentMedia = gui.currentMedia
        dirutils = self.gui.mainDir+os.sep+'utils'
        
        AddZone.setTitle('Add Zones')
        AddZone.setGeometry(QRect(0, 0, 135, 45))
        self.hlWidget = QWidget(AddZone)
        self.hlWidget.setGeometry(QRect(4, 20, 135, 25))
        zoneSizer = QHBoxLayout(self.hlWidget)
        zoneSizer.setContentsMargins(1,1,1,1)
        zoneSizer.setSpacing(2)

        policy = QSizePolicy(QSizePolicy.Preferred, QSizePolicy.Preferred)
        shortNames = ['Top_Point','Top_Line','Top_Rect','Top_Poly','Top_PolyV']
        for n in shortNames:
            but = QPushButton(self.hlWidget)
            #but.setContentsMargins(0,0,0,0)
            but.setToolTip(n)
            icon = QIcon()
            icon.addPixmap(QPixmap(dirutils+os.sep+n+'.png'), QIcon.Normal, QIcon.Off)
            but.setIcon(icon)
            but.setIconSize(QSize(25, 25))
            #but.setSizePolicy(policy)
            but.setMaximumWidth(25)
            but.setObjectName(n.split('_')[1].upper()) #_fromUtf8(n))
            but.clicked.connect(self.onShape)
            zoneSizer.addWidget(but)
            but.setFlat(True)
        
    def onShape(self):
        line = self.gui.currentLine
        obj = self.AddZone.sender()
        shap = str(obj.objectName());
        if (line == None) or (line not in self.gui.linesDic[self.gui.currentModel][self.gui.currentGroup]):
            onMessage(self.gui,"Choose one variable")
            return
        self.gui.actions('zoneStart')
        exec('self.visu.setZoneReady(\"'+shap+'\",\"'+line+'\")')

    def onZoneCreate(self, typeZone, xy):
        self.base.onZoneCreate(typeZone, xy)
            
#////////////////////////////////////////////////////:::
class Ui_ModifZone(object):
    def setupUi(self,ModifZone,gui,core):
        ModifZone.setObjectName("ModifZone")
        self.gui,self.core = gui,core
        self.ModifZone,self.visu = ModifZone,self.gui.visu
        self.base = BaseTop(gui,core)
        self.currentModel = gui.currentModel
        dirutils = self.gui.mainDir+os.sep+'utils';#WonMessage(self.gui,dirutils)
        
        ModifZone.setTitle('Modify Zones')
        ModifZone.setGeometry(QRect(0, 0, 480, 45))
        self.hlWidget = QWidget(ModifZone)
        self.hlWidget.setGeometry(QRect(4, 20, 475, 25))
        zoneSizer = QHBoxLayout(self.hlWidget)
        zoneSizer.setContentsMargins(1,1,1,1)
        zoneSizer.setSpacing(2)

        #policy = QSizePolicy(QSizePolicy.Preferred, QSizePolicy.Preferred)
        self.zmodif,self.zmove,self.zindex=0,0,[]
        self.currentZlist = None; #EV 14/08/19
    ## choice de la zone
        self.choice = QComboBox(self.hlWidget)
        #self.choice.setSizePolicy(policy)
        self.choice.setMaximumWidth(85)  
        zoneSizer.addWidget(self.choice)
        self.choice.activated['QString'].connect(self.onChoice)
    ## valeur de la zone selectionnee
        self.valZ = QPushButton(self.hlWidget);zoneSizer.addWidget(self.valZ)
        #self.valZ.setSizePolicy(policy)
        self.valZ.setMaximumWidth(85)             
        self.valZ.clicked.connect(self.onValueZ)
    ##bouton  milieu pour zone, puis deplacer, puis modif
        for n in ['move','modifPoly','supprime','supprimeAll']:       
            but = QPushButton(self.hlWidget);
            but.setObjectName(n)
            but.setToolTip(n)
            #but.setSizePolicy(policy)
            icon = QIcon();
            icon.addPixmap(QPixmap(dirutils+os.sep+'Top_'+n+'.png'), QIcon.Normal, QIcon.Off)
            but.setIcon(icon);
            but.setIconSize(QSize(25, 25))
            but.setMaximumWidth(25)
            but.setMaximumHeight(25)
            but.setFlat(True)
            zoneSizer.addWidget(but)
            but.clicked.connect(self.clk)
        version = QLabel("       version 12/10/2021 ")
        zoneSizer.addWidget(version)
        #version.SetFont(wx.Font(7, wx.DEFAULT, wx.NORMAL, wx.NORMAL))
        self.obs = Observer()
        
    def record(self,name):
        self.obs.obs_value=name # to record a movement
        
    def clk(self):
        obj = self.ModifZone.sender()
        n = obj.objectName();
        if n == 'move': self.onMoveZ()
        if n == 'modifPoly': self.onModifZone()
        if n == 'supprime': self.onDelZone()
        if n == 'supprimeAll': self.onDelAllZones()
  
    def setChoiceList(self,obj,l):
        obj.clear()
        obj.addItems(l)
    
    def updateChoiceZone(self,line):
        # mise a jour de la liste de zone pour la nouvelle variable selectionnee
        # sur un milieu donne
        dicz = self.core.diczone[self.gui.currentModel].dic
        if (line == None) or (line not in list(dicz.keys())): 
            self.setChoiceList(self.choice,[''])
            self.currentZlist = None #EV 14/08/19
            self.valZ.setText("Val : ") #EV 14/08/19
            return
        #self.valZ.SetLabel(line+' :')
        self.currentZname = None 
        self.currentZlist = dicz[line]
        namelist = self.currentZlist['name']
        self.setChoiceList(self.choice,namelist)
        self.valZ.setText("Val : ")
   
    def onChoice(self, evt):  #choice of the  zone
        self.izone = self.choice.currentIndex();
        if self.currentZlist: #EV 14/08/19
            self.currentZname = self.currentZlist['name'][self.izone]
            self.valZ.setText("Val : "+str(self.currentZlist['value'][self.izone])[:5])
            #self.valZ.Enable(True)
        
    def onValueZ(self, evt):
        self.izone = self.choice.currentIndex(); #EV 14/08/19
        if self.currentZlist: #EV 14/08/19
            dialg = zoneDialog(self, self.core,self.gui.currentModel,self.gui.currentLine,self.currentZlist,self.izone)
            result = dialg.saveCurrent() #exec_() #show()
            if result !=None:
                line = self.gui.currentLine;#onMessage(self.gui,str(self.currentZlist['value']))
                val = self.currentZlist['value'][self.izone];#print val
                self.valZ.setText("Val : "+str(val)[:5])
                med = self.currentZlist['media'][self.izone]
                xy = self.currentZlist['coords'][self.izone]
                self.visu.modifZoneAttr(line, self.izone,val,med,xy)
                #self.modifZones(line)
                self.core.makeTtable()  # oa added 16/6 for isntant modif of values for isntnatfit
                self.record('zoneEnd') # oa 16/6 for the same reason
        
    def txt2coords(self,s):
        s1, xy = s.split('\n'),[]
        for s2 in s1:
            if len(s2)>1: 
                exec('a=['+s2+']')
                xy.append(a)
        return xy
        
    def onTable(self):
        tbl = self.core.diczone[self.gui.currentModel].dic.getTableOfZones(line)
    
    def onMoveZ(self):
        """ start moving a zone"""
        if self.currentZname != None:
            #self.gui.actions('zoneStart') 
            self.visu.startMoveZone(self.gui.currentLine, self.izone)
            #self.modifZones(self.gui.currentLine)
        else :
            onMessage (self.gui,"Please select a zone")
        
    def onModifZone(self):
        """" start modification of zone after tests"""
        if self.currentZname != None:
            #self.gui.actions('zoneStart') 
            self.visu.modifZone(self.gui.currentLine, self.izone)
            #self.modifZones(self.gui.currentLine)
        else :
            onMessage (self.gui,"Please select a zone")
            
    # get back coordinates and change them
    def onModifZoneCoord(self, line, index, coord):
        #self.gui.actions('zoneEnd')
        mod = self.gui.currentModel
        self.core.diczone[mod].setValue(self.gui.currentLine,'coords',self.izone,coord)        
        self.record('zoneEnd')
        
    def onDelZone(self):
        # si pas de zone selectionnee, on ne supprime pas...
        mod, line = self.gui.currentModel, self.gui.currentLine
        if self.currentZname != None:
            self.core.diczone[mod].delZone(line,self.izone)
            self.visu.delZone(line,self.izone)
            self.updateChoiceZone(line)
        else :
            onMessage(self.gui,"Please select a zone")

    def onDelAllZones(self):
        znb = len(self.currentZlist['name'])
        ok = onQuestion(self.gui,'Caution you will destroy all zones')
        if ok : #retour==wx.ID_OK:
            for i in range(znb-1,-1,-1):
                self.core.diczone[self.gui.currentModel].delZone(self.gui.currentLine,i)
                self.visu.delZone(self.gui.currentLine,i)
