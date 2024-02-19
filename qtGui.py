# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'ui_file.ui'
#
# Created: Sat Feb 15 15:09:21 2014
#      by: PyQt5 UI code generator 4.8.6



from PyQt5.QtCore import *
from PyQt5.QtGui import *

import os
from parameters import BaseParms
from menus import *
from qtDialogs import *
from geometry import *
from qtVisu import *
from core import *
from config import *
from guiShow import guiShow

class Ui_Main(object):
    def setupUi(self,Main,iface,plugin_dir):
        Main.setObjectName("Main")
        Main.resize(200,560)
        Main.setWindowTitle("qORTi3d")
        self.gui = QWidget(Main)
        self.gui.setGeometry(QRect(5, 15, 195, 480)) #left,top,w,h
        self.gui.gtyp,self.gui.plugin_dir = 'qt',plugin_dir
        core = Core(self.gui)
        cfg = Config(core)
        self.dialogs = cfg.dialogs
        self.gui.linesDic,self.gui.linesCommDic = {},{}
        for mod in core.modelList:
            self.gui.linesDic[mod] = core.diczone[mod].getLinesDic()
            self.gui.linesCommDic[mod] = core.diczone[mod].getLinesCommDic()
        self.gui.on3D,self.gui.onRCT,self.gui.onSetMediaNb = self.on3D,self.onRCT,self.onSetMediaNb
        
        core.addin.setGui(self.gui)
        self.gui.visu = qtVisu(iface,self.gui,core)
        
        self.toolBox = QToolBox(self.gui)
        self.toolBox.setGeometry(QRect(0, 0, 195, 480))

        self.page0 = QWidget()        
        self.pFile = Ui_File()
        self.pFile.setupUi(self.page0,self.gui,core)
        self.toolBox.addItem(self.page0,"Files and tools")
        
        self.page1 = QWidget()
        self.pParameters = Ui_Parameters()
        self.pParameters.setupUi(self.page1,self.gui,core,plugin_dir)
        self.toolBox.addItem(self.page1,"Model parameters")
        
        self.page2 = QWidget()
        self.gui.varBox = Ui_Var()
        self.gui.varBox.setupUi(self.page2,self.gui,core)
        self.toolBox.addItem(self.page2,"Spatial variables")

        self.page3 = QWidget()
        self.gui.dlgShow = Ui_Show()
        self.gui.dlgShow.setupUi(self.page3,self.gui,core)
        self.toolBox.addItem(self.page3,"Results")

        QMetaObject.connectSlotsByName(Main)
            
    def on3D(self,dm): pass
    def onRCT(self,bool): pass # for compatibility with wx
    def onSetMediaNb(self,nbM,nbL): pass
        
class Ui_File(object):
    def setupUi(self,File,gui,core):
        self.menus = Menus(gui,core)
#        File.setObjectName("File")
        self.File,self.gui,self.core = File,gui,core
        
        File.resize(177, 104)
        self.gridLayoutWidget = QWidget(File)
        self.gridLayoutWidget.setGeometry(QRect(0, 0, 180, 120))
        self.gridLayout = QGridLayout(self.gridLayoutWidget)
        self.gridLayout.setMargin(0)

        llabel=['File','Import','Export','Add_in','?']
        litems=[["","New","Open","Save", "Save as"],["","Solutions","User Species"],
                ["","Current variable"],[],["","Help","Download new v.","Back to old v."]]
        self.combo= []
        for i,lab in enumerate(llabel):
            label = QLabel(self.gridLayoutWidget)
            label.setText(lab)
            self.gridLayout.addWidget(label, i, 0, 1, 1)
            self.combo.append(QComboBox(self.gridLayoutWidget))
            self.combo[i].addItems(litems[i])
            self.gridLayout.addWidget(self.combo[i], i, 1, 1, 1)
            self.combo[i].activated['QString'].connect(self.onClick)
        
        QMetaObject.connectSlotsByName(File)
        #self.core.openModel('E://ipht3d//exv2//test','iqexer1')

    def onClick(self):
        obj = self.File.sender()
        n = obj.currentText();
        #self.gui.dialogs.onMessage(self.gui,n)
        if n=="New": self.menus.OnNew()
        if n=="Open": self.menus.OnOpen()
        if n=="Save": self.menus.OnSave()
        if n=="Save as": self.menus.OnSaveAs()
        if n=='Solutions': self.menus.OnImportSolutions()
        if n=='User Species': self.menus.OnImportUserSpecies()
        if n=='Current variable': self.menus.OnExportVar()
        if n=='Help': self.menus.OnHelp()
        if n=='Download new v.':self.menus.OnDownload()
        if n=='Back to old v.': self.menus.OnBackVersion()

class Ui_Parameters(object):
    def setupUi(self, Parameters,gui,core,plugin_dir):
        self.gui,self.core,self.plugin_dir =gui, core,plugin_dir
        self.base = BaseParms(gui,core)
        Parameters.setObjectName("Parameters")
        Parameters.resize(197, 348)
        Parameters.setWindowTitle( "Parameters")
        self.dictBox={}
        skey = list(self.base.groups.keys()); skey.sort()
        for i,g in enumerate(skey): 
            self.dictBox[g] = Box(Parameters,self,g,i)
        QMetaObject.connectSlotsByName(Parameters)
        
class Box:
    def __init__(self,Parameters,parent,gr,nb):
        '''parent is the Ui_parameters class above'''
        self.box = QGroupBox(Parameters)
        self.Parameters,self.parent = Parameters,parent
        y0=20+nb*60
        self.box.setGeometry(QRect(5, y0, 170, y0+40))
        self.box.setTitle(gr)
        self.hlWidget = QWidget(self.box)
        self.hlWidget.setGeometry(QRect(9, 15, 158, 28))
        self.hl = QHBoxLayout(self.hlWidget)
        self.hl.setMargin(0)
        dirutils = parent.plugin_dir+os.sep+'utils'
        #self.parent.gui.dialogs.onMessage(self.parent.gui,os.listdir(dirutils)[0])

        butA = parent.core.addin.addButton(self,gr) # a list of buttons
        if butA !=None:
            for short,name,pos in butA : 
                if pos==1 : continue
                buta = QPushButton(self.hlWidget)
                shortName = 'Ad_'+short+'.gif'
                if shortName in os.listdir(dirutils):
                    icon = QIcon()
                    icon.addPixmap(QPixmap(dirutils+os.sep+shortName), QIcon.Normal, QIcon.Off)
                    buta.setIcon(icon)
                else :
                    buta.setText(short)
                buta.setObjectName(name)
                self.hl.addWidget(buta)
                buta.clicked.connect(self.onButton)
                parent.base.dicaction[name] = 'self.addin.doaction(\''+name+'\')'

        for i in range(len(parent.base.groups[gr])):
            n=parent.base.groups[gr][i]
            shortName = gr[2:4]+'_'+n
            but = QPushButton(self.hlWidget)
            but.setToolTip(n)
            icon = QIcon()
            icon.addPixmap(QPixmap(dirutils+os.sep+shortName+'.gif'), QIcon.Normal, QIcon.Off)
            but.setIcon(icon)
            but.setObjectName(shortName) #_fromUtf8(n))
            but.clicked.connect(self.onButton)
            self.hl.addWidget(but)

        if butA !=None:
            for short,name,pos in butA : 
                if pos==0 : continue
                buta = QPushButton(self.hlWidget)
                shortName = 'Ad_'+short+'.gif'
                if shortName in os.listdir(dirutils):
                    icon = QIcon()
                    icon.addPixmap(QPixmap(dirutils+os.sep+shortName), QIcon.Normal, QIcon.Off)
                    buta.setIcon(icon)
                else :
                    buta.setText(short)
                buta.setObjectName(name)
                self.hl.addWidget(buta)
                buta.clicked.connect(self.onButton)
                parent.base.dicaction[name] = 'self.addin.doaction(\''+name+'\')'
            
    def onButton(self):
        s = self.Parameters.sender()
        name = s.objectName()
        self.parent.base.action(name)

class Ui_Var(object):
    def setupUi(self, Var,gui,core):
        self.gui,self.core = gui,core
        Var.setObjectName("Var")
        #Var.resize(177, 104)
        Var.setWindowTitle("Var")
        self.gridLayoutWidget = QWidget(Var)
        self.gridLayoutWidget.setGeometry(QRect(0, 0, 180, 180))
        #self.gridLayoutWidget.setObjectName(_fromUtf8("gridLayoutWidget"))
        self.gridLayout = QGridLayout(self.gridLayoutWidget)
        self.gridLayout.setMargin(0)
        #self.gridLayout.setObjectName(_fromUtf8("gridLayout"))

        categories = ['Modflow series','Min3p','OpGeo']
        label = QLabel(self.gridLayoutWidget)
        label.setText("Category")
        self.gridLayout.addWidget(label, 0, 0, 1, 1)
        self.choiceC = QComboBox(self.gridLayoutWidget)
        for c in categories: self.choiceC.addItem(c)
        self.choiceC.activated['QString'].connect(self.onChoiceCategory)
        self.gridLayout.addWidget(self.choiceC, 0, 1, 1, 1)

        label = QLabel(self.gridLayoutWidget)
        label.setText("Model")
        self.gridLayout.addWidget(label, 1, 0, 1, 1)
        self.choiceM = QComboBox(self.gridLayoutWidget)
        self.choiceM.activated['QString'].connect(self.onChoiceModel)
        self.gridLayout.addWidget(self.choiceM, 1, 1, 1, 1)

        label = QLabel(self.gridLayoutWidget)
        label.setText("Group")
        self.gridLayout.addWidget(label, 2, 0, 1, 1)
        self.choiceG = QComboBox(self.gridLayoutWidget)
        self.choiceG.activated['QString'].connect(self.onChoiceGroup)
        self.gridLayout.addWidget(self.choiceG, 2, 1, 1, 1)

        label = QLabel(self.gridLayoutWidget)
        label.setText("Line")
        self.gridLayout.addWidget(label, 3, 0, 1, 1)
        self.choiceL = QComboBox(self.gridLayoutWidget)
        self.choiceL.activated['QString'].connect(self.onChoiceLine)
        self.gridLayout.addWidget(self.choiceL, 3, 1, 1, 1)

        label = QLabel(self.gridLayoutWidget)
        label.setText("Media")
        self.gridLayout.addWidget(label, 4, 0, 1, 1)
        self.choice3D = QComboBox(self.gridLayoutWidget)
        self.choice3D.activated['QString'].connect(self.onChoiceMedia)
        self.gridLayout.addWidget(self.choice3D, 4, 1, 1, 1)

        label = QLabel(self.gridLayoutWidget)
        label.setText("Backg.")
        self.gridLayout.addWidget(label, 5, 0, 1, 1)
        self.backg = QLineEdit(self.gridLayoutWidget)
        self.backg.returnPressed.connect(self.onBackOk)
        self.gridLayout.addWidget(self.backg, 5, 1, 1, 1)

        label = QLabel(self.gridLayoutWidget)
        label.setText("Type")
        self.gridLayout.addWidget(label, 6, 0, 1, 1)
        self.choiceT = QComboBox(self.gridLayoutWidget)
        typeList = ['one_value','formula','zone','edit','interpolate','import']
        self.choiceT.addItems(typeList)
        self.choiceT.activated['QString'].connect(self.onChoiceType)
        self.gridLayout.addWidget(self.choiceT, 6, 1, 1, 1)
        
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

    def onChoiceCategory(self):
        """the categories of models : modflow series, sutra, min3p"""
        #evt.skip()
        c0 = self.choiceC.currentText()
        self.chooseCategory(c0)
        
    def chooseCategory(self,c0):
        self.currentCategory = c0;#print c0,self.core.modelList
        lshort = ['Sutra','Min3p','Fipy']
        lshort2 = ['Sutr','Min3','Fipy']
        if c0 in lshort:
            lmodels = [x for x in self.core.modelList if x[:4]==c0[:4]]
        else :
            lmodels = [x for x in self.core.modelList if x[:4] not in lshort2]
        if 'Observation' not in lmodels: lmodels.append('Observation')
        self.setChoiceList(self.choiceM,lmodels)

    def onChoiceModel(self,evt):
        """contains the models : for modflow series : modflow, mt3d, pht3d..."""
        m0 = self.choiceM.currentText()
        if m0[:4] not in ['Min3','Fipy','Sutr']: #modflow series
            m1 = m0.lower();mod = m1[0].upper()+m1[1:];
        else :
            mod = m0
        self.currentModel = mod
        lmodules = self.core.addin.getUsedModulesList(mod); #print 'topbar,l84',lmodules# to set the groups
        self.setChoiceList(self.choiceG,lmodules)
        
    def onChoiceGroup(self,evt):
        curGroup = self.choiceG.currentText()
        self.choiceL.clear()
        if curGroup not in list(self.gui.linesDic[self.currentModel].keys()) : return
        lines = self.gui.linesDic[self.currentModel][curGroup]
        indx = self.testConditions(self.currentModel,lines)
        lcomm = self.gui.linesCommDic[self.currentModel][curGroup]
        for il in indx: 
            if lines[il] in self.blind[self.currentModel]: continue
            self.choiceL.addItem(lines[il]+' '+lcomm[il])
        self.currentLine = lines[0]
        
    def testConditions(self,modName,lstL):
        """ test if the lines indices given in lstL sastify the condition"""
        indexout=[];
        for i,l in enumerate(lstL):
            cond = self.core.dickword[modName].lines[l]['cond']
            if self.core.testCondition(modName,cond):
                indexout.append(i)
        return indexout

    def onChoiceLine(self,evt):
        line = str(self.choiceL.currentText().split()[0])
        self.currentLine = line
        media = self.currentMedia
        vallist = self.core.dicval[self.currentModel][line];#print 'topbar valist',vallist
        nval = len(vallist)
        nmedia = getNmedia(self.core)
        if nval<nmedia : vallist.extend([vallist[0]]*(nmedia-nval))
        self.backg.setText(str(vallist[media])) # set to the new value
        typ = self.core.dictype[self.currentModel][line][0]
        i = self.choiceT.findText(typ)
        self.choiceT.setCurrentIndex(i)
            
    def onChoiceMedia(self,evt):
        """changes the media in visualization and stores the current media"""
        media = self.choice3D.currentIndex()
        self.currentMedia = media
        line = self.currentLine
        vallist = self.core.dicval[self.currentModel][line]
        self.backg.SetValue(str(vallist[media]))
        
    def onBackOk(self):
        line = self.currentLine
        media = self.currentMedia
        self.core.setValue(self.currentModel,line,media,float(self.backg.text()))
        
    def onChoiceType(self,evt):
        line = self.currentLine
        choice = str(self.choiceT.currentText()); #print choice
        self.core.dictype[self.currentModel][line] = [choice]
        if choice =='zone': self.gui.visu.createLayer(line)
        elif choice == 'import': self.onImport(evt)
        elif choice == 'edit': self.onEdit(evt)
        elif choice == 'formula': self.onFormula(evt)
        elif choice == 'importZones': 
            self.onImportZones(evt)
            self.core.dictype[self.currentModel][line] = ['zone']
        
    def onImportZones(self,evt):
        fdialg = myFileDialogOpen()
        fileDir,fileName = fdialg.getFile(self.gui,evt,'choose zone file','*.txt')
        self.core.importZones(fileDir,fileName,self.currentModel,self.currentLine)

    def onEdit(self,evt):
        pass

    def onFormula(self,evt):
        """opens a dialog to ask for python formula and executes them
        to get the value of the given keyword in the last line
        """
        ll = self.currentLine
        formula = self.core.getFormula(self.currentModel,ll)[0]; #print formula
        dialg = textDialog(self,'input python formula',(340,300),str(formula))
        retour = dialg.getText()
        if retour != None:
            formula = retour;
            self.core.dicformula[self.currentModel][ll]=[str(formula)]
            self.core.dictype[self.currentModel][ll]=['formula']

class Ui_Show(object):
    def setupUi(self,Show,gui,core):
        self.Show = Show
        self.gui,self.core = gui,core
        self.gui.guiShow = guiShow(gui,core)
        self.groups = self.gui.guiShow.groups
        Show.setObjectName("Show")
        self.dictBox={}
        pos=0
        for ig in range(len(self.groups)): 
            for g0 in self.groups: #pour ordonner
                if self.groups[g0][0]==ig: g=g0
            names = self.groups[g][1:]  
            self.dictBox[g] = showBox(Show,self.gui.guiShow,names,g,pos)
            pos += len(names)*24+20
        Show.resize(200,560)
        QMetaObject.connectSlotsByName(Show)
        
    def getCurrentTime(self):
        combo = self.findChild(QComboBox,'Aquifer_Tstep_L')
        return combo.getText()
        
    def getNames(self,nameBox):
        combo = self.Show.findChild(QComboBox,nameBox) 
        names = [combo.itemText(i) for i in range(combo.count())]
        return names
    def setNames(self,nameBox,names,opt='strings'):
        combo = self.Show.findChild(QComboBox,nameBox) 
        combo.clear() # danger if set it is not possible to add items
        combo.addItems([str(n) for n in names])

    def uncheckContours(self):
        """used to uncheck the other contours when group is changed"""
        dic = self.gui.guiShow.dicVisu
        for n,m in [('Flow','Head'),('Flow','Wcontent'),('Transport','Tracer')]:
            self.onTickBox(n,m,'B',dic[n][m])

    def onTickBox(self,group,name,tag,bool):
        """ to change the state of a button whithout doing any action"""
        item = self.Show.findChild(QCheckBox,group+'_'+name+'_'+tag);
        if tag=='B': item.setCheckState(bool)

class showBox:
    def __init__(self,Show,parent,names,g,pos):
        self.Show,self.parent = Show,parent
        self.group = QGroupBox(Show)
        self.group.setTitle(g)
        ln = len(names)
        self.group.setGeometry(QRect(0, pos, 195, pos+ln*24))
        self.hlWidget = QWidget(self.group)
        self.hlWidget.setGeometry(QRect(5,5 ,195, ln*24))
        boxGrid = QGridLayout(self.hlWidget)
        boxGrid.setMargin(1)
        self.buts = list(range(len(names)))
        for i,n in enumerate(names):
            if type(n)==type([1,2]): #cas liste -> choix
                text = QLabel(self.hlWidget)
                text.setText(n[0])
                boxGrid.addWidget(text,i,0,1,1)
                liste = n[1]
                self.buts[i] = QComboBox(self.hlWidget)
                self.buts[i].setObjectName(g+'_'+n[0]+'_L')
                self.buts[i].addItems(liste)
                boxGrid.addWidget(self.buts[i],i,1,1,1)
                self.buts[i].activated['QString'].connect(self.onClick)        
            else : # cas simple : checkbox
                text = QLabel(self.hlWidget)
                text.setText(n)
                boxGrid.addWidget(text,i,0,1,1)
                self.buts[i] = QCheckBox(self.hlWidget)
                self.buts[i].setObjectName(g+'_'+n+'_B')
                boxGrid.addWidget(self.buts[i],i,1,1,1)
                self.buts[i].stateChanged.connect(self.onClick)
                #chk.clicked.connect(self.onClick)

    def onClick(self,value):
        """action when a box is clicked, tag L : list """
        item = self.Show.sender()
        n = item.objectName(); 
        [group,name,tag]=n.split('_');#print 'guish onclick',group,name,tag
        if tag=='L': 
            if name in ['Layer','Tstep']: 
                retour = item.currentIndex()
            else :
                retour = item.currentText() # case of list, retour is the name
        else: retour = item.isChecked() # a check box retour is True or False 
        if name in self.parent.Vtypes['Array']: self.parent.resetDicContour()
        self.parent.dicVisu[group][name]=retour
        nz,ny,nx = shape(self.parent.core.Zblock)
        if name == 'Plane': 
            exec('self.parent.'+name+'=\"'+str(retour)+'\"')
            #self.changeIcOri(retour)
            if retour =='Z' : self.setNames('Model_Layer_L',list(range(nz-1)))
            if retour =='Y' : self.setNames('Model_Layer_L',list(range(ny-1)))
            if retour =='X' : self.setNames('Model_Layer_L',list(range(nx-1)))
            self.parent.dicVisu['Model']['Layer']=0
        self.parent.onClick2(group,name,retour)

    def setNames(self,nameBox,names,opt='strings'):
        self.parent.setNames(nameBox,names)