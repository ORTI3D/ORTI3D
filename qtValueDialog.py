# -*- coding: utf-8 -*-
"""
Created on Sun Aug 02 10:22:05 2015

@author: olive
"""
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from qtDialogs import * # OA 18/9/19
from geometry import * # OA 20/3/20

# try SPH dec 2023 adapted ftm to fitin this file
#from .trySPH import *
from defSPH import *
from configSPH import *

#############################################################
#############################################################
#trySPH part
#import datetime as dt

#trySPH_enableTrySPH = True
#enableTrySPH = False

#if trySPH_enableTrySPH:
#    trySPH_saveList = []

# def trySPH_addSave(values=None,saveList=None):
#     current = {'values': values,
#                'date': dt.datetime.now().strftime('%Y%m%d'),
#                'time': dt.datetime.now().strftime('%H:%M:%S'),
#                'lastAction': 'debug: always the same string...'}
#     saveList.append(current)
#
# def trySPH_loadSave(number=None):
#     pass
################################################
################################################

class qtValueDialog(QDialog):
    def __init__(self, parent,gui,core,modName):
        QDialog.__init__(self)
        self.parent = parent

        # modif canceled, the 'gui' parameter is finally used below
        # add SPH: before this 'gui' parameter was unused (and so useless)
        # but it turns out it is convenient to propagate the 'currentModel' attribute (which is in 'gui' but not in 'parent')
        #self.gui = gui

        self.setWindowTitle(modName+' parameters')
        self.layoutWidget = QWidget(self)
        screenShape = QDesktopWidget().screenGeometry()
        self.layoutWidget.setGeometry(QRect(5, 5, screenShape.width()*.2,screenShape.height()*.5))
        self.vbox = QVBoxLayout(self.layoutWidget)
        #self.vbox.setGeometry(QRect(0, 0, 80,60))

        grpList = core.getUsedModulesList(modName)
        self.chgroups = QComboBox(self.layoutWidget)
        for i,n in enumerate(grpList):
            self.chgroups.addItem("")
            self.chgroups.setItemText(i, n)
        self.chgroups.activated['QString'].connect(self.onChoiceGroup)
        self.vbox.addWidget(self.chgroups)
        
        self.chlines = QComboBox(self.layoutWidget)
        self.chlines.activated['QString'].connect(self.onChoiceLine)
        self.vbox.addWidget(self.chlines)

        # modif SPH to propagate the model
        # #print('debug SPH 240118:',self.currentModel)
        self.boxkeys = qtBoxKeys(self,parent,gui.currentModel)# line after the modif


        self.vbox.addWidget(self.boxkeys.layoutWidget)
        self.vbox.addStretch(1)
        #bBox = QDialogButtonBox(layoutWidget)
        #bBox.setOrientation(Qt.Horizontal)
        #bBox.setStandardButtons(QDialogButtonBox.Cancel|QDialogButtonBox.Ok)
        #QObject.connect(bBox, SIGNAL("accepted()"), parent.OnSetNewVal)
        self.bBox = QHBoxLayout() #EV 22/07/2019
        butApply = QPushButton('Apply')
        butApply.clicked.connect(parent.OnSetNewVal)
        self.bBox.addWidget(butApply)
        butClose = QPushButton('Close')
        butClose.clicked.connect(self.close)
        self.bBox.addWidget(butClose)
        self.vbox.addLayout(self.bBox)
        
        QMetaObject.connectSlotsByName(self)
        
    def onChoiceGroup(self,value): self.parent.onChoiceGroup(str(value))
    def onChoiceLine(self,value): self.parent.onChoiceLine(str(value))   
     
class qtBoxKeys:

    def __init__(self, Main, parent,currentModel):
        self.currentModel = currentModel
    # end modif SPH

        self.Main,self.parent = Main,parent
        self.layoutWidget = QWidget(Main)
        self.gridLayout = QGridLayout(self.layoutWidget) #EV 22/07/2019
        self.labl,self.lValBut,self.values=[],[],[]
        
    def setVisible(self,bool):self.layoutWidget.setVisible(bool)
        
    def addButtons(self,title,comm,names,values,details,types): #EV 25/09/19
        self.values,self.types,self.title,self.comm = values,types,title,comm;
        self.nb=len(names);
        # clear the layout
        for b in self.labl: b.deleteLater()
        for b in self.lValBut: b.deleteLater()
        self.labl,self.lValBut=[],[];
        #self.parent.gui.onMessage(str(names)+' '+str(values)+' '+str(details))
        
        for i in range(self.nb):
            bname,bcontent,bselect,btype=self.parent.makeButton(names[i],values[i],details[i],types[i])
            txt = QLabel(self.layoutWidget)
            txt.setText(bname)
            if btype in ['choice','laychoice']: 
                but = QComboBox(self.layoutWidget)
                but.addItems(bcontent)
                but.setCurrentIndex(bselect)
            elif btype in ['layint','layfloat']:  # OA modif  2/10/19
                but = QPushButton('Media',self.layoutWidget)
                but.clicked.connect(self.onOpenVectDialog)
                self.vect = values; # OA 18/9/19 values contain all layers for type lay..
            elif btype =='textlong':  # OA modif  21/11/19
                scrollArea = QScrollArea(self.layoutWidget)
                scrollArea.setGeometry(QRect(50, 50, 100, 50))
                but = QTextEdit(scrollArea)
                but.setText(str(bcontent))
            else :
                but = QLineEdit(self.layoutWidget)
                but.setText(str(bcontent))
            if names[i].split('(')[0] in self.parent.blind: 
                but.setDisabled(True)        
            self.labl.append(txt)
            self.gridLayout.addWidget(txt,i,0,1,1)
            self.lValBut.append(but)
            self.gridLayout.addWidget(but,i,1,1,1)
        screenShape = QDesktopWidget().screenGeometry() # OA 2/5/20, this and two lines below fro long list of buttons
        w,h = screenShape.width()*.2,screenShape.height()*(.2+.025*self.nb)
        self.Main.layoutWidget.setGeometry(QRect(5, 5,w,h))
            
    def onOpenVectDialog(self): # OA 17/9/19
        ''' creates a vector dialog in case of layint type '''
        nr = len(self.vect);#print('qtValDlg 98',self.vect)
        nmed = getNmedia(self.parent.core)#;print('qt val 100',self.vect,nmed)
        if nr<nmed : self.vect.extend([self.vect[-1]]*(nmed-nr)) # the nb of media is higher than vect
        elif nr>nmed: self.vect = self.vect[:nmed] # the contrary
        dicV = {self.comm:{'cols':['media','value'],'rows':['']*nmed,'data':[[i,a] for i,a in enumerate(self.vect)]}}; #EV 25/09/19
        dlgVect = myNoteBook(self.parent.core.gui,self.title,dicV) #EV 25/09/19
        dicout = dlgVect.getValues()
        if dicout != None:
            self.vect = [x[1] for x in dicout[self.comm]['data']]

    def getValues(self):

# note SPH january 2024: at this 'getValues' level, the value of self.currentModel is always 'Modflow' and is not relevant
#I finally use self.parent.currentLine & self.parent.currentGroup, but there is no self.parent.currentModel nor self.core
#moreover self.parent.currentGroup is not what the main GUI call group ('Modflow series', 'Modflow USG', 'Min3p', etc.) but
#is the currentLine family. With these values I can't distinct between 'Modflow series' & 'Modflow USG'

#it seems to me : the fact that self.currentModel is always 'Modflow' shows eiter a bug or uselessness


        if modifSPH2023 and modifSPH2023_useInspect:
            import inspect
            #curframe = inspect.currentframe()
            #calframe = inspect.getouterframes(curframe,index=2)
            print('debug 240115e:',inspect.getouterframes(inspect.currentframe()))
            #print(str(inspect.getframeinfo(inspect.currentframe().f_back)))# to get just filename, line_number, function_name, lines, index of just the caller

        #try SPH 231215: init save before any modif
        if modifSPH2023:
            import copy
        #     # canceled
        #     #sph_addSave(calframe=calframe,values=self.values,backupList='dev',mode='pickle')#trySPH_saveList)
            self.previousValues = copy.deepcopy(self.values)

        #nb = len(self.values) # OA 23/9/19 removed
        for i in range(self.nb):
            but = self.lValBut[i]
            val = self.values[i]
            if self.types[i] in ['choice','laychoice']: #OA 2/10/19 added laychoice
                self.values[i] = but.currentIndex()
                continue
            elif self.types[i] in ['layint','layfloat']: # OA added 23/9/19, EV 25/09/19 added layfloat
                self.values = self.vect
                continue
            elif self.types[i]=='textlong': # added 21/11/19
                self.values = [but.document().toPlainText()]
                continue
            if but.text() not in ['formula','zone','array']:
                if self.types[i] in ['int','vecint','arrint']:
                    val=int(but.text())
                elif self.types[i] in ['float','vecfloat','arrfloat']:
                    val=float(but.text())
                else :  
                    val=str(but.text())  # OA 10/9/18 added str
            else :  
                val=str(but.text());#print(val)  # OA 10/9/18 added str
            self.values[i]=val*1

        # try SPH 231215: save after any modif of the values
        # self.sph_addSave()
 #       if modifSPH2023_debugLevel>5:
# the value in self.currentModel is ambiguous ! need to find the currentGroup
#                 print('debug 240129: at the getValues level, currentModel is',self.currentModel)
#                 print('    and currentLine is',self.parent.currentLine)
#                 print('    and currentGroup is',self.parent.currentGroup)
            #for key in trySPH_modelsAndTerms:
            #    if self.currentModel

            # TODO: ADD WARNING
            # TODO: THE pickle LOADING PROCESS IS NOT COMPLETELY SECURED ACCORDING TO INTERNET COMMENTS !!!
            # sph_addSave(currentLine=self.parent.currentLine,
            #             currentGroup=self.currentModel,
            #             initialValues=initialValues,
            #             finalValues=self.values,#backupList='dev',
            #             mode='pickle')#trySPH_saveList)

        return self.values

