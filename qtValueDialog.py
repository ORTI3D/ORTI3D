# -*- coding: utf-8 -*-
"""
Created on Sun Aug 02 10:22:05 2015

@author: olive
"""
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
        
class qtValueDialog(QDialog):
    def __init__(self, parent,gui,core,modName):
        QDialog.__init__(self)
        self.parent = parent
        self.setWindowTitle(modName+' parameters')
        layoutWidget = QWidget(self)
        screenShape = QDesktopWidget().screenGeometry()
        layoutWidget.setGeometry(QRect(5, 5, screenShape.width()*.2,400))
        self.vbox = QVBoxLayout(layoutWidget)
        self.vbox.setGeometry(QRect(0, 0, 80,60))

        grpList = core.getUsedModulesList(modName)
        self.chgroups = QComboBox(layoutWidget)
        for i,n in enumerate(grpList):
            self.chgroups.addItem("")
            self.chgroups.setItemText(i, n)
        self.chgroups.activated['QString'].connect(self.onChoiceGroup)
        self.vbox.addWidget(self.chgroups)
        
        self.chlines = QComboBox(layoutWidget)
        self.chlines.activated['QString'].connect(self.onChoiceLine)
        self.vbox.addWidget(self.chlines)
        
        self.boxkeys = qtBoxKeys(self,parent)
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
    def __init__(self,Main,parent):
        self.Main,self.parent = Main,parent
        self.layoutWidget = QWidget(Main)
        self.gridLayout = QGridLayout(self.layoutWidget) #EV 22/07/2019
        self.labl,self.lValBut,self.values=[],[],[]
        
    def setVisible(self,bool):self.layoutWidget.setVisible(bool)
        
    def addButtons(self,names,values,details,types):
        self.values,self.types = values,types;
        nb=len(names);
        # clear the layout
        for b in self.labl: b.deleteLater()
        for b in self.lValBut: b.deleteLater()
        self.labl,self.lValBut=[],[];
        #self.parent.gui.onMessage(str(names)+' '+str(values)+' '+str(details))
        
        for i in range(nb):
            bname,bcontent,bselect,btype=self.parent.makeButton(names[i],values[i],details[i],types[i])
            #print(bname,bcontent,bselect,btype)
            txt = QLabel(self.layoutWidget)
            txt.setText(bname)
            if btype == 'choice': 
                but = QComboBox(self.layoutWidget)
                but.addItems(bcontent)
                but.setCurrentIndex(bselect)
            #elif btype == 'text':
            else :
                but = QLineEdit(self.layoutWidget)
                but.setText(str(bcontent))
            if names[i].split('(')[0] in self.parent.blind: 
                but.setDisabled(True)        
            self.labl.append(txt)
            self.gridLayout.addWidget(txt,i,0,1,1)
            self.lValBut.append(but)
            self.gridLayout.addWidget(but,i,1,1,1)
            
    def getValues(self):
        nb = len(self.values)
        for i in range(nb):
            but = self.lValBut[i]
            val = self.values[i]
            #print but,val,self.labl[i],self.types[i]
            if self.types[i] in ['choice','layint']: #,'layint']: OA 6/11/18 removed layint
                self.values[i] = but.currentIndex()
                continue
            if but.text() not in ['formula','zone','array']:
                if self.types[i] in ['int','vecint','arrint']:
                    val=int(but.text())
                elif self.types[i] in ['float','vecfloat','arrfloat']:
                    val=float(but.text())
                else :  
                    val=str(but.text())  # OA 10/9/18 added str
            else :  
                val=str(but.text())  # OA 10/9/18 added str
            self.values[i]=val*1
        return self.values
