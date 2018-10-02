# -*- coding: utf-8 -*-


from PyQt5.QtCore import *
from PyQt5.QtGui import *

import os
from .parameters import BaseParms
from .qtDialogs import *
from .core import *
from .config import *

class Ui_Parameters(object):
    def setupUi(self, Parameters,gui,core,plugin_dir):
        self.gui,self.core,self.plugin_dir =gui, core,plugin_dir
        self.base = BaseParms(gui,core)
        Parameters.setObjectName("Parameters")
        #Parameters.resize(150, 500)
        Parameters.setWindowTitle( "Parameters")
        self.dictBox={}
        skey = list(self.base.groups.keys()); skey.sort()
        for i,g in enumerate(skey): 
            self.dictBox[g] = Box(Parameters,self,g,i)
        QMetaObject.connectSlotsByName(Parameters)
        
class Box(): #QGroupBox):
    def __init__(self,Parameters,parent,gr,nb):
        '''parent is the Ui_parameters class above'''
        self.box = QGroupBox(Parameters)
        self.screenShape = QDesktopWidget().screenGeometry()
        self.box.setGeometry(QRect(0, nb*70, self.screenShape.width()*0.08, 60))
        self.box.setMaximumWidth(145)
        self.Parameters,self.parent = Parameters,parent
        self.box.setTitle(gr)
        self.hl = QHBoxLayout(self.box)
        self.hl.setSpacing(0)
        self.hl.setContentsMargins(1,1,1,1)
        dirutils = parent.plugin_dir+os.sep+'utils'
        #self.parent.gui.dialogs.onMessage(self.parent.gui,os.listdir(dirutils)[0])
        #policy = QSizePolicy(QSizePolicy.Preferred, QSizePolicy.Preferred)

        butA = parent.core.addin.addButton(self,gr) # a list of buttons
        if butA !=None:
            #print os.listdir(dirutils)
            for short,name,pos in butA : 
                if pos==1 : continue
                buta = QPushButton()
                shortName = name+'.png'#;print 'a', shortName
                if shortName in os.listdir(dirutils):
                    icon = QIcon()
                    icon.addPixmap(QPixmap(dirutils+os.sep+shortName), QIcon.Normal, QIcon.Off)
                    #buta.setSizePolicy(policy)
                    buta.setIcon(icon)
                    buta.setIconSize(QSize(25,25))
                    buta.setMaximumWidth(25)
                    buta.setMaximumHeight(25)
                    buta.setFlat(True)
                else :
                    buta.setText(short)
                buta.setObjectName(name)
                buta.setToolTip(name)
                self.hl.addWidget(buta)
                buta.clicked.connect(self.onButton)
                parent.base.dicaction[name] = 'self.addin.doaction(\''+name+'\')'

        for i in range(len(parent.base.groups[gr])):
            n=parent.base.groups[gr][i]
            shortName = gr[2:4]+'_'+n;#print i,shortName
            but = QPushButton(self.box) #hlWidget)
            but.setToolTip(n)
            icon = QIcon()
            icon.addPixmap(QPixmap(dirutils+os.sep+shortName+'.png'), QIcon.Normal, QIcon.Off)
            #but.setSizePolicy(policy)
            but.setIcon(icon)
            but.setIconSize(QSize(25, 25))
            but.setMaximumWidth(25)
            but.setMaximumHeight(25)
            but.setFlat(True)
            but.setObjectName(shortName) #_fromUtf8(n))
            but.clicked.connect(self.onButton)
            self.hl.addWidget(but)

        if butA !=None:
            for short,name,pos in butA : 
                if pos==0 : continue
                buta = QPushButton(self.box) #hlWidget)
                shortName = name+'.png'
                if shortName in os.listdir(dirutils):
                    icon = QIcon()
                    icon.addPixmap(QPixmap(dirutils+os.sep+shortName), QIcon.Normal, QIcon.Off)
                    buta.setIcon(icon)
                    buta.setIconSize(QSize(25,25))
                    buta.setFixedWidth(25)
                    buta.setFixedHeight(25)
                    buta.setFlat(True)
                else :
                    buta.setText(short)
                buta.setObjectName(name)
                buta.setToolTip(name)
                self.hl.addWidget(buta)
                buta.clicked.connect(self.onButton)
                parent.base.dicaction[name] = 'self.addin.doaction(\''+name+'\')'
            
    def onButton(self):
        s = self.Parameters.sender()
        name = s.objectName()
        self.parent.base.action(name)

