# -*- coding: utf-8 -*-
"""
Created on Fri Feb 28 13:39:43 2020

@author: asus
"""

from matplotlib.backends.qt_compat import QtCore, QtWidgets
from matplotlib.backends.backend_qt5agg import (
        FigureCanvas, NavigationToolbar2QT as NavigationToolbar)
from matplotlib.figure import Figure
from .qtDialogs import *
from .geometry import *
import numpy as np
import matplotlib.ticker as ticker

class mybudget(QDialog):
    '''This dialog provides plot for budget. There are 2 types of graphs:
    - Mass balance graphs: budgets through the various sources and sinks
    - Zone budget graphs: budgets for user-defined zones 
    The result that can be represented are mass of Flow, Transport, and 
    Chemistry.
    '''
    def __init__(self,gui,core,typ,res):
        self.gui,self.core= gui,core
        self.typ,self.res= typ,res
        QDialog.__init__(self,gui) 
        self.setModal(False)
        self.setWindowTitle('Plot of results')
        screenShape = QtWidgets.QDesktopWidget().screenGeometry()
        self.setGeometry(QRect(5, 5, screenShape.width()*.75, screenShape.height()*.7))
    ## main horizontal layout
        self.horizontalLayout = QHBoxLayout(self)
        self.horizontalLayout.setContentsMargins(10, 20, 10, 10)
    ## the left panel vertical layout  
        self.verticalLayout = QVBoxLayout()
        #self.verticalLayout.setGeometry(QRect(5, 5, 250, screenShape.height()*.68))
    ## title
        if self.typ == 'M' : title = 'Mass balance Graphs'
        if self.typ == 'Z' : title = 'Zone budget Graphs'
        label = str(title +' - '+self.res)
        self.label = QtWidgets.QLabel(self)
        self.label.setMaximumSize(250, 24)
        self.label.setText(label)
        font = QFont()
        font.setPointSize(9)
        font.setBold(True)
        self.label.setFont(font)
        self.verticalLayout.addWidget(self.label, alignment=Qt.AlignHCenter)
    ## model time list
        self.tlist = self.core.getTlist2()
    ## frame 1
        self.frame = QtWidgets.QFrame(self)
        self.frame.setMaximumSize(QtCore.QSize(250, 35)) 
        self.gl = QGridLayout(self.frame)
    ## Different type of graph
        self.label_1 = QtWidgets.QLabel(self.frame)
        self.label_1.setText("Type of graph")
        self.gl.addWidget(self.label_1,0,0,1,1)
        self.plgroup = QComboBox(self)
        self.plgroup.addItems(['Percent Discrepency','Time Series','In-Out',
                               'Time Step'])
        self.plgroup.setCurrentIndex(0)
        self.plgroup.activated['QString'].connect(self.onTstep)
        self.gl.addWidget(self.plgroup,0,1,1,1)
        self.verticalLayout.addWidget(self.frame)
    ## frame 2
        self.frame2 = QtWidgets.QFrame(self)
        self.frame2.setMaximumSize(QtCore.QSize(250,35)) 
        self.gl2 = QGridLayout(self.frame2)
    ## Time for time step graph
        self.label_2 = QtWidgets.QLabel(self.frame2)
        self.label_2.setText("Time")
        self.gl2.addWidget(self.label_2,0,0,1,1)
        self.Tstep = QComboBox(self.frame2)
        self.Tstep.addItems([str(n) for n in self.tlist])
        self.Tstep.setCurrentIndex(0)
        self.gl2.addWidget(self.Tstep,0,1,1,1)
        self.verticalLayout.addWidget(self.frame2) 
        self.frame2.hide()
    ## Choice of zone to perform budget
        if self.typ == 'Z':
    ## frame 3
            self.frame3 = QtWidgets.QFrame(self)
            self.frame3.setMaximumSize(QtCore.QSize(250, 35)) 
            self.gl3 = QGridLayout(self.frame3)
            lzname=self.core.diczone['Observation'].dic['obs.1']['name']
            self.label_3 = QtWidgets.QLabel(self.frame3)
            self.label_3.setText("Zone budget zone")
            self.gl3.addWidget(self.label_3,0,0,1,1)
            self.zgroup = QComboBox(self.frame3)
            self.zgroup.addItems([str(n) for n in lzname])
            self.zgroup.setCurrentIndex(0)
            self.gl3.addWidget(self.zgroup,0,1,1,1)
            self.zgroup.activated['QString'].connect(self.updateChoices)
            self.verticalLayout.addWidget(self.frame3) 
    ## the options :need to go in the interface to search for zones and others
        self.hlayout=QHBoxLayout()
        dic=self.getChoices() #self.res,self.typ
        self.nb = myNoteBookCheck(self.gui,"Options",dic)
        self.hlayout.addWidget(self.nb)
        self.nb.layout.removeWidget(self.nb.buttonBox) 
        self.nb.buttonBox.deleteLater()
        del self.nb.buttonBox
        self.verticalLayout.addLayout(self.hlayout)
        self.nb.apply()
    ## Apply button
        self.pushButton = QPushButton(self)
        self.pushButton.setText('Apply')
        self.verticalLayout.addWidget(self.pushButton, alignment=Qt.AlignHCenter)
        self.pushButton.clicked.connect(self.buildGraph)
     ## add vertical layout   
        self.horizontalLayout.addLayout(self.verticalLayout)
     ## the right panel vertical layout
        self.verticalLayout2 = QVBoxLayout()
     ## the matplotlib figure 
        self.figure = Figure(tight_layout=True,figsize=(7.8, 3), dpi=100) # EV 04/02/20 
        self.cnv = FigureCanvas(self.figure) 
        #self._ax = self.cnv.figure.subplots()#.add_axes([0.1, 0.15, 0.7, 0.8])
    ## add matplotlib figure
        #self.horizontalLayout.addWidget(self.cnv)
        self.verticalLayout2.addWidget(self.cnv)
        self.toolbar = NavigationToolbar(self.cnv, self) # EV 04/02/20 
        self.verticalLayout2.addWidget(self.toolbar, alignment=Qt.AlignHCenter) # EV 04/02/20 
    ## Export button
        self.pushButton2 = QPushButton(self)
        self.pushButton2.setText('Export')
        self.verticalLayout2.addWidget(self.pushButton2, alignment=Qt.AlignHCenter)
        #self.pushButton2.clicked.connect(self.onExport)
    ## add vertical layout 2
        self.horizontalLayout.addLayout(self.verticalLayout2)
        QMetaObject.connectSlotsByName(self)  #OA 1/6/19
        
    def onTstep(self):
    ## show or hide time combobox for the time step graph
        if self.plgroup.currentIndex()==3:
            self.frame2.show()
        else : self.frame2.hide()

    def getObsZone(self):
    ## get the names of model observation zone
        zname=self.core.diczone['Observation'].dic['obs.1']['name']
        zone = self.zgroup.currentText()
        if zone in zname:
            dicZin={'Zones':{}};dicZout={'Zones':{}}
            zname.remove(zone)
            zIn=[zone+' from '+zname[i] for i in range(len(zname))]
            dicZin['Zones'] = list(zip(zIn,[False]*len(zIn)))
            zOut=[zname[i]+' to '+zone for i in range(len(zname))]
            dicZout['Zones'] = list(zip(zOut,[False]*len(zOut)))
        return dicZin, dicZout
    
    def getBoundaries(self):
    ## get the boundaries in the model (well, drn, storage, chd...)
        dic={'Bound':{}} ; nbound=['STORAGE','CONSTANT HEAD']
        boundTyp=['WELLS','DRAINS','RIVER LEAKAGE','ET','HEAD DEP BOUNDS',
                  'RECHARGE']
        pack=['wel','drn','riv','evt','ghb','rch']
        for i, n in enumerate (pack) :
            if self.core.diczone['Modflow'].getNbZones(n+'.1')>0:
                nbound.append(boundTyp[i])
        dic['Bound'] = list(zip(nbound,[False]*len(nbound)))
        return dic
        
    def getSpecies(self):
    ## get the names of model chemical species
        dic={'Species':{}} 
        species = self.core.addin.chem.getListSpecies() 
        dic['Species']=list(zip(species,[False]*len(species)))
        return dic
    
    def getChoices(self):
    ## return a dic in function of type of graph and result to plot
        self.dicBound = self.getBoundaries()
        dicIn={'In':{}}  ; dicOut={'Out':{}} 
        if self.typ == 'M': 
            dicIn['In']=self.dicBound['Bound']
            dicOut['Out']=self.dicBound['Bound']
        else :
            dicZin,dicZout=self.getObsZone()
            dicIn['In']=self.dicBound['Bound']+dicZin['Zones']
            dicOut['Out']=self.dicBound['Bound']+dicZout['Zones']
        if self.res =='Chemistry': 
            dicSpecies=self.getSpecies()
            dic = {**dicSpecies,**dicIn,**dicOut} 
        else : dic = {**dicIn,**dicOut} 
        return dic
    
    def updateChoices(self):
        dicZin,dicZout=self.getObsZone()
        print('dicZin',dicZin,'dicZout',dicZout)
        dicIn={'In':{}}  ; dicOut={'Out':{}} 
        dicIn['In']=self.dicBound['Bound']+dicZin['Zones']
        dicOut['Out']=self.dicBound['Bound']+dicZout['Zones']
        if self.res =='Chemistry': 
            dicSpecies=self.getSpecies()
            dic = {**dicSpecies,**dicIn,**dicOut} 
        else : dic = {**dicIn,**dicOut} 
        self.nb = myNoteBookCheck(self.gui,"Options",dic)
        self.nb.update()
        self.hlayout.addWidget(self.nb)
        self.nb.layout.removeWidget(self.nb.buttonBox) 
        self.nb.buttonBox.deleteLater()
        del self.nb.buttonBox
        #self.verticalLayout.addWidget(self.nb)
        self.nb.apply()
    
    def getValues(self):
        for k in list(self.nb.dicIn.keys()):
            if k in list(self.nb.pages.keys()):
                names,boo = list(zip(*self.nb.dicIn[k]))
                lout = []
                items = self.nb.dwidget[k]
                for item in items:
                    lout.append(item.checkState())
                self.nb.dicOut[k] = list(zip(names,lout))
        return self.nb.dicOut
    
    def getOptions(self):
        '''get the plot options from the window very simple now'''
        dicIn={'ptyp':{},'graph':{},'inlist':{},'outlist':{},'splist':{}} # 'lylist':{}
        dicIn['ptyp']=self.typ
        dicIn['graph']=str(self.plgroup.currentText())
        dic=self.getValues() 
        dicIn['inlist']=[dic['In'][i][0] for i in range(len(dic['In'])) if dic['In'][i][1]==2]
        dicIn['outlist']=[dic['Out'][i][0] for i in range(len(dic['Out'])) if dic['Out'][i][1]==2]
        dicIn['splist']=[self.res]
        if self.res=='Chemistry' :
            dicIn['splist']=[dic['Species'][i][0] for i in range(len(dic['Species'])) if dic['Species'][i][1]==2] 
        return dicIn
    
    def buildGraph(self):
        dicIn=self.getOptions()
        print('dicIn',dicIn)