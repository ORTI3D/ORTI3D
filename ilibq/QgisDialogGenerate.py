# -*- coding: utf-8 -*-

from PyQt4 import QtCore, QtGui
import qgis.core as qgis

from qtDialogs import *
from myDialogs import *

import types
from qtVisu import *
from qgisGui import *
from menus import Menus

# Import the utilities from the fTools plugin (a standard QGIS plugin),
# which provide convenience functions for handling QGIS vector layers
import sys, os, imp, platform
import fTools
path = os.path.dirname(fTools.__file__)
ftu = imp.load_source('ftu', os.path.join(path,'tools','ftools_utils.py'))

class mainDialogGenerate(QtGui.QDialog):
    def __init__(self,parent, iface,core):
        self.menu = Menus(iface,core)
        QtGui.QDialog.__init__(self)
        self.setWindowFlags(QtCore.Qt.WindowStaysOnTopHint)
        self.parent, self.iface, self.core = parent, iface, core
        self.linesDic,self.linesCommDic = {},{}
        for mod in core.modelList:
            self.linesDic[mod] = core.diczone[mod].getLinesDic()
            self.linesCommDic[mod] = core.diczone[mod].getLinesCommDic()

        self.ui = Ui_Main()
        self.ui.setupUi(self,iface,core)
        
        self.ui.pFile.comboFile.activated['QString'].connect(self.fileChange)
        self.ui.pVar.comboLine.activated['QString'].connect(self.lineChange)
        self.ui.pVar.comboType.activated['QString'].connect(self.typeChange)
        
    def showDialogAndDisconnect(self):
        self.show()
        self.iface.actionToggleEditing().triggered.disconnect(self.showDialogAndDisconnect)
    
    def fileChange(self,value):
        if value=='Open': 
            self.menu.OnOpen(self)
            self.parent.visu.clearLayers()
            self.parent.visu.initDomain()
            self.parent.visu.qgsCore2Zones()
            # set the show results part 
            tl2 = self.core.getTlist2() # pb of shape of tl2
            self.iface.onMessage(str(tl2))
            self.ui.pShow.setNames('Aquifer_Tstep_L',tl2)
            listSpec = self.core.addin.pht3d.getListSpecies() # just the names
            self.ui.pShow.setChemSpecies(listSpec)
        if value=='Save':  
            self.parent.visu.qgsZones2Core()
            self.menu.OnSave(self)
            #self.ui.pFile.comboImport.clear()
            #self.ui.pFile.comboImport.insertItems(0,['c','b'])

    def lineChange(self,line):
        self.line = line[:5]
        typeList = ['one_value','formula','zone','edit','interpolate','import']
        typ = self.core.dictype[self.model][self.line][0]
        self.ui.pVar.comboType.setCurrentIndex(typeList.index(typ))

    def typeChange(self,typ):
        self.core.dictype[self.model][self.line][0] = typ
        lines = self.linesDic[self.model][self.group]
        il = lines.index(self.line)
        layerName = lines[il]
        if typ == 'zone': self.parent.visu.createLayer(layerName)

