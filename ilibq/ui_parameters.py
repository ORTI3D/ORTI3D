# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'ui_parameters.ui'
#
# Created: Sun Feb 16 10:03:57 2014
#      by: PyQt4 UI code generator 4.8.6
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    _fromUtf8 = lambda s: s
    
from parameters import BaseParms

class Ui_Parameters(object):
    def setupUi(self, Parameters,gui,core):
        self.gui,self.core =gui, core
        self.Base = BaseParms(gui,core)
        Parameters.setObjectName(_fromUtf8("Parameters"))
        Parameters.resize(197, 348)
        Parameters.setWindowTitle(QtGui.QApplication.translate("Parameters", "Parameters", None, QtGui.QApplication.UnicodeUTF8))
        self.dictBox={}
        skey = self.Base.groups.keys(); skey.sort()
        i=0
        for g in skey: 
            self.dictBox[g] = Box(Parameters,self,g,i);i+=1

        self.retranslateUi(Parameters)
        QtCore.QMetaObject.connectSlotsByName(Parameters)

    def retranslateUi(self, Parameters): pass
        
class Box:
    def __init__(self,Parameters,parent,gr,nb):
        self.box = QtGui.QGroupBox(Parameters)
        self.parent = parent
        y0=20+nb*60
        self.box.setGeometry(QtCore.QRect(10, y0, 170, y0+40))
        self.box.setTitle(QtGui.QApplication.translate("Parameters", gr, None, QtGui.QApplication.UnicodeUTF8))
        self.box.setObjectName(_fromUtf8(gr))
        self.hlWidget = QtGui.QWidget(self.box)
        self.hlWidget.setGeometry(QtCore.QRect(9, 15, 158, 28))
        self.hlWidget.setObjectName(_fromUtf8("horizontalLayoutWidget"))
        self.hl = QtGui.QHBoxLayout(self.hlWidget)
        self.hl.setMargin(0)
        self.hl.setObjectName(_fromUtf8("horizontalLayout"))
        for i in range(len(parent.Base.groups[gr])):
            n=parent.Base.groups[gr][i]
            shortName = gr[2:4]+'_'+n
            # tries to find buttons in addins
            butA = parent.core.addin.addButton(self,gr,i) # a list of buttons
            if butA !=None:
                for short,name in butA : 
                    buta = QtGui.QPushButton(self.hlWidget)
                    buta.setText(QtGui.QApplication.translate("Parameters", short, None, QtGui.QApplication.UnicodeUTF8))
                    buta.setObjectName(_fromUtf8(name))
                    self.hl.addWidget(buta)
                    buta.clicked.connect(self.onButton)
                    
            but = QtGui.QPushButton(self.hlWidget)
            but.setToolTip(QtGui.QApplication.translate("Parameters", n, None, QtGui.QApplication.UnicodeUTF8))
            but.setText(QtGui.QApplication.translate("Parameters", n, None, QtGui.QApplication.UnicodeUTF8))
            # icon = QtGui.QIcon()
            #icon.addPixmap(QtGui.QPixmap(_fromUtf8("F:/iPHT3D/Lib2_b/utils/Ch_P.gif")), QtGui.QIcon.Normal, QtGui.QIcon.Off)
            #but.setIcon(icon)
            but.setObjectName(_fromUtf8(n))
            but.clicked.connect(self.onButton)
            self.hl.addWidget(but)
            
    def onButton(self):
        s = self.parent.gui.sender()
        name = s.objectName()
        self.parent.Base.action(name)
