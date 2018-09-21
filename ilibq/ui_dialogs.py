# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'ui_file.ui'
#
# Created: Sat Feb 15 15:09:21 2014
#      by: PyQt4 UI code generator 4.8.6
#
# WARNING! All changes made in this file will be lost!

from PyQt4.QtCore import *
from PyQt4.QtGui import *

#try:
#    _fromUtf8 = QString.fromUtf8
#except AttributeError:
#    _fromUtf8 = lambda s: s
from guiShow import *
from parameters import BaseParms

class Ui_Main(object):
    def setupUi(self,Main,gui,core):
        Main.setObjectName("Main")
        Main.resize(200,560)
        Main.setWindowTitle("Main")
        self.layoutWidget = QWidget(Main)
        self.layoutWidget.setGeometry(QRect(5, 5, 195, 480)) #left,top,w,h
        self.toolBox = QToolBox(self.layoutWidget)
        self.toolBox.setGeometry(QRect(0, 0, 195, 480))
        self.page = QWidget()
        self.pFile = Ui_File()
        self.pFile.setupUi(self.page,gui)
        self.toolBox.addItem(self.page,"Files and tools")
        
        self.page1 = QWidget()
        self.pParameters = Ui_Parameters()
        self.pParameters.setupUi(self.page1,gui,core)
        self.toolBox.addItem(self.page1,"Model parameters")
        
        self.page2 = QWidget()
        self.pVar = Ui_Var()
        self.pVar.setupUi(self.page2,gui,core)
        self.toolBox.addItem(self.page2,"Spatial variables")
        
        self.page3 = QWidget()
        self.pShow = guiShow(gui,core)
        self.pShow.dialg.setupUi(self.page3,self.pShow,gui,core)
        self.toolBox.addItem(self.page3,"View results")

        QMetaObject.connectSlotsByName(Main)

class Ui_File(object):
    def setupUi(self, File,gui):
#        File.setObjectName("File")
        File.resize(177, 104)
        self.gui = gui
        self.gridLayoutWidget = QWidget(File)
        self.gridLayoutWidget.setGeometry(QRect(0, 0, 180, 120))
        self.gridLayout = QGridLayout(self.gridLayoutWidget)
        self.gridLayout.setMargin(0)
        self.comboFile = QComboBox(self.gridLayoutWidget)
        self.comboFile.addItems(["","Open","Save", "Save as"])
        self.gridLayout.addWidget(self.comboFile, 0, 1, 1, 1)
        self.comboImport = QComboBox(self.gridLayoutWidget)
        self.comboImport.addItems(["","test"])
        self.gridLayout.addWidget(self.comboImport, 1, 1, 1, 1)
        self.label = QLabel(self.gridLayoutWidget)
        self.label.setText("Import")
        self.gridLayout.addWidget(self.label, 1, 0, 1, 1)
        self.label_2 = QLabel(self.gridLayoutWidget)
        self.label_2.setText("File")
        self.gridLayout.addWidget(self.label_2, 0, 0, 1, 1)

        self.comboImport.activated['QString'].connect(self.onClick)
        self.retranslateUi(File)
        QMetaObject.connectSlotsByName(File)

    def retranslateUi(self, File):
        pass
    def onClick(self):
        self.gui.onMessage('hello')

class Ui_Parameters(object):
    def setupUi(self, Parameters,gui,core):
        self.gui,self.core =gui, core
        self.Base = BaseParms(gui,core)
        Parameters.setObjectName("Parameters")
        Parameters.resize(197, 348)
        Parameters.setWindowTitle( "Parameters")
        self.dictBox={}
        skey = self.Base.groups.keys(); skey.sort()
        for i,g in enumerate(skey): 
            self.dictBox[g] = Box(Parameters,self,g,i)
        QMetaObject.connectSlotsByName(Parameters)
        
class Box:
    def __init__(self,Parameters,parent,gr,nb):
        self.box = QGroupBox(Parameters)
        self.parent = parent
        y0=20+nb*60
        self.box.setGeometry(QRect(5, y0, 170, y0+40))
        self.box.setTitle(gr)
        self.hlWidget = QWidget(self.box)
        self.hlWidget.setGeometry(QRect(9, 15, 158, 28))
        self.hl = QHBoxLayout(self.hlWidget)
        self.hl.setMargin(0)
        for i in range(len(parent.Base.groups[gr])):
            n=parent.Base.groups[gr][i]
            shortName = gr[2:4]+'_'+n
            # tries to find buttons in addins
            butA = parent.core.addin.addButton(self,gr,i) # a list of buttons
            if butA !=None:
                for short,name in butA : 
                    buta = QPushButton(self.hlWidget)
                    buta.setText(QApplication.translate("Parameters", short, None, QApplication.UnicodeUTF8))
                    buta.setObjectName(name)
                    self.hl.addWidget(buta)
                    buta.clicked.connect(self.onButton)
                    
            but = QPushButton(self.hlWidget)
            but.setToolTip(QApplication.translate("Parameters", n, None, QApplication.UnicodeUTF8))
            but.setText(QApplication.translate("Parameters", n, None, QApplication.UnicodeUTF8))
            # icon = QIcon()
            #icon.addPixmap(QPixmap(_fromUtf8("F:/iPHT3D/Lib2_b/utils/Ch_P.gif")), QIcon.Normal, QIcon.Off)
            #but.setIcon(icon)
            but.setObjectName(shortName) #_fromUtf8(n))
            but.clicked.connect(self.onButton)
            self.hl.addWidget(but)
            
    def onButton(self):
        s = self.parent.gui.sender()
        name = s.objectName()
        self.parent.Base.action(name)

class Ui_Var(object):
    def setupUi(self, Var,gui,core):
        Var.setObjectName("Var")
        #Var.resize(177, 104)
        Var.setWindowTitle(QApplication.translate("Var", "Var", None, QApplication.UnicodeUTF8))
        self.gridLayoutWidget = QWidget(Var)
        self.gridLayoutWidget.setGeometry(QRect(0, 0, 180, 180))
        #self.gridLayoutWidget.setObjectName(_fromUtf8("gridLayoutWidget"))
        self.gridLayout = QGridLayout(self.gridLayoutWidget)
        self.gridLayout.setMargin(0)
        #self.gridLayout.setObjectName(_fromUtf8("gridLayout"))

        label = QLabel(self.gridLayoutWidget)
        label.setText(QApplication.translate("", "Model", None, QApplication.UnicodeUTF8))
        self.gridLayout.addWidget(label, 0, 0, 1, 1)
        self.comboModel = QComboBox(self.gridLayoutWidget)
        #self.comboFile.setObjectName(_fromUtf8("comboFile"))
        models = core.modelList
        for i,n in enumerate(models):
            self.comboModel.addItem("")
            self.comboModel.setItemText(i, QApplication.translate("File", n, None, QApplication.UnicodeUTF8))
        self.gridLayout.addWidget(self.comboModel, 0, 1, 1, 1)

        label = QLabel(self.gridLayoutWidget)
        label.setText(QApplication.translate("", "Group", None, QApplication.UnicodeUTF8))
        self.gridLayout.addWidget(label, 1, 0, 1, 1)
        self.comboGroup = QComboBox(self.gridLayoutWidget)
        #self.comboImport.setObjectName(_fromUtf8("comboImport"))
        self.gridLayout.addWidget(self.comboGroup, 1, 1, 1, 1)

        label = QLabel(self.gridLayoutWidget)
        label.setText(QApplication.translate("", "Line", None, QApplication.UnicodeUTF8))
        self.gridLayout.addWidget(label, 2, 0, 1, 1)
        self.comboLine = QComboBox(self.gridLayoutWidget)
        #self.comboImport.setObjectName(_fromUtf8("comboImport"))
        self.gridLayout.addWidget(self.comboLine, 2, 1, 1, 1)

        label = QLabel(self.gridLayoutWidget)
        label.setText(QApplication.translate("", "Type", None, QApplication.UnicodeUTF8))
        self.gridLayout.addWidget(label, 3, 0, 1, 1)
        self.comboType = QComboBox(self.gridLayoutWidget)
        typeList = ['one_value','formula','zone','edit','interpolate','import']
        for i,n in enumerate(typeList):
            self.comboType.addItem("")
            self.comboType.setItemText(i, QApplication.translate("Type", n, None, QApplication.UnicodeUTF8))
        self.gridLayout.addWidget(self.comboType, 3, 1, 1, 1)

        label = QLabel(self.gridLayoutWidget)
        label.setText(QApplication.translate("", "Value", None, QApplication.UnicodeUTF8))
        self.gridLayout.addWidget(label, 4, 0, 1, 1)
        self.value = QLineEdit(self.gridLayoutWidget)
        self.gridLayout.addWidget(self.value, 4, 1, 1, 1)
        
        self.retranslateUi(Var)
        QMetaObject.connectSlotsByName(Var)

    def retranslateUi(self, Var):
        pass

