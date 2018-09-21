# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'ui_file.ui'
#
# Created: Sat Feb 15 15:09:21 2014
#      by: PyQt4 UI code generator 4.8.6
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    _fromUtf8 = lambda s: s

class Ui_File(object):
    def setupUi(self, File):
        File.setObjectName(_fromUtf8("File"))
        File.resize(177, 104)
        File.setWindowTitle(QtGui.QApplication.translate("File", "File", None, QtGui.QApplication.UnicodeUTF8))
        self.gridLayoutWidget = QtGui.QWidget(File)
        self.gridLayoutWidget.setGeometry(QtCore.QRect(10, 19, 148, 71))
        self.gridLayoutWidget.setObjectName(_fromUtf8("gridLayoutWidget"))
        self.gridLayout = QtGui.QGridLayout(self.gridLayoutWidget)
        self.gridLayout.setMargin(0)
        self.gridLayout.setObjectName(_fromUtf8("gridLayout"))
        self.comboFile = QtGui.QComboBox(self.gridLayoutWidget)
        self.comboFile.setObjectName(_fromUtf8("comboFile"))
        self.comboFile.addItem(_fromUtf8(""))
        self.comboFile.setItemText(0, QtGui.QApplication.translate("File", "Open", None, QtGui.QApplication.UnicodeUTF8))
        self.comboFile.addItem(_fromUtf8(""))
        self.comboFile.setItemText(1, QtGui.QApplication.translate("File", "Save", None, QtGui.QApplication.UnicodeUTF8))
        self.comboFile.addItem(_fromUtf8(""))
        self.comboFile.setItemText(2, QtGui.QApplication.translate("File", "Save as", None, QtGui.QApplication.UnicodeUTF8))
        self.gridLayout.addWidget(self.comboFile, 0, 1, 1, 1)
        self.comboImport = QtGui.QComboBox(self.gridLayoutWidget)
        self.comboImport.setObjectName(_fromUtf8("comboImport"))
        self.gridLayout.addWidget(self.comboImport, 1, 1, 1, 1)
        self.label = QtGui.QLabel(self.gridLayoutWidget)
        self.label.setText(QtGui.QApplication.translate("File", "Import", None, QtGui.QApplication.UnicodeUTF8))
        self.label.setObjectName(_fromUtf8("label"))
        self.gridLayout.addWidget(self.label, 1, 0, 1, 1)
        self.label_2 = QtGui.QLabel(self.gridLayoutWidget)
        self.label_2.setText(QtGui.QApplication.translate("File", "File", None, QtGui.QApplication.UnicodeUTF8))
        self.label_2.setObjectName(_fromUtf8("label_2"))
        self.gridLayout.addWidget(self.label_2, 0, 0, 1, 1)

        self.retranslateUi(File)
        QtCore.QMetaObject.connectSlotsByName(File)

    def retranslateUi(self, File):
        pass

