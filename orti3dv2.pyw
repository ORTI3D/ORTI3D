import sys
from PyQt5 import QtGui
from PyQt5 import QtWidgets
from PyQt5.QtWidgets import QApplication
from orti3dGui import *

def except_hook(cls, exception, traceback):
    sys.__excepthook__(cls, exception, traceback)

if __name__ == '__main__':

    
    app = QtGui.QApplication(sys.argv)
    font=app.font()
    font.setPointSize(8)
    app.setFont(font)
    mainWindow = orti3dGui('Open Reactive Transport interface')
    #print(mainWindow,mainWindow.widget)
    screenShape = QtWidgets.QDesktopWidget().screenGeometry()

    mainWindow.resize(screenShape.width()*.8, screenShape.height()*.8)
    #mainWindow.setGeometry(100, 100, 800, 500)

    mainWindow.show()

    sys.exit(app.exec_())
