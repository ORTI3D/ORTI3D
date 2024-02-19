import sys
from PyQt5 import QtWidgets
from orti3dGui import *

if __name__ == '__main__':

    app = QtWidgets.QApplication(sys.argv)

    mainWindow = orti3dGui('iQpht3d Open Reactive Transport')
    #print(mainWindow,mainWindow.widget)
    mainWindow.setGeometry(100, 100, 800, 500)
    mainWindow.show()

    sys.exit(app.exec_())

