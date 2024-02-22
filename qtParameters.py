# -*- coding: utf-8 -*-


from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *

import os
from parameters import BaseParms
from qtDialogs import *
from core import *
from config import *

from defSPH import *
from configSPH import *

class Ui_Parameters(object):
    def setupUi(self, Parameters,gui,core,plugin_dir):
        self.gui,self.core,self.plugin_dir =gui, core,plugin_dir
        self.base = BaseParms(gui,core)
        Parameters.setObjectName("Parameters")
        self.mainbx = QVBoxLayout(Parameters) # this and 5 lines below added OA 6/11
        title = QLabel(Parameters)
        title.setText("Parameters")
        font = QFont();font.setPointSize(10);font.setBold(True)
        title.setFont(font)
        title.setMaximumHeight(30)
        self.mainbx.addWidget(title)
        self.dictBox={}
        skey = list(self.base.groups.keys()); skey.sort()
        for i,g in enumerate(skey): 
            self.dictBox[g] = Box(Parameters,self,g,i)
        self.mainbx.addStretch(0) # OA 6/11

        ############################################################################################################
        #modif SPH 2023: add an history button
        #               if modifSPH2023_cleanSaveDirAtStart is true, also remove the previous content of the save dir
        #############################################################################################################
        if modifSPH2023:
            self.mode = modifSPH2023_backupFormat
            self.followHistory = \
                (modifSPH2023_replayScheme == sph_possibleReplayScheme.checkedRestoreFinalValues_uncheckedReplayHistory)
            #self.history = []
            SPHbox_Q = QPushButton(Parameters)
            SPHbox_Q.setText("history")
            font = QFont()
            font.setPointSize(10)
            font.setBold(False)
            SPHbox_Q.setFont(font)
            SPHbox_Q.setMaximumHeight(30)
            self.mainbx.addWidget(SPHbox_Q)
            SPHbox_Q.clicked.connect(self.sph_interpreteSaveButtonClicked)

            sph_manageBackupDirectory()

        ####################################
        #end modif SPH 2023
        ##################################

        QMetaObject.connectSlotsByName(Parameters)

    ###################################
    #modif SPH 2023 & 2024: add an history button
    ###################################
    """
    personal notations (SPH):
    1. My QWidgets all have a suffix '_Q', conversely if a parameter name ends with '_Q' then it is a QWidgets
    2. The methods connected to a QWidget generally have a name beginning with 'sph_interprete'
    3. Generally I am in lower camel case,
        the symbol '_' being only the separator with the 'sph' prefix, the 'Q' suffix or the 'modifSPH2023' prefix
    4. My method names generally begins with 'sph_'
    """
    def sph_interpreteSaveButtonClicked(self,backupList=[]):
        """
        to manage all the history interface developed by SPH:
        - a table that shows all backups, each backup can be checked
        - a replay button to reload all backups that are checked
        - a set of buttons that seem convenient to me

        Parameters
        ----------
            backupList : the string 'dev' if debugging
                this 'dev' option is now useless and this parameter 'backupList' is useless

        Returns
        -------
            nothing
        """
        backupList = sph_constructBackupList(backupList=backupList)
        #sph_makeHistoryInterface(backupList=backupList)

        if modifSPH2023 and modifSPH2023_usePickle:
            import pickle as pk
        if modifSPH2023_debugLevel>5:
            print('entering saveButtonClicked...')
            for count,backup in enumerate(backupList):
                print('backup num',str(count),': ',backup)


        # note SPH 2023: adaptation from https://stackoverflow.com/questions/46830464/applying-a-layout-to-a-qdockwidget-in-pyqt5
        # note 'self.gui' works and not 'self.mainbx' because the latter is a QVBoxLayout,
        # which can not be parent of a QDockWidget !
        docked_Q = QtWidgets.QDockWidget('history list',self.gui)
        self.gui.addDockWidget(QtCore.Qt.RightDockWidgetArea, docked_Q)# is it useful for such a 'floating' window ?
        dockedWidget_Q = QtWidgets.QWidget(self.gui)# is 'self.gui' useful here ? TBV
        docked_Q.setWidget(dockedWidget_Q)
        #test 240212
        #dockedWidget_Q.setLayout(QtWidgets.QVBoxLayout())
        mainLayout = QtWidgets.QVBoxLayout()
        dockedWidget_Q.setLayout(mainLayout)

        if modifSPH2023_debugLevel > 10:
            print('debug 240209:',sph_possibleBackupFormat,sph_possibleBackupFormat.memory.value)


        self.backupList = backupList
        ####################################
# note 240217: inverse the numbering in qtableview seems hard,
        # rather it seems to me easier to mask this automatic numbering and to make a good one

        #building the table
        #240217: add number col
        col = ['chrono','check2play', 'time', 'date', 'model group', 'line','hover4details']
        self.posCheck2Play = col.index('check2play')

        self.backupTable_Q = QtWidgets.QTableWidget(dockedWidget_Q)
        self.backupTable_Q.setRowCount(len(backupList))
        self.backupTable_Q.setColumnCount(len(col))
        self.backupTable_Q.setHorizontalHeaderLabels(col)

        # 240217: add personal number col
        #self.backupTable_Q.horizontalHeader().hide()

        if modifSPH2023_hideVerticalHeaderInBackupTable:
            self.backupTable_Q.verticalHeader().hide()



        for count,backup in enumerate(backupList):
            if modifSPH2023_debugLevel > 8:
                print('adding',backup.fileName,'...')
            #note SPH 2023: in the file names I used camel case, so I regive the true name here
            #that can be, at the moment, 'Modflow series' (noted 'ModflowSeries' in the file names), 'Modflow USG' (noted ModflowUsg) or 'Min3p'

            # test 240209
            #print('debug 240209c:',backup)

            self.backupTable_Q.setCellWidget(count, 4, QtWidgets.QLabel(
                ' '+backup.modelGroup.replace('Usg',' USG').replace('Series',' series')+' '))
            self.backupTable_Q.setCellWidget(count, 5, QtWidgets.QLabel(' '+backup.currentLine+' '))
            self.backupTable_Q.setCellWidget(count, 3, QtWidgets.QLabel(' '+backup.date+' '))
            self.backupTable_Q.setCellWidget(count, 2, QtWidgets.QLabel(' '+backup.time+' '))
            temp = QtWidgets.QCheckBox()
            temp.setChecked(True)
            self.backupTable_Q.setCellWidget(count,1, temp)
            self.backupTable_Q.setCellWidget(count,0, QtWidgets.QLabel('  '+str(len(backupList)-count)))

# test240214 for tooltip
            temp = QtWidgets.QLabel(' >> details here << ')
            #for param in backup.init
            toolTipString = ""
            currentSave = sph_loadSave(backupList=modifSPH2023_saveList,
                                       id=count,
                                       fileName=backupList[count].fileName,
                                       mode=modifSPH2023_backupFormat)
            isChanged = False
            for count2,param in enumerate(currentSave['initialValues']):
                if param != currentSave['finalValues'][count2]:
                    isChanged = True
                    toolTipString += 'par. n°'+str(count2+1)+': '+str(param)+'-->'+str(currentSave['finalValues'][count2])+'\n'
            if not isChanged:
                toolTipString += 'no change'
            temp.setToolTip(toolTipString)

            self.backupTable_Q.setCellWidget(count, 6, temp)

        horizontalHeader = self.backupTable_Q.horizontalHeader()
        for count,temp in enumerate(col):
            horizontalHeader.setSectionResizeMode(count, QtWidgets.QHeaderView.ResizeToContents)
        dockedWidget_Q.layout().addWidget(self.backupTable_Q)
        if len(backupList) == 0:
            temp = QtWidgets.QLabel('no backup yet !')
            temp.setFont(QFont('Times', 16))
            dockedWidget_Q.layout().addWidget(temp)
######################

        #window.backupList_Q.itemDoubleClicked.connect(window.interpreteBackupDoubleClick)
        #window.backupList_Q.itemClicked.connect(window.interpreteBackupClick)
        #self.backupTable_Q.itemClicked.connect(self.sph_interpreteBackupClicked)# does not work ?! TBV !
        #upButton_Q = QtWidgets.QPushButton()
        #upButton_Q.clicked.connect(window.interpreteUpButtonClick)
        #downButton_Q = QtWidgets.QPushButton()
        #upButton_Q.clicked.connect(window.interpreteAddButtonClick)
        #addButton_Q = QtWidgets.QPushButton()

# test 240212
        buttonLayout = QtWidgets.QHBoxLayout()
        leftButtons_L = QtWidgets.QVBoxLayout()
        rightButtons_L = QtWidgets.QVBoxLayout()
        mainLayout.addLayout(buttonLayout)
        buttonLayout.addLayout(leftButtons_L)
        buttonLayout.addLayout(rightButtons_L)
        #dockedWidget_Q.layout().setlayout(QtWidgets.QHBoxLayout())


        replayButton_Q = QtWidgets.QPushButton()
        replayButton_Q.setText('replay')
        replayButton_Q.setFont(QFont('Times', 16))
        # setFonttSize(16)#setStylesheet("font-size: 16px")

        # test 240212
        #dockedWidget_Q.layout().addWidget(replayButton_Q)
        leftButtons_L.addWidget(replayButton_Q)

        #just for debugging
        if modifSPH2023_debugLevel>10:
            self.replayMode = False
            replayButton_Q.clicked.connect(self.sph_interpreteLoad)
            print('WARNING: this is a fast debug where replay is equiv to a simple load')
        else:
            replayButton_Q.clicked.connect(self.sph_interpreteReplay)

        checkAllButton_Q = QtWidgets.QPushButton()
        checkUpToHereButton_Q = QtWidgets.QPushButton()
        uncheckAllButton_Q = QtWidgets.QPushButton()
        refreshButton_Q = QtWidgets.QPushButton()

        # test 240212 canceled
        removeButton_Q = QtWidgets.QPushButton()
        # removeButton_Q = QtWidgets.QPushButton(dockedWidget_Q)

        removeButton_Q.setText('remove')





        checkAllButton_Q.setText('check all')
        uncheckAllButton_Q.setText('uncheck all')
        checkUpToHereButton_Q.setText('check up to here')
        refreshButton_Q.setText('refresh')

        # test 240212" \
        #dockedWidget_Q.layout().addWidget(checkAllButton_Q)
        #dockedWidget_Q.layout().addWidget(uncheckAllButton_Q)
        #dockedWidget_Q.layout().addWidget(checkUpToHereButton_Q)
        rightButtons_L.addWidget(refreshButton_Q)
        rightButtons_L.addWidget(checkAllButton_Q)
        rightButtons_L.addWidget(uncheckAllButton_Q)
        rightButtons_L.addWidget(checkUpToHereButton_Q)
        # test 240212
        # dockedWidget_Q.layout().addWidget(removeButton_Q)
        rightButtons_L.addWidget(removeButton_Q)

        checkAllButton_Q.clicked.connect(self.sph_interpreteCheckAll)
        uncheckAllButton_Q.clicked.connect(self.sph_interpreteUncheckAll)
        checkUpToHereButton_Q.clicked.connect(self.sph_interpreteCheckUpToHere)
        removeButton_Q.clicked.connect(self.sph_interpreteRemoveButton)
        refreshButton_Q.clicked.connect(self.sph_interpreteRefreshButton)

        # window.vBox_Q.addWidget(window.backupList_Q)
        # # vBox_Q.addWidget(addButton_Q)
        # window.vBox_Q.addWidget(window.removeButton_Q)
        # window.vBox_Q.addWidget(window.infoButton_Q)
        #vBox_Q = QtWidgets.QVBoxLayout()
        ##toto = QtWidgets.QWidget
        ##toto.setLayout(vBox_Q)
        ##window.setWidget(toto)
        # window.vBox_Q = QtWidgets.QVBoxLayout(window)
        #window.layout_Q = QtWidgets.QVBoxLayout()
        #window.setLayout(window.layout_Q)
        #window.setWindowTitle('save list')
        docked_Q.setVisible(True)
        docked_Q.setFloating(True)
        docked_Q.setGeometry(400, 300, 650, 650)

# TODO 240202:
#    - when unchecked, back to initial values -> ok
# - remove all saves wnen closing -> rather I remove all when opening...

    def sph_interpreteRefreshButton(self):
        #self.backupList = sph_constructBackupList()
        self.sph_interpreteSaveButtonClicked()

    def sph_interpreteRemoveButton(self):
        """
        action when the 'remove' button is pressed
        note only one remove action can be done at a time, allowing the simple use of 'pop' method
        """
        result = QtWidgets.QMessageBox.question(self.gui,
                                                    "Confirm Removing..",
                                                    "Are you sure you want to definitely remove this save ?",
                                                    QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No)
        if result == QtWidgets.QMessageBox.Yes:
            item = self.backupTable_Q.currentRow()
            # add 240219
            if modifSPH2023_backupFormat == sph_possibleBackupFormat.pickle:
                os.remove(os.path.join(os.getcwd(),'lib','ilibq',modifSPH2023_saveDir,self.backupList[item].fileName))
            self.backupList.pop(item)
            self.backupTable_Q.removeRow(item)
            # add 240209, corrected 240212
            if modifSPH2023_backupFormat == sph_possibleBackupFormat.memory:
                modifSPH2023_saveList.pop(item)
            self.sph_interpreteRefreshButton()

    def sph_interpreteLoad(self):
        """
        action when the 'load' button is pressed
        also used as the 'replay' action for debugging...
        """
        if self.replayMode:
            item = self.internalCount
        else:
            item = self.backupTable_Q.currentRow()
        if modifSPH2023_debugLevel>15:
            print('debug: loading save',item+1)
            print('debug 240110:',hasattr(self,'values'),hasattr(self,'boxkeys'))
        # adapted 240209
        # rk1: in pickle mode 'id' & backupList are useless ; in memory mode 'fileName' is useless
        # rk2: modifSPH2023_saveList is self.backupList (memory ref). TBV
        if modifSPH2023_debugLevel>5:
            print('debug 240212d:',self.backupList)
            for count,test in enumerate(self.backupList):
                print('   ***',count,test.currentLine)
            print('debug 240212e:',modifSPH2023_saveList)
        currentSave = sph_loadSave(backupList=modifSPH2023_saveList,
                                   id=item,
                                   fileName=self.backupList[item].fileName,
                                   mode=modifSPH2023_backupFormat)

        #print('*** debug 240123:', inspect.getmembers(self.core))
# looks like a list encapsulating a tuple with the second element being a dictionary
        #print('*** debug 240123:', type(self.core))#[1]['Modflow'])
        if modifSPH2023_debugLevel>20:
            print('*** debug 240216:', self.core.KwList)  # this looks like a dict of all possible 'currentLine' names
        # with their positions
        #print('*** debug 240123c:', self.core.dicval)# this looks like a dict of all values. that's what we need to overwrite when loading !!!
    #dialg = genericDialog(self.gui,'value',lst0)
        #values = dialg.getValues()
        #print('debug 240122b:',self.core.ttable)
        #print('debug 240122:',self.core.getTlist2())

        #print('debug 240125d:',self.core.dicval,hasattr(self.core,'dicaddin'),(self.core.dicaddin['Model'])['group'])
        if currentSave['backupType'] == sph_backupType.leftColumn:

            #for groupKey in self.core.dicval.keys():
            #    if currentSave['currentLine'] in self.core.dicval[groupKey].keys():
            #        print('debug 240219f, yep:',groupKey)

            for groupKey in self.core.dicval.keys():
                if groupKey.startswith('Opgeo'):#according to OA, this model was given up
                    # but it has few repeated keywords that mess things up
                    continue
                if currentSave['currentLine'] in self.core.dicval[groupKey].keys():

                    if modifSPH2023_debugLevel>5:
                        print('### ### debug 240129h:',currentSave['currentLine'])

                    d = (self.core.dicval[groupKey])#[currentSave['currentLine']]

                    #print('debug 240219f:',self.core.dicval)

                    if modifSPH2023_debugLevel>5:
                        print('debug 240129g:',d)
                #d = self.core.dicval[currentSave['currentModel']]


            if modifSPH2023_debugLevel>5:
                print(currentSave['modelGroup'],'.',currentSave['currentLine'],'is affected with initialValues:',currentSave['initialValues'])
                print('   changed to',currentSave['finalValues'])

            if self.loadingMode == sph_loadingMode.final:
                if self.followHistory:
                    # add 240219: if followHistory (for instance in the current replay mode)
                    # then final values are always set in case of change with the history updated
                    for count,param in enumerate(currentSave['initialValues']):
                        if param != currentSave['finalValues'][count]:
                            (d[currentSave['currentLine']])[count] = currentSave['finalValues'][count]
                            self.history.append( (currentSave['currentLine'],count) )
                else:#no loop on parameters to gain time
                    d[currentSave['currentLine']] = currentSave['finalValues']
                (self.core.dicaddin['Model'])['group'] = currentSave['modelGroup']
            elif self.loadingMode == sph_loadingMode.initial:
                if self.followHistory:
                    # add 240219: if followHistory (for instance in the current replay mode)
                    # then initiaValues are set only if this param is not already in the history (then it is a default value)
                    for count,param in enumerate(currentSave['initialValues']):
                        if (currentSave['currentLine'],count) not in self.history:

                            if modifSPH2023_debugLevel > 20:
                                print('debug 240219e:',count, d,d[currentSave['currentLine']],currentSave['initialValues'])

                            (d[currentSave['currentLine']])[count] = currentSave['initialValues'][count]
                            self.history.append( (currentSave['currentLine'],count) )
                else: #no loop on parameters to gain time
                    d[currentSave['currentLine']] = currentSave['initialValues']
        else:# case introduced 240212
            pass#TODO


    def sph_interpreteReplay(self):
        """
        action of the 'replay' button
        """
        self.replayMode=True
        if modifSPH2023_replayScheme == sph_possibleReplayScheme.checkedRestoreFinalValues_uncheckedReplayHistory:
            self.history = []
        for count,b in enumerate(self.backupList):

            if modifSPH2023_debugLevel>20:
                print('debug 240219b:',self.history)

            self.internalCount = len(self.backupList)-1-count
            temp = self.backupTable_Q.cellWidget(self.internalCount,self.posCheck2Play)
            if modifSPH2023_debugLevel>5:
                print('debug 240124:',self.internalCount,temp,temp.text(),temp.isChecked())
            if temp.isChecked():
                self.loadingMode = sph_loadingMode.final

                if modifSPH2023_debugLevel > 20:
                    print('debug 240219c:', self.history)

            else:
                # modif SPH 240202: loading scheme changed
                # self.loadingMode = sph_loadingMode.noLoading
                if modifSPH2023_replayScheme == sph_possibleReplayScheme.checkedRestoreFinalValues_uncheckedDoNothing:
                    self.loadingMode = sph_loadingMode.noLoading
                else:
                    self.loadingMode = sph_loadingMode.initial

                    if modifSPH2023_debugLevel > 20:
                        print('debug 240219d:', self.history)

            # if self.loadingMode is noLoading, then the call of self.sph_interpreteLoad is useless.
            if self.loadingMode != sph_loadingMode.noLoading:
                self.sph_interpreteLoad()

            if modifSPH2023_debugLevel>20:
                print('debug 240219:',self.history)

    def sph_interpreteCheckAll(self):
        """
        actionof the "check all' button
        """
        for count,b in enumerate(self.backupList):
            self.backupTable_Q.cellWidget(count, self.posCheck2Play).setChecked(True)

    def sph_interpreteUncheckAll(self):
        """
        action of the "uncheck all' button
        """
        for count, b in enumerate(self.backupList):
            self.backupTable_Q.cellWidget(count, self.posCheck2Play).setChecked(False)

    def sph_interpreteCheckUpToHere(self):
        """
        action of the 'check up to here' button
        """
        for count, b in enumerate(self.backupList):
            recount = len(self.backupList)-1-count
            condition = ( recount>=self.backupTable_Q.currentRow() )
            self.backupTable_Q.cellWidget(recount, self.posCheck2Play).setChecked(condition)

    ####################################
    #end modif SPH 2023
    ##################################

def boutonVisible(wdow,nomBut,bool): #EV 05/08/19
    it=wdow.findChild(QPushButton,nomBut)
    if it :
        it.setEnabled(bool)
        
def boutonIcon(wdow,nomBut,file):  # OA 22/8/19
    it=wdow.findChild(QPushButton,nomBut)
    if it :
        it.setIcon(QIcon(wdow.mainDir+os.sep+wdow.u_dir+os.sep+file))
        
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
        parent.mainbx.addWidget(self.box) # OA 6/11/18
        #self.parent.gui.dialogs.onMessage(self.parent.gui,os.listdir(dirutils)[0])
        #policy = QSizePolicy(QSizePolicy.Preferred, QSizePolicy.Preferred)

        butA = parent.core.addin.addButton(self,gr) # a list of buttons
        if butA !=None:
            #print os.listdir(dirutils)
            for short,name,pos in butA : 
                if pos==1 : continue
                buta = QPushButton()
                shortName = name+'.png'#;print ('a', shortName)
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
            shortName = gr[2:4]+'_'+n;#print (i,shortName)
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
