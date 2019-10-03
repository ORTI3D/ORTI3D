# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'ui_file.ui'
#
# Created: Sat Feb 15 15:09:21 2014
#      by: PyQt5 UI code generator 4.8.6

from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
import os
from .qtDialogs import *
from .geometry import *
from .core import *
from .config import *
from .guiShow import guiShow
        
def selectComboValue(wdow,comboName,txt): #OA 6/9/19
    combo = wdow.findChild(QComboBox,comboName)
    if combo :
        indx = combo.findText(txt)
        if indx >= 0: combo.setCurrentIndex(indx)

class Ui_Show(object):
    def setupUi(self,Show,gui,core):
        from .multiPlot import multiPlot
        self.Show = Show
        self.gui,self.core = gui,core
        self.guiShow = guiShow(gui,core)
        self.groups = self.guiShow.groups
        Show.setObjectName("Show")
        self.tWidget = QWidget(Show)
        self.screenShape = QDesktopWidget().screenGeometry()
        self.tWidget.setGeometry(QRect(0, 0, self.screenShape.width()*0.105, 45)) 
        self.dictBox={}
        wd0 = QDesktopWidget().screenGeometry().width()/100 # OA 23/2/19
        if gui.gtyp =='qgis': pos=0  # ifelse added 1/6/19 to use this in qgis too
        else : 
            self.makeTop(self,gui)
            pos=40
        for ig in range(len(self.groups)): 
            for g0 in self.groups: #pour ordonner
                if self.groups[g0][0]==ig: g=g0
            names = self.groups[g][1:]  
            self.dictBox[g] = showBox(Show,self,names,g,pos,ig)
            pos += len(names)*wd0*1.3+wd0*2.5 # OA 23/2/19
        selectComboValue(self.gui,'Model_Plane_L','Z')
        QMetaObject.connectSlotsByName(Show)
        
    def makeTop(self,parent,gui):  # OA 1/6/19 separated form above  
        topHBox = QHBoxLayout(parent.tWidget)
        title = QLabel(parent.Show) # this and 5 lines below added OA 6/11
        title.setText("Results")
        font = QFont();font.setPointSize(10);font.setBold(True)
        title.setFont(font)
        title.setMaximumHeight(40)
        topHBox.addWidget(title)
        self.swiPlane = QPushButton()
        self.swiPlane.setIcon(QIcon(gui.u_dir+os.sep+'Vis_OriZ.png'))
        self.swiPlane.setIconSize(QSize(25, 25))
        self.swiPlane.setMaximumWidth(25)
        self.swiPlane.setFlat(True)
        topHBox.addWidget(self.swiPlane)
        self.icCont = QIcon(gui.u_dir+os.sep+'Vis_SwiCont.png')
        self.icImg = QIcon(gui.u_dir+os.sep+'Vis_SwiImg.png')
        self.icImgCont = QIcon(gui.u_dir+os.sep+'Vis_SwiImgCont.png')
        self.swiCont = QPushButton()
        self.swiCont.setIcon(self.icCont)
        self.swiCont.clicked.connect(parent.switchImg) 
        self.swiCont.setIconSize(QSize(25, 25))
        self.swiCont.setMaximumWidth(25)
        self.swiCont.setFlat(True)   
        topHBox.addWidget(self.swiCont)
        topHBox.addStretch(0)
        return topHBox
        
    def switchPlane(self,plane): # OA added 22/8/19
        '''plane is  aletter in x,Y,Z'''
        self.swiPlane.setIcon(QIcon(self.gui.u_dir+os.sep+'Vis_Ori'+plane+'.png'))
        
    def switchImg(self,evt):
        if self.guiShow.swiImg=='Contour':
            self.swiCont.setIcon(self.icImg)
            self.guiShow.swiImg='Image'
        elif self.guiShow.swiImg=='Image':
            self.swiCont.setIcon(self.icImgCont)
            self.guiShow.swiImg='ImageContour'
        elif self.guiShow.swiImg=='ImageContour':
            self.swiCont.setIcon(self.icCont)
            self.guiShow.swiImg='Contour'
        self.guiShow.redraw() # OA 22/8/19

    def getCurrentTime(self):
        combo = self.Show.findChild(QComboBox,'Model_Tstep_L')
        return combo.currentText()
        
    def getNames(self,nameBox):
        combo = self.Show.findChild(QComboBox,nameBox) 
        names = [combo.itemText(i) for i in range(combo.count())]
        return names

    def setNames(self,nameBox,names,opt='strings'):
        combo = self.Show.findChild(QComboBox,nameBox) 
        combo.clear() # danger if set it is not possible to add items
        combo.addItems([str(n) for n in names])

    def uncheckContours(self,group,name,bool):
        """used to uncheck the other contours when group is changed"""
        dic = self.guiShow.dicVisu
        for n,m in [('Flow','Head'),('Flow','Wcontent'),('Transport','Tracer')]:
            self.onTickBox(n,m,'B',False) #modif OA 9/6/19
        if group in ['Flow','Transport']:self.onTickBox(group,name,'B',bool) #added OA 9/6/19

    def onTickBox(self,group,name,tag,bool):
        """ to change the state of a button whithout doing any action"""
        if name == None: return  # added OA 9/6/19
        item = self.Show.findChild(QCheckBox,group+'_'+name+'_'+tag);
        if tag=='B': item.setChecked(bool)
            
    def onClick(self,value):
        """action when a box is clicked, tag L : list """
        item = self.Show.sender()
        n = str(item.objectName()); 
        [group,name,tag]=n.split('_');#QgsMessageLog.logMessage('guish onclick'+str(group)+str(name)+str(tag))
        if tag=='L': 
            if name in ['Layer','Tstep']: 
                retour = item.currentIndex()
            else :
                retour = str(item.currentText()) # case of list, retour is the name
        else: retour = item.isChecked() # a check box retour is True or False 
        if name in self.guiShow.Vtypes['Array']: self.guiShow.resetDicContour()
        #print(retour,self.guiShow.dicVisu)
        self.guiShow.dicVisu[group][name] = retour
        if name == 'Plane': 
            a = shape(self.core.Zblock)
            if len(a) == 3 : nz,ny,nx = a
            else : nz,nn = a;nz +=1; nx,ny = 1,1 # for ogs these are nod layers
            exec('self.guiShow.'+name+'=\"'+str(retour)+'\"')
            self.switchPlane(retour) # OA 22/8/19
            if retour =='Z' : self.setNames('Model_Layer_L',list(range(nz-1)))
            if retour =='Y' : self.setNames('Model_Layer_L',list(range(ny-1)))
            if retour =='X' : self.setNames('Model_Layer_L',list(range(nx-1)))
            self.guiShow.dicVisu['Model']['Layer']=0
        self.guiShow.onClick2(group,name,retour)
            
    def OnChange(self):
        """ change caracteristics of an object"""
        obj = self.Show.sender()
        n = str(obj.objectName());
        [group,name,tag]=n.split('_')
        change = self.guiShow.change
        #item2=self.FindWindowByName(group+'_'+name+'_L');
        #if item2 != None: name=item2.GetStringSelection()
        color = self.guiShow.getGlist(group,name)['color']
        value = self.guiShow.getGlist(group,name)['value']
        if name in list(change.keys()): # case different than contours
            if color == None : color = (0,0,0)
            lst0=[(name,'Color',color)]
            if change[name] != None:
                lst0.append((change[name][0],'Text',change[name][1]))
            dialg = genericDialog(self.gui,name,lst0)
            lst1 = dialg.getValues()
            if lst1 != None:
                color = lst1[0]
                if len(lst1)>1: value=lst1[1]
            else : return
        else: # cas contour
            dialg = dialogContour(self.gui, "Contours",value,color)
            value = dialg.GetStrings()
#            if value != None:
#                c = dlgContour.coul;
#                color=[(c[0].Red(),c[0].Green(),c[0].Blue()),(c[1].Red(),c[1].Green(),c[1].Blue()),
#                     (c[2].Red(),c[2].Green(),c[2].Blue()),int(c[3])];#print 'in change',color
#            else : return
        if name == 'Species': Cgroup,Cname,name = self.guiShow.getCurrentContour() # OA 1/10/19
        self.guiShow.setGlistParm(group,name,'value',value)
        self.guiShow.setGlistParm(group,name,'color',color)
        #self.onTickBox(group,name,tag,True)
        self.gui.visu.changeObject(group,name,value,color)

    def onPlot(self): #EV 19/12/18 
        item = self.Show.findChild(QComboBox,'Observation_Type_L')
        typ = item.currentText()[0]
        if typ=='T' : typ = 'B' 
        elif typ=='C': typ ='X'
        item2 = self.Show.findChild(QComboBox,'Observation_Result_L')
        res = item2.currentText()
        if self.core.diczone['Observation'].dic: #EV 14/08/19
            m = multiPlot(self.gui,self.core,typ,res)
            m.show()
        else : #EV 14/08/19
            mess=onMessage(self.gui,'Create an observation zone to plot the resuts.')
            return mess
        
    # def onObservation(self,group,tstep):
    #     if group not in ['Flow','Transport','Chemistry']: 
    #         onMessage(self.gui,'choose one variable')
    #         return
    #     item = self.Show.findChild(QComboBox,'Observation_Type_L')
    #     typ = item.currentText()[0]  # B P or X
    #     if group=='Chemistry': 
    #         lesp=self.getNames('Chemistry_Species_L')
    #         lesp.extend(self.getNames('Chemistry_User_L'))
    #     elif group=='Flow': lesp=['Head','Flux','Wcontent']
    #     else : lesp=['Transport']
    #     data = list(zip(lesp,['Check']*len(lesp),[False]*len(lesp)))
    #     #dialog to choose species to plot
    #     if len(lesp)>1: 
    #         dlg = genericDialog(self.gui,'species',data)
    #         lst1=dlg.getValues();#onMessage(self.gui,str(lst1))
    #         if lst1 != None:
    #             lst2=[]
    #             for i in range(len(lst1)): 
    #                 if lst1[i] : lst2.append(str(lesp[i]))
    #             lesp = lst2
    #         else :return
    #     #dialog to choose layers if in 3D
    #     layers = 0
    #     if self.core.addin.getDim()=='3D':
    #         data = [('Layers (1, 3-6)','Text','')]
    #         dlg = genericDialog(self.gui,'Select',data)
    #         d = dlg.getValues() #dialog to choose type of graph
    #         if d != None:
    #             layers=d[0]
    #         else :return
    #     # dialog for type of graph, for flow this dialog is useless
    #     if typ != 'X':
    #         if group in ['Chemistry','Transport']: lst0 = ['Value','Weighted value','Mass Discharge (Md)','Mass Flux (J)']
    #         else : lst0 = ['Value','Weighted value','Flow (Q)','Darcy Flux (q)']
    #         data=[('Type','Choice',('Value',lst0))]
    #         dlg = genericDialog(self.gui,'Select',data)
    #         d = dlg.getValues() #dialog to choose details
    #         if d != None:
    #             val=d[0]
    #             typ+=str(lst0.index(val));
    #         else : return
    #     else :
    #         typ+='0' # for Xy graphs it is always the value
    #         
    #     item = self.Show.findChild(QComboBox,'Observation_Zone_L')
    #     znam = item.currentText()
    #     #print('plot',typ,tstep,group,znam,lesp,layers)
    #     dist,val,labl = self.core.onPtObs(typ,tstep,group,znam,lesp,layers);#print 'guishow 263',val
    #     #plt = plot(self.gui,-1)
    #     typ1='-'
    #     if typ[0]=='X': typ1 = '+';
    #     plt = plotxy(self.gui,dist,val,labl[1:],znam,labl[0],"val",typ1);
    #     plt.getValues()
    #     if typ[0]=='X' and lesp[0]=='Head': self.dicplots['X_head']= plt# to be able to recall the dialog
    #     if typ[0]=='X' and lesp[0]=='Transport': self.dicplots['X_tracer']= plt
    #     #plt.Raise()

class showBox:
    def __init__(self,Show,parent,names,gr,pos,ig):
        self.Show,self.parent = Show,parent
        self.group = QGroupBox(Show)
        self.group.setTitle(str(ig+1)+'.'+gr)
        ln = len(names)
        wd0 = QDesktopWidget().screenGeometry().width()/100 # OA 23/2/19
        self.group.setGeometry(QRect(10, pos, wd0*10.5, wd0*2.4+ln*wd0*1.3)) # OA 23/2/19 to see it on laptop
        self.hlWidget = QWidget(self.group)
        self.hlWidget.setGeometry(QRect(3,15 ,wd0*10, wd0*1.4+ln*wd0*1.3)) # OA 23/2/19 to see it on laptop
        boxGrid = QGridLayout(self.hlWidget)
        boxGrid.alignment()
        boxGrid.setContentsMargins(0,0,0,0)
        boxGrid.setSpacing(0)
        self.buts = list(range(len(names)))
        policy = QSizePolicy(QSizePolicy.Preferred, QSizePolicy.Preferred)
        for i,n in enumerate(names):
            if type(n)==type([1,2]): #cas liste -> choix
                name = gr+'_'+n[0]+'_L';
                text = QLabel(self.hlWidget)
                text.setText(n[0])
                boxGrid.addWidget(text,i,0,1,1)
                liste = n[1]
                self.buts[i] = QComboBox(self.hlWidget)
                self.buts[i].setObjectName(gr+'_'+n[0]+'_L')
                self.buts[i].addItems(liste)
                boxGrid.addWidget(self.buts[i],i,1,1,1)
                self.buts[i].activated['QString'].connect(parent.onClick)        
            else : # cas simple : checkbox
                name = gr+'_'+n+'_B';
                text = QLabel(self.hlWidget)
                text.setText(n)
                boxGrid.addWidget(text,i,0,1,1)
                self.buts[i] = QCheckBox(self.hlWidget)
                self.buts[i].setObjectName(gr+'_'+n+'_B')
                boxGrid.addWidget(self.buts[i],i,1,1,1)
                self.buts[i].clicked.connect(parent.onClick)
            self.buts[i].setSizePolicy(policy)
            self.buts[i].setMaximumHeight(18)
            if gr not in ['Model','Observation'] and parent.gui.gtyp != 'qgis': # OA 1/6/19 no C for qgis
                but = QPushButton('C',self.hlWidget)
                but.setObjectName(name[:-2]+'_C')
                but.setSizePolicy(policy)
                but.setMaximumWidth(16)
                but.setMaximumHeight(18)
                but.clicked.connect(parent.OnChange)
                boxGrid.addWidget(but, i,2,1,1)
            #if n in ['Map','Visible','Type','Zone']:but.Enable(False)
            #if n[0] in ['Plane','Layer','Tstep','Units']:but.Enable(False)
            #parent.numCtrl += 2

    def setNames(self,nameBox,names,opt='strings'):
        self.parent.setNames(nameBox,names)
