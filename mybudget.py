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
import csv

class myBudget(QDialog):
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
        self.setGeometry(QRect(5, 5, screenShape.width()*.75, 
                               screenShape.height()*.7))
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
        self.plgroup.addItems(['Percent Discrepency','In-Out','Time Series',
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
            self.lzname=self.core.diczone['Observation'].dic['obs.1']['name']
            self.label_3 = QtWidgets.QLabel(self.frame3)
            self.label_3.setText("Zone budget zone")
            self.gl3.addWidget(self.label_3,0,0,1,1)
            self.zgroup = QComboBox(self.frame3)
            self.zgroup.addItems([str(n) for n in self.lzname])
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
        self.verticalLayout.addWidget(self.pushButton, 
                                      alignment=Qt.AlignHCenter)
        self.pushButton.clicked.connect(self.buildGraph)
     ## add vertical layout   
        self.horizontalLayout.addLayout(self.verticalLayout)
     ## the right panel vertical layout
        self.verticalLayout2 = QVBoxLayout()
     ## the matplotlib figure 
        self.figure = Figure(tight_layout=True,figsize=(7.8, 3), dpi=100)  
        self.cnv = FigureCanvas(self.figure) 
        #self._ax = self.cnv.figure.subplots()#.add_axes([0.1, 0.15, 0.7, 0.8])
    ## add matplotlib figure
        #self.horizontalLayout.addWidget(self.cnv)
        self.verticalLayout2.addWidget(self.cnv)
        self.toolbar = NavigationToolbar(self.cnv, self) # EV 04/02/20 
        self.verticalLayout2.addWidget(self.toolbar, alignment=Qt.AlignHCenter) 
    ## Export button
        self.pushButton2 = QPushButton(self)
        self.pushButton2.setText('Export')
        self.verticalLayout2.addWidget(self.pushButton2, 
                                       alignment=Qt.AlignHCenter)
        #self.pushButton2.clicked.connect(self.onExport)
    ## add vertical layout 2
        self.horizontalLayout.addLayout(self.verticalLayout2)
        QMetaObject.connectSlotsByName(self)  #OA 1/6/19

############################ Initialize the dialog ############################
        
    def onTstep(self):
    ## show or hide time combobox for the time step graph
        if self.plgroup.currentIndex()==3:
            self.frame2.show()
        else : self.frame2.hide()

    def getObsZone(self):
    ## get the names of model observation zone
        zname=self.lzname#self.core.diczone['Observation'].dic['obs.1']['name']
        zone = self.zgroup.currentText()
        dicZin={'Zones':{}};dicZout={'Zones':{}}
        zIn=[zone+' from '+zname[i] for i in range(len(zname)) 
            if zname[i] != zone]
        zIn.append('Total IN')
        dicZin['Zones'] = list(zip(zIn,[False]*len(zIn)))
        zOut=[zone+' to '+zname[i] for i in range(len(zname))
            if zname[i] != zone]
        zOut.append('Total OUT')
        dicZout['Zones'] = list(zip(zOut,[False]*len(zOut)))
        return dicZin, dicZout
    
    def getBoundaries(self):
    ## get the boundaries in the model (well, drn, storage, chd...)
        dic={'Bound':[]} ; nbound=['STORAGE','CONSTANT HEAD']
        boundTyp=['WELLS','DRAINS','RIVER LEAKAGE','ET','HEAD DEP BOUNDS',
                  'RECHARGE']
        pack=['wel','drn','riv','evt','ghb','rch']
        lmod=self.core.getUsedModulesList('Modflow')
        for i, n in enumerate (pack) :
            if n.upper() in lmod:
                if self.core.diczone['Modflow'].getNbZones(n+'.1')>0:
                    nbound.append(boundTyp[i])
        dic['Bound'] = list(zip(nbound,[False]*len(nbound)))
        return dic
    
    def getMt3dBound(self):
        dic={};
        nbound=['TOTAL IN','TOTAL OUT','SOURCES IN','SINKS OUT',
           'TOTAL MASS IN AQUIFER']
        dic['Bound'] = list(zip(nbound,[False]*len(nbound)))
        return dic
        
    def getSpecies(self):
    ## get the names of model chemical species
        dic={'Species':[]} 
        species = self.core.addin.chem.getListSpecies() 
        dic['Species']=list(zip(species,[False]*len(species)))
        return dic
    
    def getChoices(self):
    ## return a dic in function of type of graph and result to plot
        dicIn={'In':[]}  ; dicOut={'Out':[]} ; dicIO={'In / Out':[]}
        if self.typ == 'M': 
            if self.res!='Transport':
                self.dicBound = self.getBoundaries()
                dicIn['In']=self.dicBound['Bound']+[('TOTAL IN',False)]
                dicOut['Out']=self.dicBound['Bound']+[('TOTAL OUT',False)]
            else:
                self.dicBound = self.getMt3dBound()
                dicIO['In / Out']=self.dicBound['Bound']
                dic = {**dicIO,**dicIn,**dicOut}
                return dic
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
        self.nb.hide()
        dicZin,dicZout=self.getObsZone()
        #print('dicZin',dicZin,'dicZout',dicZout)
        dicIn={'In':[]}  ; dicOut={'Out':[]} 
        dicIn['In']=self.dicBound['Bound']+dicZin['Zones']
        dicOut['Out']=self.dicBound['Bound']+dicZout['Zones']
        if self.res =='Chemistry': 
            dicSpecies=self.getSpecies()
            dic = {**dicSpecies,**dicIn,**dicOut} 
        else : dic = {**dicIn,**dicOut} 
        self.nb = myNoteBookCheck(self.gui,"Options",dic)
        self.hlayout.addWidget(self.nb)
        self.nb.layout.removeWidget(self.nb.buttonBox) 
        self.nb.buttonBox.deleteLater()
        del self.nb.buttonBox
        #self.verticalLayout.addWidget(self.nb)
        self.nb.apply()
        self.nb.show()

###################### Get the option and plot the graphs #####################
    
    def getValues(self):
        ''' Get the value of the noteBookCheck '''
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
        '''Get the plot options from the dialog and put in in dicIn'''
        dicIn={'ptyp':[],'graph':[],'zone':[],'inlist':[],'outlist':[],
               'iolist':[],'splist':[]} # 'lylist':{}
        dicIn['ptyp']=self.typ
        dicIn['graph']=str(self.plgroup.currentText())
        if self.typ == 'Z' : dicIn['zone'] = self.zgroup.currentText()
        dic=self.getValues() 
        print('dic',dic)
        for i in range(len(dic['In'])):
            if dic['In'][i][1]==2:
                if len(dic['In'][i][0].split(' from '))==2:
                    dicIn['inlist'].append(dic['In'][i][0].split(' from ')[1])
                else: dicIn['inlist'].append(dic['In'][i][0])
        for i in range(len(dic['Out'])):
            if dic['Out'][i][1]==2:
                if len(dic['Out'][i][0].split(' to '))==2:
                    dicIn['outlist'].append(dic['Out'][i][0].split(' to ')[1])
                else: dicIn['outlist'].append(dic['Out'][i][0])
        for i in range(len(dic['In / Out'])):
            if dic['In / Out'][i][1]==2:
                dicIn['iolist'].append(dic['In / Out'][i][0])
        #print('dicin',dicIn['inlist']); print('dicout',dicIn['outlist'])
        dicIn['splist']=[self.res]
        if self.res=='Chemistry' :
            dicIn['splist']=[dic['Species'].index(dic['Species'][i]) 
            for i in range(len(dic['Species'])) if dic['Species'][i][1]==2] 
        #print('splist',dicIn['splist'])
        return dicIn
    
    def buildGraph(self):
        dicIn=self.getOptions()
        self.ptyp,self.graph=dicIn['ptyp'],dicIn['graph']
        self.zone,self.splist=dicIn['zone'],dicIn['splist']
        self.inlist,self.outlist= dicIn['inlist'],dicIn['outlist']
        ## Zone budget
        if dicIn['ptyp'] == 'Z':
            self.writeZoneFile(self.inlist,self.outlist,self.zone)
            self.writeInFile()
            self.core.runZonebud()
            self.xy=self.readZBfile(self.graph,self.zone,self.inlist,
                          self.outlist)
            if self.res == 'Flow' : self.plotData(self.xy,self.graph)
        else :
            if self.res == 'Chemistry' : 
                self.xy=self.readPHT3Dfile(self.graph,self.inlist,
                                           self.outlist,self.splist)
            if self. == 'Transport' :
                self.xy=self.readMT3Dfile(self.graph,self.iolist)
            self.plotData(self.xy,self.graph)
    
    def plotData(self,xy,graph):
        '''Plot data as scatter plot or vertical barchart'''
        self.figure.clf()
        self._ax=self.figure.add_subplot(1,1,1)
      ### Scatter plot
        if graph != 'Time Step': 
            if graph != 'Time Series': self._ax.plot(xy['x'],xy['y'])
            else:
                for i, color in enumerate(xy['Cin']):
                    self._ax.plot(xy['x'],xy['yin'][i],c=color,marker='v')
                for i, color in enumerate(xy['Cout']):
                    self._ax.plot(xy['x'],xy['yout'][i],c=color,marker='^')
            self._ax.legend(xy['lab'])
            #self._ax.set_title(self.zolist[i],fontweight="bold", size=9)
            #self._ax.legend(self.llabel,fontsize = 8,loc='best')
            #self._ax.set_ylabel(aylabel, fontsize = 8) 
            #self._ax.set_xlabel(self.axlabel, fontsize = 8)
            self._ax.ticklabel_format(useOffset=False, 
                                      style='sci',scilimits=(-4,4),axis='both'
                                      ,useMathText=True)
            self._ax.tick_params(axis='both', labelsize=8)
       ### Vertical barchart
        else :
            x = np.arange(len(xy['lab']))  # the label locations
            width = 0.35  # the width of the bars
            barIn = self._ax.bar(x - width/2, xy['yin'], width, label='IN')
            barOut = self._ax.bar(x + width/2, xy['yout'], width, label='OUT')
            #self._ax.set_ylabel('')
            #self._ax.set_title('')
            self._ax.set_xticks(x)
            self._ax.set_xticklabels(xy['lab'])
            self._ax.legend()
            self.autolabel(barIn)
            self.autolabel(barOut)
        self._ax.figure.canvas.draw()
                        
    def autolabel(self,rects):
        '''Attach a text label above each bar in *rects*, 
        displaying its height'''
        for rect in rects:
            height = round(rect.get_height(),1) # ;print(height)
            self._ax.annotate('{}'.format(height),
                        xy=(rect.get_x() + rect.get_width() / 2, height),
                        xytext=(0, 3),  # 3 points vertical offset
                        textcoords="offset points",
                        ha='center', va='bottom')
            
################################ Zone budget ##################################        
    
    def writeZoneFile(self,inlist,outlist,zone):
        ''' File with the location of the zones selected by the user. 
        Read by zonbud.exe'''
        nx,ny,xvect,yvect = getXYvects(self.core)
        nmedia = getNmedia(self.core)
        nlayers = getNlayers(self.core)
        lilay = getNlayersPerMedia(self.core)
        self.fullPath = self.core.fileDir+os.sep+self.core.fileName
        f1=open(self.fullPath +'.zone','w')
        f1.write(' %0i %1i %1i ' %(nlayers, ny, nx) +'\n')
        inout=sort(inlist+outlist)
        zname=self.lzname #self.core.diczone['Observation'].dic['obs.1']['name']
        zlist=[zone]
        for i in range(len(zname)):
            if zname[i] in inout:zlist.append(zname[i])
        #print('z',zlist)
        m0 = ones((nlayers,ny,nx)) ; lay=0
        for im in range(nmedia):
            mat=zone2grid(self.core,'Observation','obs.1',im,opt=zlist,iper=0)
            for il in range(int(lilay[im])): ## several layers can exist in each media
                m0[lay]=mat[-1::-1]
                lay +=1
        for mlay in m0:
            f1.write('INTERNAL ('+str(nx)+'I5)\n')
            np.savetxt(f1, mlay, fmt='%4i')
        f1.close()
    
    def writeInFile(self):
        '''File read by zonbud.exe'''
        fDir, fName = self.core.fileDir, self.core.fileName
        f1=open(fDir+'zonbud.in','w')
        s= 'zonbud CSV2 \n'
        s+= fName+'.budget \n'
        s+= 'ZONEBUDGET run \n'
        s+=fName+'.zone \nA'
        f1.write(s);f1.close()
        
    def readZBfile(self,graph,zone,inlist,outlist):
        '''Read zonbud.csv2 file following the graph to plot :
            Percent Discrepency: Percent error vs time
            In-Out: In - out vs time
            Time Series: In Zone or BC and Out Zone or BC vs time
            Time Step: In Zone or BC and Out Zone or BC for 1 time
            Return xy dictionnarie'''
      ### read file and put data in dict td
        td={'title':[],'data':[]}
        df = csv.reader(open(self.core.fileDir+os.sep+'zonbud.2.csv', mode='r')
            ,skipinitialspace=True)
        i=0
        for row in df:
            if i ==0 : 
                row=[w.strip() for w in row]
                td['title']=row
                i+=1
            else : td['data'].append(row)
        
      ### put data in dict xy following type of graph
        xy={'lab':[],'x':[],'y':[]}
        #idz=self.zgroup.currentIndex()
        zname=self.lzname
        idz = zname.index(zone)
        #print('idz+1',str(idz+1))
        
        if graph == 'Percent Discrepency' :
            for row in td['data']:
                if row[td['title'].index('ZONE')]==str(idz+1):
                    xy['x'].append(float(row[td['title'].index('TOTIM')]))
                    xy['y'].append(float(row[td['title'].index('Percent Error')]))
            xy['lab'].append('Percent Discrepency')
            #print ('perc dis',xy)
        
        if graph == 'In-Out':
            for row in td['data']:
                if row[td['title'].index('ZONE')]==str(idz+1):
                    xy['x'].append(float(row[td['title'].index('TOTIM')]))
                    xy['y'].append(float(row[td['title'].index('IN-OUT')]))
            xy['lab'].append('In-Out')
            #print ('in-out',xy)
        
        if graph == 'Time Series':
            ind,lab=self.getIndLab(td,inlist,outlist,zone)
            xy={'lab':[],'x':[],'yin':[],'yout':[],'Cin':[],'Cout':[]}
            i,j,k=0,0,0; lin=[]; lout=[]
            for key, value in ind.items(): 
                yi=[];yo=[] 
                if value[0]:
                    for row in td['data']:
                        if row[td['title'].index('ZONE')]==str(idz+1):
                            if k==0: xy['x'].append(float(
                                    row[td['title'].index('TOTIM')]))
                            yi.append(float(row[int(value[0])]))
                    xy['yin'].append(yi)
                    xy['Cin'].append('C'+str(i))
                    lin.append(lab[key][0])
                    k+=1
                i+=1
                if value[1]:
                    for row in td['data']:
                        if row[td['title'].index('ZONE')]==str(idz+1):
                            if k==0: xy['x'].append(float(
                                    row[td['title'].index('TOTIM')]))
                            yo.append(float(row[int(value[1])]))
                    xy['yout'].append(yo)
                    xy['Cout'].append('C'+str(j))
                    lout.append(lab[key][1])
                    k+=1
                j+=1
            xy['lab']=lin+lout
            #print ('TSeries',xy)
        
        if graph == 'Time Step':
            ind,lab=self.getIndLab(td,inlist,outlist,zone)
            xy={'lab':[],'yin':[],'yout':[]};i=0
            t=self.Tstep.currentText()
            for row in td['data']:
                if row[td['title'].index('ZONE')]==str(idz+1):
                    if float(row[td['title'].index('TOTIM')])==float(t):
                        for key, value in ind.items(): 
                            if value[0]:
                                xy['yin'].append(float(row[int(value[0])]))
                            else : xy['yin'].append(np.nan)
                            if value[1]:
                                xy['yout'].append(float(row[int(value[1])]))
                            else: xy['yout'].append(np.nan)
                            if str(list(lab.keys())[i]).isdigit():
                                xy['lab'].append(zname[list(lab.keys())[i]-1])
                            else:xy['lab'].append(list(lab.keys())[i])
                            i+=1
            #print ('TStep',xy)
        return xy
    
    def getIndLab(self,td,inlist,outlist,zone):
        ''' Return the label and column index for each zone entry and exit 
        and BC. Organized in pairs, 'in' and 'out' index or label are indicated 
        for each zone or BC. Takes the value False if only 'in' or only 'out' 
        is selected for a BC or a zone.
        '''
        inout= list(set(inlist+outlist))
        dicb={}; inT,outT=False,False
        for i in inout :
            if i in ['Total IN','Total OUT']:
                if i in ['Total IN']:inT=True
                if i in ['Total OUT']:outT=True
                dicb['Total']=[inT,outT]
            else :
                if i in inlist :
                    dicb[i]=[True]
                else : dicb[i]=[False]
                if i in outlist :
                    dicb[i].append(True)
                else : dicb[i].append(False)
        ind={};lab={};zname=self.lzname
        for key, value in dicb.items(): 
            i=0
            if key in zname : key = zname.index(key)+1
            lab[key]=[]
            for j,title in enumerate(td['title']) :
                if str(key) in title and i==0 :
                    if value[0]==True: 
                        ind[key]=[j]
                        if str(key)== 'Total': lab[key].append(title)
                        elif len(title.split('FROM '))==1: 
                            lab[key].append(title+' (in)')
                        else : lab[key].append(zone+' from '+str(zname[key-1]))
                    else : 
                        ind[key]=[False]
                        lab[key].append(False)
                    i+=1 ; continue
                if str(key) in title and i==1 :
                    if value[1]==True: 
                        ind[key].append(j)
                        if str(key)== 'Total': lab[key].append(title)
                        elif len(title.split('TO '))==1: 
                            lab[key].append(title+' (out)')
                        else : lab[key].append(zone+' to ' +str(zname[key-1]))
                    else : 
                        ind[key].append(False)
                        lab[key].append(False)
        return ind, lab

############################ Mass balance Mt3dms ##############################
        
    def readMT3Dfile(self,graph,inlist,outlist):
      ### read PHT3D.XMAS file and put data in dict td   
      tr=2
    
############################# Mass balance Pht3d ##############################  

    def readPHT3Dfile(self,graph,inlist,outlist,splist):
      ### read PHT3D.XMAS file and put data in dict td
        td={'title':[],'data':[]}
        sind=format((int(splist[0])+1),'03d');  #print('sind',sind)
        fName=self.core.fileDir+os.sep+'PHT3D'+str(sind)+'.XMAS'
        file = np.loadtxt(fName,skiprows = 2)
        td['data']=file
        f1=open(fName,'r'); #utils//
        tiCol=[]
        title=list(filter(None, f1.readline().split('  ')))
        for i in range(len(title)):
            if title[i]==' HEAD-DEPENDENT BOUNDARY':
                title[i]='HEAD DEP BOUNDS'
            if i == 0: tiCol.append(title[i])
            else : 
                tiCol.append(title[i].strip()+' (IN)')
                tiCol.append(title[i].strip()+' (OUT)')
        f1.close()
        td['title']=tiCol
        inout=list(set(inlist+outlist))
      
        ### put data in dict xy following type of graph
        
        if graph == 'Percent Discrepency' :
            xy={'lab':[],'x':[],'y':[]}
            for row in td['data']:
                vin=row[td['title'].index('TOTAL (IN)')]
                vout=row[td['title'].index('TOTAL (OUT)')]
                pdisc=(float(vin)+float(vout))/(0.5*(float(vin)-float(vout)))*100
                xy['y'].append(pdisc)
                xy['x'].append(row[td['title'].index('TIME')])
            #print('Pd',xy)
                
        if graph == 'Time Series':
            xy={'lab':[],'x':[],'yin':[],'yout':[],'Cin':[],'Cout':[]}           
            xy['x']=td['data'][:,0]
            i=0 ; labin,labout=[],[]
            for title in td['title']:
                for l1 in inout:
                    if l1 in inlist:
                        if l1 in title and 'IN' in title :
                            ind=td['title'].index(title)
                            xy['yin'].append(td['data'][:,ind])
                            labin.append(title)
                            xy['Cin'].append('C'+str(i))
                    if l1 in outlist:
                        if l1 in title and 'OUT' in title :
                            ind2=td['title'].index(title)
                            xy['yout'].append(td['data'][:,ind2])
                            labout.append(title)
                            xy['Cout'].append('C'+str(i))
                            i+=1
                    else : i+=1
            xy['lab']=labin+labout            
            #print('Tseries',xy)
            
        if graph == 'Time Step':
            xy={'lab':[],'yin':[],'yout':[]} 
            t=self.Tstep.currentText()
            for i,row in enumerate(td['data']):
                if float(row[td['title'].index('TIME')])==float(t):
                    tind=i
            for title in td['title']:       
                for l1 in inout: 
                    if l1 in inlist and l1 in outlist:
                        if l1 in title and 'IN' in title :
                            ind=td['title'].index(title)
                            xy['yin'].append(td['data'][tind][ind])
                        if l1 in title and 'OUT' in title :
                            ind=td['title'].index(title)
                            xy['yout'].append(-td['data'][tind][ind])
                            xy['lab'].append(title.replace(' (OUT)',''))
                    else :
                        if l1 in outlist:
                            if l1 in title and 'OUT' in title :
                                ind=td['title'].index(title)
                                xy['yout'].append(-td['data'][tind][ind])
                                xy['lab'].append(title.replace(' (OUT)',''))
                                xy['yin'].append(np.nan)
                        if l1 in inlist:
                            if l1 in title and 'IN' in title :
                                ind=td['title'].index(title)
                                xy['yin'].append(td['data'][tind][ind])
                                xy['lab'].append(title.replace(' (IN)',''))
                                xy['yout'].append(np.nan) 
            #print('Tstep',xy)
        return xy