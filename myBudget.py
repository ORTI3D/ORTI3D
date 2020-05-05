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
#import seaborn as sns

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
        self.frame.setMaximumSize(QtCore.QSize(250, 85)) 
        self.gl = QGridLayout(self.frame)
    ## Different type of graph (frame 1)
        self.label_1 = QtWidgets.QLabel(self.frame)
        self.label_1.setText("Type of graph")
        self.gl.addWidget(self.label_1,0,0,1,1)
        self.plgroup = QComboBox(self)
        tgraph=['Percent Discrepancy','In-Out','Time Series','Time Step']
        if self.res != 'Flow' and self.typ=='M' : 
            tgraph=['Percent Discrepancy','Time Series','Time Step']
        if self.res != 'Flow' and self.typ=='Z' : 
            tgraph=['Time Series','Time Step']
        self.plgroup.addItems(tgraph)
        self.plgroup.setCurrentIndex(0)
        self.plgroup.activated['QString'].connect(self.onTstep)
        self.gl.addWidget(self.plgroup,0,1,1,1)
    ## Cumulative or not
        self.label_12 = QtWidgets.QLabel(self.frame)
        self.label_12.setText("View")
        self.gl.addWidget(self.label_12,1,0,1,1)
        self.vgroup = QComboBox(self)
        if self.res != 'Flow' : view=['Cumulative mass','Rate']
        else : view=['Cumulative volume','Rate']
        self.vgroup.addItems(view)
        if self.res != 'Flow' : self.vgroup.setCurrentIndex(0)
        else : self.vgroup.setCurrentIndex(1)
        self.vgroup.activated['QString'].connect(self.onView)
        self.gl.addWidget(self.vgroup,1,1,1,1)
        #self.verticalLayout.addWidget(self.frame)
        #if self.typ == 'Z' and self.res != 'Flow':
            #self.frame.hide()
    ## Species for chemistry 
        if self.res == 'Chemistry':
            self.label_13 = QtWidgets.QLabel(self.frame)
            self.label_13.setText("Species")
            self.gl.addWidget(self.label_13,2,0,1,1)
            self.cgroup = QComboBox(self)
            self.species=self.getSpecies()
            self.cgroup.addItems(self.species)
            self.gl.addWidget(self.cgroup,2,1,1,1)
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
        self.lzname=self.core.diczone['Observation'].dic['obs.1']['name']
        if self.typ == 'Z' and self.res == 'Flow':
    ## frame 3
            self.frame3 = QtWidgets.QFrame(self)
            self.frame3.setMaximumSize(QtCore.QSize(250, 35)) 
            self.gl3 = QGridLayout(self.frame3)
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
        self.pushButton2.clicked.connect(self.onExport)
    ## add vertical layout 2
        self.horizontalLayout.addLayout(self.verticalLayout2)
        QMetaObject.connectSlotsByName(self)  #OA 1/6/19

############################ Initialize the dialog ############################
        
    def onTstep(self):
    ## show or hide time combobox for the time step graph
        if self.plgroup.currentText()=='Time Step':
            self.frame2.show()
        else : self.frame2.hide()
    
    def onView(self):
    ## Recalculate if self.xy exist (if "apply" has been already clicked)
        try : 
            self.xy != None
            self.buildGraph()
        except : pass

    def getObsZone(self):
    ## get the names of model observation zone
        zname=self.lzname#self.core.diczone['Observation'].dic['obs.1']['name']
        dicZin={'Zones':{}};dicZout={'Zones':{}}
        if self.res=='Flow' :
            zone = self.zgroup.currentText()
            zIn=[zone+' from '+zname[i] for i in range(len(zname)) 
                if zname[i] != zone]
            zIn.append('TOTAL IN')
            dicZin['Zones'] = list(zip(zIn,[False]*len(zIn)))
            zOut=[zone+' to '+zname[i] for i in range(len(zname))
                if zname[i] != zone]
            zOut.append('TOTAL OUT')
            dicZout['Zones'] = list(zip(zOut,[False]*len(zOut)))
        else : dicZin['Zones'] = list(zip(zname,[False]*len(zname)))
        return dicZin, dicZout
    
    def getBoundaries(self):
    ## get the boundaries in the model (well, drn, storage, chd...)
        dic={'Bound':[]} ; nbound=['STORAGE','CONSTANT HEAD']
        lmod=self.core.getUsedModulesList('Mt3dms')
        if self.res == 'Chemistry':
            rec=['REACTION (PHREEQC)','MASS STORAGE (SOLUTE)']
            if 'RCT' in lmod:
                rec.append('MASS STORAGE (SORBED)')
            nbound=nbound + rec
        boundTyp=['WELLS','DRAINS','RIVER LEAKAGE','ET','HEAD DEP BOUNDS',
                  'RECHARGE']
        pack=['wel','drn','riv','evt','ghb','rch']
        lmodf=self.core.getUsedModulesList('Modflow')
        for i, n in enumerate (pack) :
            if n.upper() in lmodf:
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
        #dic={'Species':[]} 
        species = self.core.addin.chem.getListSpecies() 
        #dic['Species']=list(zip(species,[False]*len(species)))
        return species[:-2]#dic
    
    def getChoices(self):
    ## return a dic in function of type of graph and result to plot
        dic={'In':[],'Out':[]}
        if self.typ == 'M': 
            if self.res!='Transport':
                self.dicBound = self.getBoundaries()
                dic['In']=self.dicBound['Bound']+[('TOTAL IN',False)]
                dic['Out']=self.dicBound['Bound']+[('TOTAL OUT',False)]
            else :
                self.dicBound = self.getMt3dBound()
                dic['In / Out']=self.dicBound['Bound']
        else :
            dicZin,dicZout=self.getObsZone()
            if self.res=='Flow':
                self.dicBound = self.getBoundaries()
                dic['In']=self.dicBound['Bound']+dicZin['Zones']
                dic['Out']=self.dicBound['Bound']+dicZout['Zones']
            else : 
                dic['Zones']=dicZin['Zones']
        #if self.res =='Chemistry': 
            #dicSpecies=self.getSpecies()
            #dic = {**dicSpecies,**dic} 
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
               'iolist':[],'zlist':[],'splist':[]} # 'lylist':{}
        dicIn['ptyp']=self.typ
        dicIn['graph']=str(self.plgroup.currentText())
        if self.typ == 'Z' and self.res == 'Flow' : 
            dicIn['zone'] = self.zgroup.currentText()
        dic=self.getValues() 
        #print('dic',dic)
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
        if self.res=='Transport' and self.typ=='M':
            for i in range(len(dic['In / Out'])):
                if dic['In / Out'][i][1]==2:
                    dicIn['iolist'].append(dic['In / Out'][i][0])
        if self.res in ['Transport','Chemistry'] and self.typ=='Z':
            for i in range(len(dic['Zones'])):
                if dic['Zones'][i][1]==2:
                    dicIn['zlist'].append(dic['Zones'][i][0])
        if self.res =='Transport' :dicIn['splist']=['Tracer']
        if self.res=='Chemistry' :
            dicIn['splist']=self.cgroup.currentIndex()
            #dicIn['splist']=[dic['Species'].index(dic['Species'][i]) 
            #for i in range(len(dic['Species'])) if dic['Species'][i][1]==2] 
        #print('dicIn',dicIn)
        return dicIn
    
    def buildGraph(self):
        dicIn=self.getOptions()
        self.ptyp,self.graph=dicIn['ptyp'],dicIn['graph']
        self.zone,self.splist=dicIn['zone'],dicIn['splist']
        self.inlist,self.outlist= dicIn['inlist'],dicIn['outlist']
        self.iolist=dicIn['iolist'] ; self.zlist=dicIn['zlist']
        ## Zone budget
        if dicIn['ptyp'] == 'Z':
            if self.res == 'Flow':
                self.writeZoneFile(self.inlist,self.outlist,self.zone)
                self.writeInFile()
                self.core.runZonebud()
                self.xy=self.readZBfile(self.graph,self.zone,self.inlist,
                              self.outlist)
            else : 
                self.xy=self.getZoneMass(self.graph,self.zlist,self.splist)
        ## Mass balance
        else :
            if self.res == 'Chemistry' : 
                self.xy=self.readPHT3Dfile(self.graph,self.inlist,
                                           self.outlist,self.splist)
            if self.res == 'Transport' :
                self.xy=self.readMT3Dfile(self.graph,self.iolist)
            if self.res == 'Flow' :
                self.zone=[]
                self.writeZoneFile(self.inlist,self.outlist,self.zone)
                self.writeInFile()
                self.core.runZonebud()
                self.xy=self.readZBfile(self.graph,self.zone,self.inlist,
                          self.outlist)
        #print('xy',self.xy)
        self.plotData(self.xy,self.graph)
    
    def plotData(self,xy,graph):
        '''Plot data as scatter plot or vertical barchart'''
        self.figure.clf()
        #sns.set()
        self._ax=self.figure.add_subplot(1,1,1)
      ### Scatter plot
        if graph!='Time Step' :
            if graph!='Time Series' : ## % discrepency & in-out
                self._ax.plot(xy['x'],xy['y'],c='red',marker='o')
            else :
                if self.res=='Transport' and self.typ=='M': ## MB - transport
                    for i, color in enumerate(xy['C']):
                        self._ax.plot(xy['x'],xy['y'][i],c=color,
                                      marker=xy['m'][i])
                else: ## ZB - flow, trans, chem & MB - flow, chem
                    for i, color in enumerate(xy['Cin']):
                        self._ax.plot(xy['x'],xy['yin'][i],c=color,marker='v')
                    for i, color in enumerate(xy['Cout']):
                        self._ax.plot(xy['x'],xy['yout'][i],c=color,marker='^')
            self._ax.legend(xy['lab'])
            #self._ax.set_title(self.zolist[i],fontweight="bold", size=9)
            #self._ax.legend(self.llabel,fontsize = 8,loc='best')
            self._ax.set_ylabel(xy['ylab'], fontsize = 8) 
            self._ax.set_xlabel(xy['xlab'], fontsize = 8)
            self._ax.ticklabel_format(useOffset=False, 
                                      style='sci',scilimits=(-4,4),axis='both'
                                      ,useMathText=True)
            self._ax.tick_params(axis='both', labelsize=8)
       ### Vertical barchart
        if graph=='Time Step' :
            x = np.arange(len(xy['lab']))  # the label locations
            width = 0.35  # the width of the bars
            if self.res=='Transport' and self.typ=='M': ## MB - transport
                bar = self._ax.bar(x,xy['y'],color=['C'+str(c) for c in x])
                self.autolabel(bar)
            elif self.res in ['Transport','Chemistry'] and self.typ=='Z': ## ZB - chem & transport
                if xy['ysorb'] :x2=x + width/2 ; x=x - width/2 
                if xy['lab2'] : label = xy['lab2'][0] #; print('lab',label)
                barsol = self._ax.bar(x, xy['ysol'], width)#, label=label)
                self.autolabel(barsol)
                if xy['ysorb'] :
                    if xy['lab2'] : label = xy['lab2'][1]
                    barsor = self._ax.bar(x2, xy['ysorb'], width)#, label=label)
                    self.autolabel(barsor)
                if self.xy['lab']:self._ax.legend(label)
            else : ## ZB - flow & MB - Flow, chem
                barIn = self._ax.bar(x - width/2, xy['yin'], 
                                     width, label='IN')
                barOut = self._ax.bar(x + width/2, xy['yout'], 
                                      width, label='OUT')
                self.autolabel(barIn)
                self.autolabel(barOut)
                if self.xy['lab']:self._ax.legend()
            #self._ax.set_ylabel('')
            #self._ax.set_title('')
            self._ax.set_xticks(x)
            self._ax.set_xticklabels(xy['lab'])
            self._ax.set_ylabel(xy['ylab'], fontsize = 8) 
        self._ax.figure.canvas.draw()
                        
    def autolabel(self,rects):
        '''Attach a text label above each bar in *rects*, 
        displaying its height'''
        for rect in rects:
            height = rect.get_height()#round(rect.get_height(),1) # ;print(height)
            self._ax.annotate('{0:.2g}'.format(height),
                        xy=(rect.get_x() + rect.get_width() / 2, height),
                        xytext=(0, 3),  # 3 points vertical offset
                        textcoords="offset points",
                        ha='center', va='bottom')
            
####################### Zone budget & Mass balance flow #######################        
    
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
        if self.typ != 'M':
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
            Percent Discrepancy: Percent error vs time
            In-Out: In - out vs time
            Time Series: In Zone or BC and Out Zone or BC vs time
            Time Step: In Zone or BC and Out Zone or BC for 1 time
            Return xy dictionnary'''
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
        zname=self.lzname
        if self.typ=='M':idz=0
        else : idz = zname.index(zone) ## index on selected zone for zone budget
        utime=self.core.getUnits('Modflow','dis.8',0)[:-1]
        if utime == '': utime='T'
        ulength=self.core.getUnits('Modflow','dis.4',0)
        if ulength=='':ulength='L'
        cview=self.vgroup.currentIndex()
        #print('idz+1',str(idz+1))
        tlist=np.insert(self.tlist,0,0) ## Cumulative or not
        tsteps=[tlist[i+1] - tlist[i] for i in range(len(tlist)-1)] ## Cumulative or not
        
        if graph == 'Percent Discrepancy' :
            for row in td['data']:
                if row[td['title'].index('ZONE')]==str(idz+1):
                    xy['x'].append(float(row[td['title'].index('TOTIM')]))
                    xy['y'].append(float(row[td['title'].index('Percent Error')]))
            xy['lab'].append('Percent Discrepancy')
            xy['ylab']='Discrepancy (%)'
            xy['xlab']='Time ('+utime+')'
            #print ('perc dis',xy)
        
        if graph == 'In-Out':
            y=[]
            for row in td['data']:
                if row[td['title'].index('ZONE')]==str(idz+1):
                    xy['x'].append(float(row[td['title'].index('TOTIM')]))
                    #xy['y'].append(float(row[td['title'].index('IN-OUT')]))
                    y.append(float(row[td['title'].index('IN-OUT')]))
            if cview != 0 : xy['y'] = y
            else:
                for t in range(len(y)):
                    if t==0 : y0=y[t]*tsteps[t]
                    else:y0=y0+y[t]*tsteps[t]
                    xy['y'].append(y0)
            xy['lab'].append('In-Out')
            if cview != 0 : xy['ylab']='Rates ('+ulength+'\u00b3/'+utime+')' 
            else : xy['ylab']='Volume ('+ulength+'\u00b3)'
            xy['xlab']='Time ('+utime+')'
            #print ('in-out',xy)
        
        if graph in ['Time Series','Time Step']:
            ind,lab=self.getIndLab(td,inlist,outlist,zone)
            xy={'lab':[],'x':[],'yin':[],'yout':[],'Cin':[],'Cout':[]}
            i,j,k=0,0,0; lin=[]; lout=[] ## i&j for color, k for xy['x'], lin&lout for legend
            for key, value in ind.items(): ## loop in user choices
                yi=[];yo=[];y2=[]  ## y2 array for cumulative
                if value[0]: ## for sources (in)
                    for row in td['data']: ## loop in zone budget out file
                        if row[td['title'].index('ZONE')]==str(idz+1): ## selected zone for zone budget
                            if k==0: xy['x'].append(float(
                                    row[td['title'].index('TOTIM')])) ## get time
                            yi.append(float(row[int(value[0])]))
                    if cview != 0 :xy['yin'].append(yi) ## not cumulative (default)
                    else : ## cumulative
                        for t in range(len(yi)):
                            if t==0 : y0=yi[t]*tsteps[t]
                            else:y0=y0+yi[t]*tsteps[t]
                            y2.append(y0)
                        xy['yin'].append(y2)
                    xy['Cin'].append('C'+str(i))
                    lin.append(lab[key][0])
                    k+=1
                i+=1
                y2=[]  ## y2 array for cumulative
                if value[1]: ## for sinks (out)
                    for row in td['data']:
                        if row[td['title'].index('ZONE')]==str(idz+1):
                            if k==0: xy['x'].append(float(
                                    row[td['title'].index('TOTIM')]))
                            yo.append(float(row[int(value[1])]))
                    if cview != 0 :xy['yout'].append(yo)
                    else :  
                        for t in range(len(yo)):
                            if t==0 : y0=yo[t]*tsteps[t]
                            else:y0=y0+yo[t]*tsteps[t]
                            y2.append(y0)
                        xy['yout'].append(y2)
                    xy['Cout'].append('C'+str(j))
                    lout.append(lab[key][1])
                    k+=1
                j+=1
            xy['lab']=lin+lout
            if cview != 0 : xy['ylab']='Rates ('+ulength+'\u00b3/'+utime+')'
            else :  xy['ylab']='Volume ('+ulength+'\u00b3)'
            xy['xlab']='Time ('+utime+')'
            #print ('TSeries',xy)
            
            if graph == 'Time Step':
                xy2={'lab':[],'yin':[],'yout':[]}
                i,j,k=0,0,0
                tind=self.Tstep.currentIndex()
                for key, value in ind.items(): 
                    if value[0]:
                        xy2['yin'].append(xy['yin'][j][tind])
                        j+=1
                    else : xy2['yin'].append(np.nan)
                    if value[1]:
                        xy2['yout'].append(xy['yout'][k][tind])
                        k+=1
                    else: xy2['yout'].append(np.nan)
                    if str(list(lab.keys())[i]).isdigit():
                        xy2['lab'].append(zname[list(lab.keys())[i]-1])
                    else:xy2['lab'].append(list(lab.keys())[i])
                    i+=1
                if cview != 0:xy2['ylab']='Rates ('+ulength+'\u00b3/'+utime+')'
                else :  xy2['ylab']='Volume ('+ulength+'\u00b3)'
                return xy2
        return xy
        '''
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
            xy['ylab']='Rates ('+ulength+'$^{3}$/'+utime+')'
            #print ('TStep',xy)'''
    
    def getIndLab(self,td,inlist,outlist,zone):
        ''' Return the label and column index for each user selected zone 
        and BC. Organized in pairs, 'in' and 'out' index or label are indicated 
        for each zone or BC. Takes the value False if only 'in' or only 'out' 
        is selected for a BC or a zone.
        '''
        inout= list(set(inlist+outlist))
        dicb={}; inT,outT=False,False
        for i in inout :
            if i in ['TOTAL IN','TOTAL OUT']: #'Total IN','Total OUT',
                if i in ['TOTAL IN']:inT=True
                if i in ['TOTAL OUT']:outT=True
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
                        if str(key)== 'Total': lab[key].append('TOTAL IN')
                        elif len(title.split('FROM '))==1: 
                            lab[key].append(title+' IN')
                        else : lab[key].append(zone+' from '+str(zname[key-1]))
                    else : 
                        ind[key]=[False]
                        lab[key].append(False)
                    i+=1 ; continue
                if str(key) in title and i==1 :
                    if value[1]==True: 
                        ind[key].append(j)
                        if str(key)== 'Total': lab[key].append('TOTAL OUT')
                        elif len(title.split('TO '))==1: 
                            lab[key].append(title+' OUT')
                        else : lab[key].append(zone+' to ' +str(zname[key-1]))
                    else : 
                        ind[key].append(False)
                        lab[key].append(False)
        return ind, lab

######################### Zone budget Mt3dms & Pht3d ##########################

    def getZoneMass(self,graph,zlist,splist):
        group=self.res ; zlayers='all' ; iper = 0 ; 
        lmod=self.core.getUsedModulesList('Mt3dms')
        utime=self.core.getUnits('Modflow','dis.8',0)[:-1]
        if utime == '': utime='T'
        if group == 'Chemistry' : 
            sp=self.species[splist] ; umass ='mol' 
        else : sp=splist ; umass='M'
        cview=self.vgroup.currentIndex()
        tlist=np.insert(self.tlist,0,0) ## Cumulative or not
        tsteps=[tlist[i+1] - tlist[i] for i in range(len(tlist)-1)] ## Cumulative or not
        
        if graph == 'Time Series':
            lin=[];lout=[];x=''
            xy={'lab':[],'x':[],'yin':[],'yout':[],'Cin':[],'Cout':[]}
            for i, zname in enumerate(zlist):
                #print('B4',iper,group,zname,[sp],zlayers)
                x,ysol,label =  self.core.onPtObs('B4',iper,group,zname,
                                                  [sp],zlayers,ss='')
                if cview!=0 :
                    ysol=ysol[:,0]
                    ysol[1:] -= ysol[:-1].copy()
                    ysol=ysol/tsteps
                    xy['yin'].append(ysol)
                else : xy['yin'].append(ysol[:,0])
                if group == 'Chemistry' : 
                    lin.append(zname+'_'+str(sp)+' (SOLUTE)')
                else : lin.append(zname+'_Tracer (SOLUTE)')
                xy['Cin'].append('C'+str(i))
                if 'RCT' in lmod:
                    x,ysor,label =  self.core.onPtObs('B4',iper,group,zname,
                                                      [sp],zlayers,ss='S')
                    if cview!=0 :
                        ysor=ysor[:,0]
                        ysor[1:] -= ysor[:-1].copy()
                        ysor=ysor/tsteps
                        xy['yout'].append(ysor)
                    else : xy['yout'].append(ysor[:,0])
                    if group == 'Chemistry' : 
                        lin.append(zname+'_'+str(sp)+' (SORBED)')
                    else : lin.append(zname+'_Tracer (SORBED)')
                    xy['Cout'].append('C'+str(i))
            xy['x']=x
            xy['lab']=lin+lout
            if cview != 0 : xy['ylab']='Rates ('+umass+'/'+utime+')'
            else :  xy['ylab']='Mass ('+umass+')'
            xy['xlab']='Time ('+utime+')'
            #print('time series',xy)
        
        if graph == 'Time Step':
            #t=self.Tstep.currentText()
            tind=self.Tstep.currentIndex()
            t=self.tlist[tind] ; t0=self.tlist[tind-1]
            if tind == 0 : t0=0
            tstep=t-t0
            xy={'lab':[],'lab2':[],'ysol':[],'ysorb':[]}
            for i, zname in enumerate(zlist):
                x,ysol,label =  self.core.onPtObs('B4',iper,group,zname,
                                            [sp],zlayers,ss='')
                ind = list(x).index(float(t))
                if cview!=0 : 
                    y=ysol[:,0]
                    y[1:] -= y[:-1].copy()
                    y=y/tstep
                else : y=ysol[:,0]
                xy['ysol'].append(y[ind])
                xy['lab'].append(zname)
                if group == 'Chemistry' : 
                    xy['lab2'].append([str(sp)+' (SOLUTE)'])
                else : xy['lab2'].append(['Tracer (SOLUTE)'])
                if 'RCT' in lmod:
                    x,ysor,label =  self.core.onPtObs('B4',iper,group,zname,
                                                [sp],zlayers,ss='S')
                    if cview!=0 :
                        y=ysor[:,0]
                        y[1:] -= y[:-1].copy()
                        y=y/tstep
                    else : y=ysor[:,0]
                    xy['ysorb'].append(y[ind])
                    if group == 'Chemistry' : 
                        xy['lab2'].append([str(sp)+' (SORBED)'])
                    else :xy['lab2'].append(['Tracer (SORBED)'])
            if cview != 0 : xy['ylab']='Rates ('+umass+'/'+utime+')'
            else :  xy['ylab']='Mass ('+umass+')'
            #print('t step',xy)
        
        return xy

############################ Mass balance Mt3dms ##############################
        
    def readMT3Dfile(self,graph,iolist):
      ### read MT3D.MAS file and put data in dict td   
        td={'title':[],'data':[]}
        td['title']=['TIME','TOTAL IN','TOTAL OUT','SOURCES IN','SINKS OUT',
                   'TOTAL MASS IN AQUIFER','Percent Discrepancy']
        fName=self.core.fileDir+os.sep+'MT3D001.MAS'
        td['data'] = np.loadtxt(fName,skiprows = 2)
        td['data']=np.delete(td['data'], [5,8], 1) 
        utime=self.core.getUnits('Modflow','dis.8',0)[:-1]
        if utime == '': utime='T'
        umass='M'
        cview=self.vgroup.currentIndex()
        
        ### put data in dict xy following type of graph
        
        if graph == 'Percent Discrepancy' :
            xy={'lab':[],'x':[],'y':[]}
            xy['y']=td['data'][:,6]
            xy['x']=td['data'][:,0]
            xy['lab'].append('Percent Discrepancy')
            xy['ylab']='Discrepancy (%)'
            xy['xlab']='Time ('+utime+')'
            #print('pd',xy)
        
        if graph == 'Time Series':
            xy={'lab':[],'x':[],'y':[],'C':[],'m':[]}  
            C=['C0','C0','C1','C1','C2']  
            m=['v','^','v','^','o']       
            xy['x']=td['data'][:,0]
            tlist=np.insert(xy['x'],0,0)
            tsteps=[tlist[i+1] - tlist[i] for i in range(len(tlist)-1)]
            for l1 in iolist:
                if l1 in td['title'] :
                    ind=td['title'].index(l1)
                    y=abs(td['data'][:,ind])
                    if cview!=0 : 
                        y[1:] -= y[:-1].copy()
                        y=y/tsteps
                    xy['y'].append(y)
                    xy['lab'].append(l1)
                    xy['C'].append(C[ind-1])
                    xy['m'].append(m[ind-1])
            if cview != 0 : xy['ylab']='Rates ('+umass+'/'+utime+')'
            else :  xy['ylab']='Mass ('+umass+')'
            xy['xlab']='Time ('+utime+')'
            #print('tser',xy)

        if graph == 'Time Step':
            xy={'lab':[],'y':[]}
            t=self.Tstep.currentText()
            for i,row in enumerate(td['data']):
                if float(row[td['title'].index('TIME')])==float(t):
                    tind=i
                    if tind == 0 : t0=0
                    else : t0 = float(td['data'][i-1,0])
                    tstep=float(t)-t0
            for l1 in iolist:
                if l1 in td['title'] :
                    ind=td['title'].index(l1)
                    y=abs(td['data'][:,ind])
                    if cview!=0 : 
                        y[1:] -= y[:-1].copy()
                        y=y/tstep
                    xy['y'].append(y[tind])
                    xy['lab'].append(l1)
            if cview != 0 : xy['ylab']='Rates ('+umass+'/'+utime+')'
            else :  xy['ylab']='Mass ('+umass+')'
            #print('tstep',xy)
        
        return xy
    
############################# Mass balance Pht3d ##############################  

    def readPHT3Dfile(self,graph,inlist,outlist,splist):
      ### read PHT3D.XMAS file and put data in dict td
        td={'title':[],'data':[]}
        sind=format((int(splist)+1),'03d');  #print('sind',sind)
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
                tiCol.append(title[i].strip()+' IN')
                tiCol.append(title[i].strip()+' OUT')
        f1.close()
        td['title']=tiCol
        inlist=[i.replace('TOTAL IN', 'TOTAL') for i in inlist]
        outlist=[i.replace('TOTAL OUT', 'TOTAL') for i in outlist]
        inout=list(set(inlist+outlist))
        
        utime=self.core.getUnits('Modflow','dis.8',0)[:-1]
        if utime == '': utime='T'
        umass ='mol'
        cview=self.vgroup.currentIndex()
      
        ### put data in dict xy following type of graph
        
        if graph == 'Percent Discrepancy' :
            xy={'lab':[],'x':[],'y':[]}
            for row in td['data']:
                vin=row[td['title'].index('TOTAL IN')]
                vout=row[td['title'].index('TOTAL OUT')]
                pdisc=(float(vin)+float(vout))/(0.5*(float(vin)-float(vout)))*100
                xy['y'].append(pdisc)
                xy['x'].append(row[td['title'].index('TIME')])
            xy['lab'].append('Percent Discrepancy')
            xy['ylab']='Discrepancy (%)'
            xy['xlab']='Time ('+utime+')'
            #print('Pd',xy)
                
        if graph == 'Time Series':
            xy={'lab':[],'x':[],'yin':[],'yout':[],'Cin':[],'Cout':[]}           
            xy['x']=td['data'][:,0]
            i,j=0,0 ; labin,labout=[],[]
            tlist=np.insert(self.tlist,0,0)
            tsteps=[tlist[i+1] - tlist[i] for i in range(len(tlist)-1)]
            for title in td['title']:
                for l1 in inout:
                    if l1 in inlist:
                        if l1 in title and 'IN' in title :
                            ind=td['title'].index(title)
                            y=abs(td['data'][:,ind])
                            if cview != 0 :
                                y[1:] -= y[:-1].copy()
                                y=y/tsteps
                            xy['yin'].append(y)
                            labin.append(title)
                            xy['Cin'].append('C'+str(i))
                            j=i ; i+=1
                    if l1 in outlist:
                        if l1 in title and 'OUT' in title :
                            ind2=td['title'].index(title)
                            y=abs(td['data'][:,ind2])
                            if cview != 0 :
                                y[1:] -= y[:-1].copy()
                                y=y/tsteps
                            xy['yout'].append(y)
                            labout.append(title)
                            xy['Cout'].append('C'+str(j))
                            j+=1
                            if j>i: i+=1
            xy['lab']=labin+labout
            if cview != 0 : xy['ylab']='Rates ('+umass+'/'+utime+')'
            else :  xy['ylab']='Mass ('+umass+')'
            xy['xlab']='Time ('+utime+')'            
            #print('Tseries',xy)
            
        if graph == 'Time Step':
            xy={'lab':[],'yin':[],'yout':[]} 
            tind=self.Tstep.currentIndex()
            t=self.tlist[tind] ; t0=self.tlist[tind-1]
            if tind == 0 : t0=0
            tstep=t-t0
            for title in td['title']:       
                for l1 in inout: 
                    if l1 in inlist and l1 in outlist:
                        if l1 in title and 'IN' in title :
                            ind=td['title'].index(title)
                            y=abs(td['data'][:,ind])
                            if cview != 0 : 
                                y[1:] -= y[:-1].copy()
                                y=y/tstep
                            xy['yin'].append(y[tind])
                        if l1 in title and 'OUT' in title :
                            ind=td['title'].index(title)
                            y=abs(td['data'][:,ind])
                            if cview != 0 : 
                                y[1:] -= y[:-1].copy()
                                y=y/tstep
                            xy['yout'].append(y[tind])
                            xy['lab'].append(title.replace(' (OUT)',''))
                    else :
                        if l1 in outlist:
                            if l1 in title and 'OUT' in title :
                                ind=td['title'].index(title)
                                y=abs(td['data'][:,ind])
                                if cview != 0 : 
                                    y[1:] -= y[:-1].copy()
                                    y=y/tstep
                                xy['yout'].append(y[tind])
                                xy['lab'].append(title.replace(' (OUT)',''))
                                xy['yin'].append(np.nan)
                        if l1 in inlist:
                            if l1 in title and 'IN' in title :
                                ind=td['title'].index(title)
                                ind=td['title'].index(title)
                                y=abs(td['data'][:,ind])
                                if cview != 0 : 
                                    y[1:] -= y[:-1].copy()
                                    y=y/tstep
                                xy['yin'].append(y[tind])
                                xy['lab'].append(title.replace(' (IN)',''))
                                xy['yout'].append(np.nan) 
            if cview != 0 : xy['ylab']='Rates ('+umass+'/'+utime+')'
            else :  xy['ylab']='Mass ('+umass+')'
            #print('Tstep',xy)
        return xy

############################### Export feature ################################  

    def onExport(self):
        dlg = myFileDialog('Save')
        fDir,fName = dlg.getsetFile(self.gui,'Save','*.txt')
        if fDir == None: return
        f1 = open(fDir+os.sep+fName+'.txt','w', encoding="utf-8")
        if self.ptyp  == 'Z':
            if self.res == 'Flow':
                f1.write('Zone budget zone: '+self.zone+'\n')
        if self.graph!='Time Step' :
            #print(self.xy['lab'])
            lab = '\t'.join(self.xy['lab'])
            if self.graph!='Time Series' : ## % discrepency & in-out
                f1.write('Time\t'+lab+'\n') ##lab
                f1.write(self.xy['xlab'].split()[1]+'\t'+
                         str(self.xy['ylab'])+'\n') ## unit
                arr=self.xy['x']+self.xy['y']
                arr=np.array(arr).reshape(2,len(self.xy['x']))
                savetxt(f1,arr.T)
            else :
                if self.res=='Transport' and self.typ=='M': ## MB - transport
                    f1.write('Time\t'+lab+'\n') ##lab
                    ylab='\t'.join([str(self.xy['ylab'])]*len(self.xy['lab']))
                    f1.write(self.xy['xlab'].split()[1]+'\t'+str(ylab)+'\n') ## unit
                    arr=[self.xy['x']]+self.xy['y']
                    arr=np.array(arr).reshape(
                            len(self.xy['lab'])+1,len(self.xy['x']))
                    savetxt(f1,arr.T)
                else :  ## ZB - flow, trans, chem & MB - flow, chem
                    f1.write('Time\t'+lab+'\n') ##lab
                    ylab='\t'.join([str(self.xy['ylab'])]*len(self.xy['lab']))
                    f1.write(self.xy['xlab'].split()[1]+'\t'+str(ylab)+'\n') ## unit
                    arr=[self.xy['x']]+self.xy['yin']+self.xy['yout']
                    arr=np.array(arr).reshape(
                            len(self.xy['lab'])+1,len(self.xy['x']))
                    savetxt(f1,arr.T)
        else : ## Time Step
            #print('lab',self.xy['lab'])
            lab = '\t'.join(self.xy['lab'])
            f1.write('Time: '+self.Tstep.currentText()+'\n')
            if self.res=='Transport' and self.typ=='M': ## MB - transport
                ylab='\t'.join([str(self.xy['ylab'])]*len(self.xy['lab']))
                f1.write(lab+'\n')
                f1.write(ylab+'\n')
                savetxt(f1,[self.xy['y']])
            elif self.res in ['Transport','Chemistry'] and self.typ=='Z': ## ZB - chem & transport
                if self.xy['ysorb'] : 
                    a=list(zip(self.xy['ysol'],self.xy['ysorb']))
                    a1= [y for x in a for y in x]
                    i=2 ; form=['solute','sorbed']
                else : a1=self.xy['ysol'] ; i=1 ; form=['solute']
                f1.write(lab+'\n')
                ylab='\t'.join([str(self.xy['ylab'])]*len(self.xy['lab'])*i)
                out='\t'.join(form*len(self.xy['lab']))
                f1.write(out+'\n')
                f1.write(ylab+'\n')
                np.savetxt(f1,[a1])
            else : ## ZB - flow & MB - Flow, chem
                llab,ly=[],[]
                for i in range (len(self.xy['lab'])):
                    if str(self.xy['yin'][i]) != 'nan' :
                        llab.append(self.xy['lab'][i]+'(in)')
                        ly.append(self.xy['yin'][i])
                    if str(self.xy['yout'][i])!= 'nan' :
                        llab.append(self.xy['lab'][i]+'(out)')
                        ly.append(self.xy['yout'][i])
                lab = '\t'.join(llab)
                f1.write(lab+'\n') 
                ylab='\t'.join([str(self.xy['ylab'])]*len(llab))
                f1.write(ylab+'\n')  
                np.savetxt(f1,[ly])
        f1.close()
            
 

"""  
            if self.res=='Transport' and self.typ=='M': ## MB - transport
                bar = self._ax.bar(x,xy['y'],color=['C'+str(c) for c in x])
                self.autolabel(bar)
            elif self.res in ['Transport','Chemistry'] and self.typ=='Z': ## ZB - chem & transport
                if xy['ysorb'] :x2=x + width/2 ; x=x - width/2 
                barsol = self._ax.bar(x, xy['ysol'], width, label='SOLUTE')
                self.autolabel(barsol)
                if xy['ysorb'] :
                    barsor = self._ax.bar(x2, xy['ysorb'], width, label='SORBED')
                    self.autolabel(barsor)
                self._ax.legend()
            else : ## ZB - flow & MB - Flow, chem
                barIn = self._ax.bar(x - width/2, xy['yin'], 
                                     width, label='IN')
                barOut = self._ax.bar(x + width/2, xy['yout'], 
                                      width, label='OUT')
                
                
        if self.ptyp == 'XY':
            f1.write('Well '+'Type '+'Layer '+'Time '+'Observed '+'Simulated'+'\n')
            for n in range(len(self.label_all)): 
                    f1.write(self.label_all[n].split('_')[0]+' ')
                    f1.write(self.label_all[n].split('_')[1]+' ')
                    f1.write(self.label_all[n].split('_')[2][3:]+' ')
                    f1.write(str(self.obs_time[n])+' ')
                    f1.write(str(self.yobs_all[n])+' ')
                    f1.write(str(self.ysim_all[n])+'\n')
            f1.close()
        elif self.ptyp[0]=='B':
            f1.write(self.axlabel.split(' ')[0])
            if self.pOrder=='Zones': zslist=self.zolist
            else : zslist= self.splist
            for i in zslist:
                for n in self.llabel: 
                    if n.split('(')[1]!='obs)':f1.write(' '+i+'_'+n)
            f1.write('\n')
            nt,ny,nz = np.shape(self.arryy)
            if nt!=1 : 
                self.arryy=array(self.arryy).transpose(1,0,2).reshape(ny,-1) #EV 14/02/19
            else : 
                self.arryy=array(self.arryy).reshape(ny,-1) #EV 27/02/2019
            arr = zeros((ny,nt*nz+1,))
            arr[:,0]=self.x ; 
            arr[:,1:]=(self.arryy)
            savetxt(f1,arr)
            f1.close()
        elif self.ptyp[0] in ['P','V']: 
            for i in range(len(self.zolist)):
                if self.ptyp[0]=='P':f1.write(self.axlabel.split(' ')[0])
                else : f1.write(self.aylabel.split(' ')[0])
                for n in self.llabel: 
                    if n.split('(')[1]!='obs)':f1.write(' '+self.zolist[i]+'_'+n) 
                f1.write('\n')
                for j in range(len(self.arrx[i])):
                    arr=np.insert(self.arryy[i][j],0,self.arrx[i][j])
                    np.savetxt(f1,arr, newline=" ")
                    f1.write('\n')
                f1.write('\n')
                f1.write('\n')
            f1.close()
"""