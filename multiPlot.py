from matplotlib.backends.qt_compat import QtCore, QtWidgets
from matplotlib.backends.backend_qt5agg import (
        FigureCanvas, NavigationToolbar2QT as NavigationToolbar)
from matplotlib.figure import Figure
from .qtDialogs import *
from .geometry import *
#from .qtShow import * # OA 1/6/19
import numpy as np
import matplotlib.ticker as ticker

class multiPlot(QDialog):
    '''This dialog provides multiplot, there are mainly three types of plots :
    - Time series
    - Horizontal profiles
    - Calibration graph
    The result that can be represented are flow (head, Wcontent, Darcy flux), 
    Transport (Tracer concentration, flux and mass), Chemistry (chemical species 
    concentrations, flux and mass).
    They can also be compared to data.
    '''
    def __init__(self,gui,core,typ,res):
        self.gui,self.core= gui,core
        self.typ,self.res= typ,res
        QDialog.__init__(self,gui) # OA 6/11/18 added gui (self,gui)
        self.setModal(False)
        self.setWindowTitle('Plot of results')
        screenShape = QtWidgets.QDesktopWidget().screenGeometry()
        self.setGeometry(QRect(5, 5, screenShape.width()*.75, screenShape.height()*.7))
    ## main horizontal layout
        self.horizontalLayout = QHBoxLayout(self)
        self.horizontalLayout.setContentsMargins(10, 20, 10, 10)
    ## the left panel vertical layout  
        self.verticalLayout = QVBoxLayout()
        #self.verticalLayout.setGeometry(QRect(5, 5, 250, screenShape.height()*.68))
    ## title
        if self.typ == 'B' : title = 'Time-series graph'
        if self.typ == 'P' : title = 'Profile graph'
        if self.typ == 'X' : title = 'Calibration graph'
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
    ## frame 
        self.frame = QtWidgets.QFrame(self)
        self.frame.setMaximumSize(QtCore.QSize(250, 60)) 
        self.gl = QGridLayout(self.frame)
    ## type of result 
        if self.res != 'W content':
            self.label_0 = QtWidgets.QLabel(self.frame)
            #self.label_0.setMaximumSize(QtCore.QSize(40, 20))
            self.label_0.setText("Type of result")
            self.gl.addWidget(self.label_0,0,0,1,1)
            self.rgroup = QComboBox(self)
            if self.res == 'Flow' : self.rgroup.addItems(['Head','W content']) #EV 02/03/20
            if self.res in ['Transport','Chemistry'] : 
                self.rgroup.addItems(['Concentration','Weighted concentration'])
            self.rgroup.setCurrentIndex(0)
            self.gl.addWidget(self.rgroup,0,1,1,1)
    ## Plot order combo box for chemistry
        if self.res=='Chemistry' and self.typ=='B':
            self.label_1 = QtWidgets.QLabel(self.frame)
            #self.label_1.setMaximumSize(QtCore.QSize(40, 20))
            self.label_1.setText("Plot order")
            self.gl.addWidget(self.label_1,1,0,1,1)
            self.plgroup = QComboBox(self)
            self.plgroup.addItems(['By zone','By species'])
            self.plgroup.setCurrentIndex(0)
            self.gl.addWidget(self.plgroup,1,1,1,1)
            #!self.verticalLayout.addWidget(self.plgroup) 
    ## Time combo box for profile and calibration graph
        if self.typ=='P' or self.typ=='X':
            #!self.frame = QtWidgets.QFrame(self)
            #!self.frame.setMaximumSize(QtCore.QSize(120, 38))
            #!self.horizontalLayout2 = QHBoxLayout(self.frame)
            #!self.horizontalLayout2.setContentsMargins(0, 0, 0, 0)
            #!self.horizontalLayout2.setSpacing(0)
            self.label_2 = QtWidgets.QLabel(self.frame)
            self.label_2.setMaximumSize(QtCore.QSize(40, 20))
            self.label_2.setText("Time")
            self.gl.addWidget(self.label_2,1,0,1,1)
            #!self.horizontalLayout2.addWidget(self.label_2)
            if self.typ=='X':
                self.tlist=list(self.tlist)
                self.tlist.insert(0, 'All')
            self.Tstep = QComboBox(self.frame)
            self.Tstep.addItems([str(n) for n in self.tlist])
            self.Tstep.setCurrentIndex(0)
            self.gl.addWidget(self.Tstep,1,1,1,1)
            #!self.horizontalLayout2.addWidget(self.Tstep)
            #!self.verticalLayout.addWidget(self.frame) 
        self.verticalLayout.addWidget(self.frame) 
    ## choise plot obsData
        if self.typ=='B': # only for time series graphs for now (or typ=='P') 
            self.checkBox = QCheckBox(self)
            self.checkBox.setText("Observed Data")    
            self.verticalLayout.addWidget(self.checkBox)
    ## the options :need to go in the interface to search for zones and others
        dic=self.getChoices(self.res,self.typ)
        self.nb = myNoteBookCheck(self.gui,"Options",dic)
        self.nb.layout.removeWidget(self.nb.buttonBox) #EV 06/12
        self.nb.buttonBox.deleteLater()
        del self.nb.buttonBox
        #self.nb.setGeometry(QRect(5, 5, 250,270))
        self.verticalLayout.addWidget(self.nb)
        self.nb.apply()
    ## Apply button
        self.pushButton = QPushButton(self)
        self.pushButton.setText('Apply')
        self.verticalLayout.addWidget(self.pushButton, alignment=Qt.AlignHCenter)
        self.pushButton.clicked.connect(self.buildPlot)
     ## add vertical layout   
        self.horizontalLayout.addLayout(self.verticalLayout)
     ## the right panel vertical layout
        self.verticalLayout2 = QVBoxLayout()
     ## the matplotlib figure 
        self.figure = Figure(tight_layout=True,figsize=(7.8, 3), dpi=100) # EV 04/02/20 
        self.cnv = FigureCanvas(self.figure) 
        #self._ax = self.cnv.figure.subplots()#.add_axes([0.1, 0.15, 0.7, 0.8])
    ## add matplotlib figure
        #self.horizontalLayout.addWidget(self.cnv)
        self.verticalLayout2.addWidget(self.cnv)
        self.toolbar = NavigationToolbar(self.cnv, self) # EV 04/02/20 
        self.verticalLayout2.addWidget(self.toolbar, alignment=Qt.AlignHCenter) # EV 04/02/20 
    ## Export button
        self.pushButton2 = QPushButton(self)
        self.pushButton2.setText('Export')
        self.verticalLayout2.addWidget(self.pushButton2, alignment=Qt.AlignHCenter)
        self.pushButton2.clicked.connect(self.onExport)
    ## add vertical layout 2
        self.horizontalLayout.addLayout(self.verticalLayout2)
        QMetaObject.connectSlotsByName(self)  #OA 1/6/19
  
    def getObsZone(self):
    ## get the names of model observation zone
        dic={'Zones':{}}
        zname=self.core.diczone['Observation'].dic['obs.1']['name']
        dic['Zones'] = list(zip(zname,[False]*len(zname)))
        return dic
    
    def getLayers(self):
    ## get the number of model layers
        dic={'Layers':{}}
        nblay=getNlayers(self.core)
        lnblay =  [str(x) for x in range(nblay)]
        dic['Layers'] = list(zip(lnblay,[False]*nblay))
        return dic
        
    def getSpecies(self):
    ## get the names of model chemical species
        dic={'Species':{}} 
        species = self.core.addin.chem.getListSpecies() 
        dic['Species']=list(zip(species,[False]*len(species)))
        return dic
    
    def getChoices(self,res,typ):
    ## return a dic in function of type of graph and result to plot
        dicObsZone=self.getObsZone()
        dicLayers=self.getLayers()
        if res =='Chemistry': 
            dicSpecies=self.getSpecies()
            if typ=='X' or len(dicLayers['Layers'])==1: dic = {**dicObsZone,**dicSpecies} #EV 26/08/19
            else : dic = {**dicObsZone,**dicLayers,**dicSpecies}
        else : 
            if typ=='X' or len(dicLayers['Layers'])==1: dic = dicObsZone #EV 26/08/19
            else : dic = {**dicObsZone,**dicLayers}
        return dic
    
    def getValues(self): #EV 14/08/19
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
        '''get the plot options from the window very simple now'''
        dicIn={'ptyp':{},'plotOrder':{},'zolist':{},'splist':{},'lylist':{}} # 'lylist':{}
        ptyp=self.typ
        rtyp = int(self.rgroup.currentIndex()) #EV 02/03/20
        dicIn['ptyp']=ptyp+str(rtyp)
        if self.rgroup.currentText()=='W content': dicIn['ptyp']=ptyp+'0'
        #if ptyp=='B' : dicIn['ptyp']='B0'
        #if ptyp=='P' : dicIn['ptyp']='P0'
        #if ptyp=='X' : dicIn['ptyp']='XY'
        #dic=self.nb.getValues() 
        dic=self.getValues() #EV 14/08/19
        nblay=getNlayers(self.core) #EV 26/08/19
        if nblay==1 : dic['Layers']= [('0', int(2))] #EV 26/08/19
        dicIn['zolist']=[dic['Zones'][i][0] for i in range(len(dic['Zones'])) if dic['Zones'][i][1]==2]
        if ptyp!='X':
            lylist=[dic['Layers'][i][0] for i in range(len(dic['Layers'])) if dic['Layers'][i][1]==2]
            dicIn['lylist']=','.join(lylist)
        #dicIn['splist']=[self.res]
        dicIn['splist']=[self.rgroup.currentText()]
        if self.res=='Transport' :dicIn['splist']=['Tracer']
        if self.res=='Chemistry' :
            if ptyp=='X' or ptyp=='P': dicIn['plotOrder']='Zones'
            else :
                plotOrder=int(self.plgroup.currentIndex())
                if plotOrder==0 : dicIn['plotOrder']='Zones'
                if plotOrder==1 : dicIn['plotOrder']='Species'
            dicIn['splist']=[dic['Species'][i][0] for i in range(len(dic['Species'])) if dic['Species'][i][1]==2]
        else :dicIn['plotOrder']='Zones' 
        #print('dicIn',dicIn)
        return dicIn
    
    def buildPlot(self): #,dicIn
        '''this method build the piece of code that, when executed will make the graphs
        it has a dicIn as an input, which contains
        - type of plots (ptyp: 'time','space','XY')
        - lists of zones (zolist), species (splist), layers (lylist)
        - ordering by : plotOrder a list of two in : zones or species 
        '''
    ## Get the plot options
        dicIn = self.getOptions() #;print('dicIn',dicIn)
        self.ptyp,self.pOrder,self.zolist,self.splist,lylist =dicIn['ptyp'],dicIn['plotOrder'],dicIn['zolist'],dicIn['splist'],dicIn['lylist']
        if not dicIn['zolist'] : #EV 26/08/2019
            mess=onMessage(self.gui,'Choose zone(s) to plot the results.')
            return mess
        if not dicIn['splist'] : #EV 26/08/2019
            mess=onMessage(self.gui,'Choose specie(s) to plot the results.')
            return mess
        if not dicIn['lylist'] and self.ptyp[0]!='X': #EV 26/08/2019
            mess=onMessage(self.gui,'Choose layer(s) to plot the results.')
            return mess
    ## Calculates nb of plots, ncol and nrows for subplot
        if self.pOrder=='Zones': nplots = len(self.zolist) 
        if self.pOrder=='Species': nplots = len(self.splist)
        ncols = int(ceil(sqrt(nplots)))
        nrows = int(ceil(nplots/ncols))
    ## sets some plot type parameters
        #!if self.ptyp=='B0': ##time series
        if self.ptyp[0]=='B': ##time series
            iper = self.tlist ; self.axlabel='Time [T]' 
        #!if self.ptyp=='P0': ##Profile
        if self.ptyp[0]=='P': ##Profile
            curTime = int(self.Tstep.currentIndex())
            iper = curTime ; self.axlabel='Distance [L]' 
    ## sets some plot result parameters
        if 'Head' in self.splist : group = 'Flow' ; aylabel = 'Head [L]'
        elif 'Wcontent' in self.splist : group = 'Flow' ; aylabel = 'Wcontent' # OA 21/2/2019
        elif 'Tracer' in self.splist : group = 'Transport' ; aylabel = 'Concentration [NL$^{-3}$]'
        else : group = 'Chemistry' ; aylabel = 'Concentration [NL$^{-3}$]'
    ## build the plots
        self.figure.clf()
    ## Calibration graph
        if self.ptyp[0] == 'X':
            ptypXY='P0' ; self.llabel=[] ; time=[]
            self.yobs_all = [] ; self.ysim_all = [] ; self.obs_time = [] ; self.label_all=[] ## for export
            self._ax=self.figure.add_subplot(1,1,1)
            if int(self.Tstep.currentIndex())==0: curTime='All' ; time=self.tlist[1:] 
            else : curTime = int(self.Tstep.currentIndex()) ;time.append(self.tlist[curTime]);self.iper=curTime
            for i in range(len(self.splist)):
                for j in range(len(self.zolist)):
                    yobs_array = [] ; ysim_array = [] 
                    for t in range(len(time)):
                        if curTime=='All':self.iper=t
                        xobs,yobs,lobs=self.getDataObs(self.splist[i],self.zolist[j],time[t])
                        if len(yobs)!=0:
                            yobs_array.append(yobs[0])
                            x,y,label =  self.core.onPtObs(ptypXY,self.iper,group,self.zolist[j],[self.splist[i]],lobs)
                            ysim_array.append(float(y)) 
                            self.ysim_all.append(float(y))
                            self.yobs_all.append(yobs[0]) 
                            self.obs_time.append(time[t])
                            self.label_all.append(self.zolist[j]+'_'+self.splist[i]+'_lay'+str(lobs))                             
                    myplot=self._ax.plot(yobs_array,ysim_array,marker='x',linestyle = 'None')
                    self.llabel.append(self.zolist[j]+'_'+self.splist[i]+'_lay'+str(lobs)) 
            self.llabel.append('1:1')
            if not self.yobs_all : #EV 26/08/2019
                mess=onMessage(self.gui,'There is no observation data for these zone(s).') 
                return mess
            myplot2=self._ax.plot([min(self.yobs_all),max(self.yobs_all)],[min(self.yobs_all),max(self.yobs_all)],'k')
            self._ax.set_title('Simulated vs Observed Data: Time = '+str(curTime)+' [T]',fontweight="bold", size=9)
            self._ax.legend(self.llabel,fontsize = 8,loc='upper center', bbox_to_anchor=(0.5, -0.1),ncol=4)
            self._ax.set_ylabel('Simulated '+aylabel, fontsize = 8) 
            self._ax.set_xlabel('Observed '+aylabel, fontsize = 8)
            self._ax.ticklabel_format(useOffset=False, style='sci',scilimits=(-4,4),axis='both',useMathText=True)
            self._ax.tick_params(axis='both', labelsize=8)
            self._ax.figure.canvas.draw()
    ## Time series and profile graph
        else :
            a = lylist.split(',');layers = [int(x) for x in a] #EV 26/08/19
            if self.pOrder=='Zones': ## for Flow, transport and chemistry
                self.arryy=[] ; self.arrx=[] ## for export
                for i in range(nplots):
                    self._ax=self.figure.add_subplot(nrows,ncols,i+1)
                    self.x,yy,label =  self.core.onPtObs(self.ptyp,iper,group,self.zolist[i],self.splist,lylist) 
                    self.llabel=[label[i+1]+'(sim)' for i in range(len(label)-1)]
                    myplot=self._ax.plot(self.x,yy)
                    self.arryy.append(yy) ; self.arrx.append(self.x) ## for export
                    if self.ptyp[0]=='B': ## observed data for time series only
                        if self.checkBox.isChecked():
                            for j in range(len(self.splist)):
                                for l in range(len(layers)):
                                    xobs,yobs,lobs=self.getDataObs(self.splist[j],self.zolist[i],layers[l])
                                    if len(yobs)!=0 :
                                        color=j*len(layers)+layers.index(lobs) #EV 26/08/19
                                        myplot2=self._ax.plot(xobs,yobs,c='C'+str(color),marker='x',linestyle = 'None') #EV 26/08/19
                                        if group!='Chemistry':self.llabel.append('lay'+str(lobs)+'(obs)')
                                        else : self.llabel.append(str(self.splist[j])+'_lay'+str(lobs)+'(obs)')
                    self._ax.set_title(self.zolist[i],fontweight="bold", size=9)
                    self._ax.legend(self.llabel,fontsize = 8,loc='best')
                    self._ax.set_ylabel(aylabel, fontsize = 8) 
                    self._ax.set_xlabel(self.axlabel, fontsize = 8)
                    self._ax.ticklabel_format(useOffset=False, style='sci',scilimits=(-4,4),axis='both',useMathText=True)
                    self._ax.tick_params(axis='both', labelsize=8)
                    self._ax.figure.canvas.draw()
            if self.pOrder=='Species': ## for time series chemistry only
                self.arryy=[]
                for i in range(nplots):
                    self._ax=self.figure.add_subplot(nrows,ncols,i+1)
                    self.yyarray= np.empty((len(self.tlist),0)) ; self.llabel0=[] ; lobslab=[]
                    for j in range(len(self.zolist)):
                        self.x,yy,label =  self.core.onPtObs(self.ptyp,iper,group,self.zolist[j],[self.splist[i]],lylist) 
                        self.yyarray = np.append(self.yyarray,yy,axis=1)
                        label=[str(self.zolist[j]+'_lay'+lylist.split(',')[x]+'(sim)') for x in range(len(lylist.split(',')))]
                        self.llabel0.append(label)
                        if self.ptyp[0]=='B': ## observed data for time series only
                            if self.checkBox.isChecked():
                                for l in range(len(layers)):
                                    xobs,yobs,lobs=self.getDataObs(self.splist[i],self.zolist[j],layers[l])
                                    if len(yobs)!=0 : 
                                        lobslab.append(self.zolist[j]+'_lay'+str(lobs)+'(obs)')
                                        color=j*len(layers)+layers.index(lobs) #EV 26/08/19
                                        myplot2=self._ax.scatter(xobs,yobs,c='C'+str(color)) #EV 26/08/19
                    self.arryy.append(self.yyarray)
                    self.llabel0.append(lobslab)
                    self.llabel = [item for sublist in self.llabel0 for item in sublist]
                    myplot=self._ax.plot(self.x,self.yyarray)
                    self._ax.set_title(self.splist[i],fontweight="bold", size=9)
                    self._ax.legend(self.llabel,fontsize = 8,loc='best')
                    self._ax.set_ylabel(aylabel, fontsize = 8)
                    self._ax.set_xlabel(self.axlabel, fontsize = 8)
                    self._ax.ticklabel_format(useOffset=False, style='plain',axis='both')
                    self._ax.tick_params(axis='both', labelsize=8)
                    self._ax.figure.canvas.draw()
            
    def getDataObs(self,splist,zname,opt):
        '''read the data file and returns it
        splist: head, tracer or chemical species
        zname: name of the zone
        opt: time for '''
        xobs=[];yobs=[];lobs=0 
        if splist in ('Head','Tracer'):
            dicName=str('obs'+splist)
            dic = self.core.dicaddin[dicName]['data']
            ispe=3
        else :
            dicName='obsChemistry'
            dic = self.core.dicaddin[dicName]['data']
            ispe = self.core.dicaddin[dicName]['cols'].index(splist)
        for i in range(len(dic)):
            if dic[i][0]==zname:
                if self.ptyp[0]=='B':
                    if int(dic[i][1])==int(opt):
                        if dic[i][ispe]!='': # EV 05/03/19
                            xobs.append(float(dic[i][2]))     
                            yobs.append(float(dic[i][ispe]))
                            lobs=int(dic[i][1])
                else :
                    if float(dic[i][2])==float(opt):
                        if dic[i][ispe]!='': # EV 05/03/19
                            xobs.append(float(dic[i][2]))     
                            yobs.append(float(dic[i][ispe]))
                            lobs=int(dic[i][1])   
        return xobs,yobs,lobs
    
    def onExport(self):
        dlg = myFileDialog('Save')
        dlg = myFileDialog('Save')
        fDir,fName = dlg.getsetFile(self.gui,'Save','*.txt')
        if fDir == None: return
        f1 = open(fDir+os.sep+fName+'.txt','w')
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
        elif self.ptyp[0]=='P': 
            for i in range(len(self.zolist)):
                f1.write(self.axlabel.split(' ')[0])
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
