# -*- coding: utf-8 -*-

from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.Qt import QFrame
import qwt
from .config import *
from .geometry import *
from PyQt5.QtWidgets import *

def onMessage(gui,text):  QMessageBox.information(gui,"Info",text)
def onQuestion(gui,text):  
    q1 = QMessageBox.question(gui,"Info",text,QMessageBox.Yes,QMessageBox.No)
    if q1 == QMessageBox.Yes: return 'Yes'
    elif q1 == QMessageBox.No: return 'No'

class textDialog(QDialog):
    def __init__(self,gui, title, tsize, text):
        QDialog.__init__(self)
        self.setWindowTitle(title)
        self.glWidget = QWidget(self)
        self.glWidget.setGeometry(QRect(5, 5, tsize[0],tsize[1]))
        #scrollArea = QScrollArea(self.glWidget)
        self.gl = QVBoxLayout(self.glWidget)
        self.txed = QTextEdit(self.glWidget) #scrollArea)
        self.txed.setText(text)
        self.gl.addWidget(self.txed)
        #scrollArea.setWidget(self.txed)                
        self.buttonBox = QDialogButtonBox(self)
        self.buttonBox.setOrientation(Qt.Horizontal)
        self.buttonBox.setStandardButtons(QDialogButtonBox.Cancel|QDialogButtonBox.Ok)
        self.buttonBox.accepted.connect(self.accept)
        self.buttonBox.rejected.connect(self.reject)
        self.gl.addWidget(self.buttonBox)
        QMetaObject.connectSlotsByName(self)

    def showDialogAndDisconnect(self):
        self.show()
        self.gui.actionToggleEditing().triggered.disconnect(self.showDialogAndDisconnect)

    def getText(self):
        self.exec_()
        return str(self.txed.document().toPlainText())
        
class genericDialog(QDialog):
    def __init__(self, gui, title, data):
        self.gui,self.data = gui,data
        QDialog.__init__(self)
        self.setWindowTitle(title)
        self.glWidget = QWidget(self)
        y0=0
        for name,typ,value in self.data:
            if typ=='Textlong': y0+=1
        nb = len(self.data)
        self.screenShape = QDesktopWidget().screenGeometry()
        self.glWidget.setGeometry(QRect(5, 5, self.screenShape.width()*.15, nb*20+y0*160+30))
        #self.glWidget.setObjectName(_fromUtf8("gridLayoutWidget"))
        self.gl = QGridLayout(self.glWidget)
        self.gl.setMargin(2)
        self.gl.setSpacing(2)
        i=0;
        self.item = [0]*nb
        for name,typ,value in self.data:
            #y0=10+i*20
            txt = QLabel(self.glWidget)
            txt.setText(name) 
            self.gl.addWidget(txt,i,0,1,1)
            if typ == 'Choice':
                self.item[i] = QComboBox(self.glWidget)
                chlist = value[1];
                j=0
                for n in chlist:
                    self.item[i].addItem("")
                    self.item[i].setItemText(j,n)
                    j+=1
                self.item[i].setCurrentIndex(value[1].index(value[0]))
                self.gl.addWidget(self.item[i],i,1,1,1)
            elif typ=='Check':
                self.item[i] = QCheckBox(self.glWidget)
                self.item[i].setCheckState(value)
                self.gl.addWidget(self.item[i],i,1,1,1)
            elif typ=='Text':
                self.item[i] = QLineEdit(self.glWidget)
                self.item[i].setText(str(value))
                self.gl.addWidget(self.item[i],i,1,1,1)
            elif typ=='Textlong':
                scrollArea = QScrollArea(self.glWidget)
                #scrollArea.setGeometry(QRect(50, y0, 100, 50))
                scrollArea.setMaximumWidth(self.screenShape.width()*.05)
                scrollArea.setFixedHeight(160)
                self.item[i] = QTextEdit(scrollArea)
                #y0 += 30
                if type(value)==type([5,6]): 
                    value = '\n'.join([str(v) for v in value]) # cretaes a text with line returns
                self.item[i].setText(str(value))
                scrollArea.setWidget(self.item[i])                
                self.gl.addWidget(scrollArea,i,1,1,1)
            i+=1
        self.buttonBox = QDialogButtonBox(self)
        self.buttonBox.setOrientation(Qt.Horizontal)
        self.buttonBox.setStandardButtons(QDialogButtonBox.Cancel|QDialogButtonBox.Ok)

        self.glWidget2 = QWidget(self)
        self.glWidget2.setGeometry(QRect(5, nb*20+y0*160+30, self.screenShape.width()*.15, 40))
        self.gl2 = QGridLayout(self.glWidget2)
        self.gl2.setMargin(0)
        self.gl2.addWidget(self.buttonBox) #,nb,1,2,1)
        self.buttonBox.accepted.connect(self.accept1)
        self.buttonBox.rejected.connect(self.reject1)
        QMetaObject.connectSlotsByName(self)

    def showDialogAndDisconnect(self):
        self.show()
        self.gui.actionToggleEditing().triggered.disconnect(self.showDialogAndDisconnect)

    def accept1(self): 
        self.close(); self.state = 'accept'
    def reject1(self): 
        self.close(); self.state = 'reject'

    def getValues(self):
        self.exec_()
        nb = len(self.data)
        val = [0]*nb
        for i in range(nb):
            typ = self.data[i][1]
            if typ == 'Choice': val[i] = str(self.item[i].currentText())
            if typ == 'Text': val[i] = str(self.item[i].text())
            if typ == 'Check': val[i] = self.item[i].checkState()
            if typ == 'Textlong': 
                v0 = str(self.item[i].document().toPlainText())
                val[i] = v0.split('\n')            
        if self.state =='accept' : return val
        else : return None
        
class myFileDialog:
    def __init__(self,opt='Open'):
        self.opt = opt
    def getsetFile(self,gui,title,filt):
        dlg = QFileDialog()
        if self.opt=='Open':
            fileName = dlg.getOpenFileName(gui,title,' ',filter=filt)
        elif self.opt in ['New','Save']:
            dlg.setFileMode(QFileDialog.AnyFile)
            fileName = dlg.getSaveFileName(gui,title,' ',filter=filt)
        fileName = fileName[0] # OA 1/10 for python 3
        fName = fileName.split('/')[-1]
        fDir = fileName.replace(fName,'')
        fName = fName.split('.')[0]
        #onMessage(gui,str(fDir)+' '+str(fName))
        #settings = QSettings()
        #settings.setValue("last_file", fDir)# QVariant(QString(fDir)))
        return str(fDir),str(fName)

class myNoteBookCheck(QDialog):
    def __init__(self, gui,title, dicIn):
        QDialog.__init__(self)
        self.setWindowTitle(title)
        self.gui,self.pages,self.layouts,self.dicIn = gui,{},{},dicIn
        self.dicOut = self.dicIn.copy()
        layout = QVBoxLayout(self)
        self.setGeometry(QRect(40, 40, 480,550))
        glWidget = QWidget(self)
        nb = QTabWidget(glWidget)
        nb.setGeometry(QRect(5, 5, 450,520))
        self.dwidget = {}
        for n in list(dicIn.keys()):
            if dicIn[n]==None:continue
            if len(dicIn[n])==0: continue
            nbChk = len(dicIn[n])
            pg = QWidget(nb)
            lay = QGridLayout(pg)
            self.layouts[n]= lay
            lay.setMargin(0);lay.setSpacing(0)
            self.dwidget[n] = [0]*nbChk
            for i in range(nbChk):
                ic = mod(i,4);il = i/4
                ch = QCheckBox(nb); self.dwidget[n][i] = ch
                ch.setText(dicIn[n][i][0])
                ch.setCheckState(dicIn[n][i][1])
                lay.addWidget(ch,il,ic)
            scroll = QScrollArea()
            scroll.setWidget(pg)
            scroll.setWidgetResizable(True)
            scroll.setFixedHeight(400)
            self.pages[n] = pg
            nb.addTab(scroll,n)
            pg.show()
        layout.addWidget(glWidget)
        buttonBox = QDialogButtonBox(self)
        buttonBox.setStandardButtons(QDialogButtonBox.Cancel|QDialogButtonBox.Ok)
        buttonBox.accepted.connect(self.accept1)
        buttonBox.rejected.connect(self.reject1)
        layout.addWidget(buttonBox)
        QMetaObject.connectSlotsByName(self)
        
    def accept1(self): 
        self.close(); self.state = 'accept'
    def reject1(self): 
        self.close(); self.state = 'reject'
        
    def showDialogAndDisconnect(self):
        self.show()
        self.gui.actionToggleEditing().triggered.disconnect(self.showDialogAndDisconnect)
            
    def getValues(self):
        self.exec_()
        for k in list(self.dicIn.keys()):
            if k in list(self.pages.keys()):
                names,boo = list(zip(*self.dicIn[k]))
                lout = []
                items = self.dwidget[k]
                for item in items:
                    lout.append(item.checkState())
                self.dicOut[k] = list(zip(names,lout))
        if self.state == 'accept': return self.dicOut
        else : return None

class myNoteBook(QDialog):
    def __init__(self, gui,title, dicIn):
        QDialog.__init__(self)
        self.setWindowTitle(title)
        self.gui,self.pages,self.dicIn = gui,{},dicIn
        self.dicOut = self.dicIn.copy()
        layout = QVBoxLayout(self)
        self.screenShape = QDesktopWidget().screenGeometry()
        self.setGeometry(QRect(40, 60, self.screenShape.width()*.42,self.screenShape.height()*.6))
        glWidget = QWidget(self)
        nb = QTabWidget(glWidget)
        nb.setGeometry(QRect(5, 20, self.screenShape.width()*.4,self.screenShape.height()*.50))
        for n in list(dicIn.keys()):
            if dicIn[n]==None:continue
            if len(dicIn[n]['rows'])==0 and n!='Species': continue
            pg = myNBpanelGrid(gui,nb,dicIn[n]) #,size=(450,500))
            self.pages[n] = pg
            nb.addTab(pg,n)
            pg.show()
        layout.addWidget(glWidget)
        buttonBox = QDialogButtonBox(self)
        buttonBox.setStandardButtons(QDialogButtonBox.Cancel|QDialogButtonBox.Ok)
        buttonBox.accepted.connect(self.accept)
        buttonBox.rejected.connect(self.reject)
        layout.addWidget(buttonBox)
        QMetaObject.connectSlotsByName(self)
        
    def showDialogAndDisconnect(self):
        self.show()
        self.gui.actionToggleEditing().triggered.disconnect(self.showDialogAndDisconnect)
            
    def getValues(self):
        self.exec_()
        for k in list(self.dicIn.keys()):
            if k in list(self.pages.keys()):
                #onMessage(self.gui,k)
                self.dicOut[k] = self.pages[k].getValues()
        return self.dicOut
            
class myNBpanelGrid(QTableWidget):       
    def __init__(self,gui,parentWidget,dicIn,size=(60,40)):
        QTableWidget.__init__(self,5,3)
        self.screenShape = QDesktopWidget().screenGeometry()
        rows,cols,data = dicIn['rows'],dicIn['cols'],dicIn['data'];#onMessage(gui,str(data))
        self.dicIn,self.dicOut = dicIn, dicIn.copy()
        self.setColumnCount(len(cols))
        self.setRowCount(len(rows))
        self.setHorizontalHeaderLabels(cols)
        self.setVerticalHeaderLabels(rows)
        self.setFont(QFont('Arial',pointSize=9))
        self.type = ['Text']*len(cols)
        for il,line in enumerate(data):
            self.setRowHeight(il,self.screenShape.height()*.03)
            for ic,item in enumerate(line): #print il,ic,item
                #onMessage(gui,str(item)+str(type(item)))
                if type(item) == type(True):
                    if il==0: self.type[ic] = 'Check'
                    a = QCheckBox()
                    a.setCheckState(item)
                    self.setCellWidget(il,ic, a)
                elif type(item)==type(5.6) and len(str(item))>8:
                    self.setItem(il,ic, QTableWidgetItem(str(item)[:7]))
                else  :
                    self.setItem(il,ic, QTableWidgetItem(str(item)))
        self.resizeColumnsToContents()
        self.resize(size[0],size[1])
        QMetaObject.connectSlotsByName(self)
        self.clip = QApplication.clipboard()
        #cpCells = CopySelectedCellsAction(gui,self)
        
    def getValues(self):
        #self.exec_()
        dicOut = {'data':[],'cols':self.dicIn['cols'],'rows':self.dicIn['rows']}
        for il in range(self.rowCount()):
            #if self.item(il,0)==None: continue
            l0=[]
            for ic in range(self.columnCount()):
                #print 'qtdl 269',il,ic,self.item(il,ic),self.type[ic]
                if self.type[ic]=='Text':
                    l0.append(str(self.item(il,ic).text()))
                else :
                    cell = self.cellWidget(il,ic)
                    l0.append(cell.checkState()!=0)
            dicOut['data'].append(l0)
        return dicOut

    def keyPressEvent(self, e):
        # from http://retrofocus28.blogspot.com/2013/09/pyqt-qtablewidget-copy-and-past-excel.html
        if (e.modifiers() & Qt.ControlModifier):
            selected = self.selectedRanges()
            if e.key() == Qt.Key_V:#past
                first_row = selected[0].topRow()
                first_col = selected[0].leftColumn()
                #copied text is split by '\n' and '\t' to paste to the cells
                rowText = self.clip.text().split('\n')
                self.setRowCount(len(rowText)-1)
                for r, row in enumerate(rowText):
                    for c, text in enumerate(row.split('\t')):
                        self.setItem(first_row+r, first_col+c, QTableWidgetItem(text))

            elif e.key() == Qt.Key_C: #copy
                s = ""
                for r in range(selected[0].topRow(),selected[0].bottomRow()+1):
                    for c in range(selected[0].leftColumn(),selected[0].rightColumn()+1):
                        try:
                            s += str(self.item(r,c).text()) + "\t"
                        except AttributeError:
                            s += "\t"
                    s = s[:-1] + "\n" #eliminate last '\t'
                self.clip.setText(s)
            
'''##########################" FOR ZONES  ###############################""
'''

class zoneDialog(QDialog):
    def __init__(self, parent, core,model,line, curzones, nb):
        QDialog.__init__(self)
        self.gui, self.core, self.model, self.line,  = parent.gui, core,model,line
        self.nb = nb
        #self.test = 'erty'
        self.glWidget = QWidget(self)
        self.screenShape = QDesktopWidget().screenGeometry()
        self.glWidget.setGeometry(QRect(5, 5, self.screenShape.width()*0.4,self.screenShape.height()*.5))
        self.gl = QGridLayout(self.glWidget)
        if 'names' in list(core.dickword[model].lines[line].keys()): self.typO=1 # several values
        else : self.typO = 0
        self.zpanel = zonePanel(self, curzones.copy(), nb, self.typO)
        self.gl.addWidget(self.zpanel)
        self.buttonBox = QDialogButtonBox(self)
        self.buttonBox.setOrientation(Qt.Horizontal)
        self.buttonBox.setStandardButtons(QDialogButtonBox.Cancel|QDialogButtonBox.Ok)
        self.buttonBox.accepted.connect(self.accept1)
        self.buttonBox.rejected.connect(self.reject1)
        self.gl.addWidget(self.buttonBox)
        #onMessage(self.gui,'base zdialg \n'+str(self)+'\n'+str(self.zones))

    def showDialogAndDisconnect(self):
        self.show()
        self.gui.actionToggleEditing().triggered.disconnect(self.showDialogAndDisconnect)
        
    def saveCurrent(self):
        # saving the entered zone
        self.exec_()
        zp = self.zpanel;
        curdic = self.core.diczone[self.model]
        #curdic.addZone(self.line)
        curzones = curdic.dic[self.line]
        #onMessage(self.gui,'in save current \n'+str(curzones))
        curzones['name'][self.nb] = str(zp.name.text())
        media = zp.media.text()
        if '-' in media :
            m1 = media.split('-')
            m2 = list(range(int(m1[0]),int(m1[1])+1))
        else :
            m2 = int(media)
        curzones['media'][self.nb] = m2;#onMessage(self.gui,str(zp.coords.getValues()['data']))
        curzones['coords'][self.nb]= self.corrCoords(zp.coords.getValues()['data'])
        val0 = ''
        if self.typO: val0 = self.getOpt()
        #print 'qtdlg 366',self.typS
        if self.typS==1: # for transient values
            v0 = zp.valBox.getValues()['data']
            val = ''
            for b in v0 : val+=str(b[0])+' '+str(b[1])+'\n'
        else :
            val = str(zp.valBox.document().toPlainText())            
        first = val.split('\n')[0]
        if first.count('.')>2:   #for pht3d zones 1.0.0.0
            val = val.replace('.','')
            val = val.replace(' ','')
            val = val.replace('\n','')
        curzones['value'][self.nb]= val0+val
        if self.state == 'accept': return 'OK'
        else : return None
        
    def accept1(self): 
        self.close()
        self.state = 'accept'
    def reject1(self): 
        self.close()
        self.state = 'reject'
           
    def corrCoords(self,lcooI):
        '''change coordinates if they are out of the domain'''
        g = self.core.dicaddin['Grid']
        ex = (float(g['x1'])-float(g['x0']))/1e5
        ey = (float(g['y1'])-float(g['y0']))/1e5
        lcoord = [];#onMessage(self.gui,str(lcooI))
        for b in lcooI:
            try : x = clip(float(b[0]),float(g['x0'])+ex,float(g['x1'])-ex)
            except ValueError : continue
            try: y = clip(float(b[1]),float(g['y0'])+ey,float(g['y1'])-ey)
            except ValueError : continue
            if len(b)==2: lcoord.append((x,y))
            else : lcoord.append((x,y,float(b[2])))
        return lcoord
        
    def getOpt(self):
        data = self.zpanel.GOdata.getValues()['data']; #print 'wxdialog 248',data
        s = '$'
        for i in range(len(data)):
            s += data[i][0]+' \n'
        return s+'$'
     
class zonePanel(QWidget):
    def __init__(self,parent, curzones, nb, typO):
        self.parent,self.core,self.model,self.line = parent,parent.core,parent.model,parent.line
        QWidget.__init__(self)
        self.mLayout = QGridLayout(self)
        self.screenShape = QDesktopWidget().screenGeometry()
        ################# zone name
        txtName = QLabel('  Zone name  ')
        self.name = QLineEdit(self)
        self.name.setText(curzones['name'][nb])
        self.mLayout.addWidget(txtName,0,0)
        self.mLayout.addWidget(self.name,0,1)
        #################" zone media
        txtMedia = QLabel('  Zone media  ')
        zm = '0'
        if curzones['media'][nb]!='': 
            zm = curzones['media'][nb]
        if type(zm)==type([5]): zm=str(zm[0])+'-'+str(zm[-1])
        self.media = QLineEdit(self)
        self.media.setText(str(zm))
        self.mLayout.addWidget(txtMedia,1,0)
        self.mLayout.addWidget(self.media,1,1)
        #if parent.core.addin.getDim() !='3D': self.media.Enable(False)
        ############# zone type : transient or not
        txtTransient = QLabel('  Temporal  ')
        self.transient = QComboBox(self)
        self.transient.addItems(['Steady','Transient'])
        self.transient.activated['QString'].connect(self.onTransient)
        self.mLayout.addWidget(txtTransient,2,0)
        self.mLayout.addWidget(self.transient,2,1)
        #############   zone coords
        txtCoords = QLabel('  Zone coords ')
        coo = curzones['coords'][nb]
        data = [list(x) for x in coo];
        nrow = len(data)
        cols = [' X ',' Y ']
        if len(data[0])==3: cols = ['   X   ','  Y   ','Z']
        dicCoo = {'cols':cols,'rows':['']*nrow,'data':data}
        self.coords = myNBpanelGrid(self.core.gui,self,dicCoo,size=(200,self.screenShape.height()*.3))
        self.mLayout.addWidget(txtCoords,3,0)
        self.mLayout.addWidget(self.coords,3,1)
        ##################" zone value
        s = '       '
        val = str(curzones['value'][nb]);
        if '$' in val : val = val.split('$')[-1]
        self.value = val
        typS = 0 # classical
        if len(val.split('\n'))>1: typS = 1 # transient
        parent.typS = typS
        #print 'qtdlg 443', typS
        if parent.line in ['ph.3','inic.1']: # for chemical solutions
            a = str(val).rjust(5,'0')
            val =  str(int(a[:2]))+'.'+a[-3]+'.'+a[-2]+'.'+a[-1]
            units = QLabel(s+'solution . mineral . exchange . surface')
        ################# options linename, units, several values
        lineName = QLabel(s+parent.core.dickword[self.model].lines[self.line]['comm'])
        units = QLabel(s+parent.core.getUnits(self.model,self.line,0))
        self.mLayout.addWidget(lineName,0,2,1,2)
        self.mLayout.addWidget(units,1,2)
        if typO:
            self.GOtitles,self.GOdata = self.gridOpt(self.line,str(curzones['value'][nb]))
            self.mLayout.addWidget(self.GOdata,2,3)
        txtValue = QLabel('  Zone value  ')
        self.mLayout.addWidget(txtValue,3,2)
        self.valBox = self.changeValueBox(typS)
        self.mLayout.addWidget(self.valBox,3,3)
        self.transient.setCurrentIndex(typS)
        self.setLayout(self.mLayout)
        #self.curzones=curzones
        #onMessage(self.parent.gui,str(parent.curzones))
        QMetaObject.connectSlotsByName(self)
        
    def onTransient(self,evt):
        self.valBox.deleteLater()  # OA 10/9/18
        self.valBox = self.changeValueBox(self.transient.currentIndex())
        self.mLayout.addWidget(self.valBox,3,3) # OA 10/9/18
        self.parent.typS = 1 - self.parent.typS  # OA 10/9/18
        
    def changeValueBox(self,typS):
        '''changes the value sizer btw one value (permanent) and a list (transient)'''
        #2for item in self.valSizer.children(): self.valSizer.removeItem(item)
        if typS == 0: # text
            valBox = QTextEdit(self) #size=size3)
            valBox.setText(self.value)
        elif typS == 1: # grid
            if self.value[0]=='$': a,b,val1 = self.value.split('$')
            else : 
                val1 = self.value.split('\n');
                while '' in val1: val1.remove('')
                #print val1
            if len(val1)>1: data = [v.split() for v in val1]
            else : data = [[i,0] for i in range(3)]# no values
            nr = len(data)
            dicV = {'cols':['t','val'],'rows':['']*nr,'data':data};#print 'qtdialg, 483',dicV
            valBox = myNBpanelGrid(self.core.gui,self,dicV,size=(120,180))

        #valBox.cellClicked.connect(self.current)
        #self.valSizer.addWidget(valBox)
        return valBox
        
    def gridOpt(self,line,val):
        titles = self.core.dickword[self.model].lines[line]['names']
        if line == 'mnwt.2a': 
            return titles,self.gridMNW(titles, val)
        else: 
            return titles,self.gridSimple(titles, line,val)

    def gridSimple(self,titles,line,val):
        '''creates specific dialog part for serveal values in a zone'''
        if '$' in val:
            a,P1,P2 = val.split('$')
            Parm1 = [[x] for x in P1.split('\n')];#print 'Nwell',Parm1
        else : 
            P1 = self.core.dickword[self.model].lines[line]['default']
            Parm1 = [[str(x)] for x in P1]
        dicV = {'cols':[' '*20],'rows':titles,'data':Parm1}
        self.gopt = myNBpanelGrid(self.core.gui,self,dicV,size=(180,20*len(titles)))
        return self.gopt
    '''
    def gridMNW(self,titles,val):
        #creates specific dialog parts for multi node wells
        #now no implementation of Qlimit, pumploc, ppflag, pumpcap
        #print val
        #Parm1 = [[0]*len(t.split(',')) for t in titles]
        if '$' in val:
            a,P1,P2 = val.split('$')
            Parm1 = [[x] for x in P1.split('\n')];#print 'Nwell',Parm1
        else : 
            Parm1 = [['0'] for x in titles]
        lossV = Parm1[1][0] # store the value of lossV
        Parm1[1][0]=['NONE','THIEM','SKIN','GENERAL','SPECIFYcwc']
        dicV = {'cols':[' '*20],'rows':titles,'data':Parm1}
        self.gopt = myNBpanelGrid(self.parent.gui,self,dicV,size=(180,180))
        self.gopt.setChoice(1,0,lossV)
        if '$' in val:
            self.changeCellR(0)
            self.changeCellR(1) # to see the effect of the current selection
        self.Bind(wgrid.EVT_GRID_CELL_CHANGE,self.changeCell)
        return titles,self.gopt
        
    def changeCell(self,evt):
        ir = evt.GetRow();#self.curRow = ir
        self.changeCellR(ir)
        
    def changeCellR(self,ir):
        """different data input for type of MN well"""
        if self.gopt.GetRowLabelValue(ir)=='Losstype': # Loss Type
            for i in range(2,6): self.gopt.SetRowLabelValue(i,' ')
            Ltype = self.gopt.GetCellValue(ir,0).strip()
            dlist = {'NONE':[],'THIEM':['Rw'],'SKIN':['Rw','Rskin','Kskin'],\
                'GENERAL':['Rw','B','C','P'],'SPECIFYcwc':['CWC']}
            #print tlist.keys()
            r = dlist[str(Ltype)]
            for i in range(len(r)):
                self.gopt.SetRowLabelValue(2+i,r[i])
        elif self.gopt.GetRowLabelValue(ir)=='Nnodes': # nnodes
            for i in range(6,8): self.gopt.SetRowLabelValue(i,' ')
            Nnodes = int(self.gopt.GetCellValue(ir,0))
            if Nnodes >0 : self.gopt.SetRowLabelValue(6,'LAY')
            elif Nnodes <0 : 
                self.gopt.SetRowLabelValue(6,'Ztop')
                self.gopt.SetRowLabelValue(7,'Zbot')
        self.gopt.Refresh()
 '''       
    def current(self,row,col):
        #print  evt.GetRow()
        if row ==self.valBox.rowCount()-1:
            self.valBox.insertRow(row+1);self.valBox.rowString.append('')
        #evt.Skip()

#//////////////////////////////////////////////////////////////////////
class dialogContour(QDialog):

    def __init__(self, parent, title, valeur, col):
        """ liste contient les attributs actuels des contours : val : [0]min, 1:max,
        2: intervalles, [3]decimales, 4:log, 5:user puis couleurs et transparence"""
        self.listCtrl,self.valeur,self.parent = [], valeur,parent;
        QDialog.__init__(self)
        self.setWindowTitle("Contour")
        self.glWidget = QWidget(self)
        self.glWidget.setGeometry(QRect(5, 5, 250, 10*35+50))
        self.gl = QGridLayout(self.glWidget)
        self.gl.setMargin(2)
        self.gl.setSpacing(2)
        # boite pour calcul automatique
        txt = QLabel(self.glWidget);txt.setText('Automatic')
        self.gl.addWidget(txt,0,0,1,1)
        self.auto = QCheckBox(self.glWidget);self.auto.setCheckState(True)
        self.gl.addWidget(self.auto,0,1,1,1)
        if valeur!=None:
            self.auto.setCheckState(valeur[4]=='auto')
            self.listuser=valeur[5]
        # intervals
        label = ['Mini', 'Maxi', 'Interval','Decimals']
        if valeur==None: valeur=[0.,10.,1.,2,False,False]
        for i in range(4):
            txt = QLabel(self.glWidget);txt.setText(label[i])
            self.gl.addWidget(txt,i+1,0,1,1)
            Vl = QLineEdit(self.glWidget);Vl.setText(str(valeur[i]))
            self.listCtrl.append(Vl)
            self.gl.addWidget(Vl,i+1,1,1,1)
            
        txt = QLabel(self.glWidget);txt.setText('log')
        self.gl.addWidget(txt,5,0,1,1)
        self.log = QCheckBox(self.glWidget);self.log.setCheckState(valeur[4]=='log')
        self.gl.addWidget(self.log,5,1,1,1)

        self.butlist = QPushButton(self.glWidget)
        self.butlist.clicked.connect(self.onListUser)
        self.gl.addWidget(self.butlist,6,0,1,1)
        self.user = QCheckBox(self.glWidget)
        self.user.setCheckState(valeur[4]=='fix')
        self.gl.addWidget(self.user,6,1,1,1)
        
        self.buttonBox = QDialogButtonBox(self)
        self.buttonBox.setOrientation(Qt.Horizontal)
        self.buttonBox.setStandardButtons(QDialogButtonBox.Cancel|QDialogButtonBox.Ok)
        self.buttonBox.accepted.connect(self.accept1)
        self.buttonBox.rejected.connect(self.reject1)
        self.gl.addWidget(self.buttonBox)

        QMetaObject.connectSlotsByName(self)

    def accept1(self): self.close(); self.state = 'accept'
    def reject1(self): self.close(); self.state = 'reject'

    def GetStrings(self):
        """renvoie les valeurs des boites et ajoute la liste user a la fin """
        self.exec_()
        v=self.valeur;#print v
        if v==None: v=[0.,10.,1.,2,'auto',None]
        for i in range(4):
            v[i]=self.listCtrl[i].text()
            try: v[i] = float(v[i])
            except ValueError:
                self.parent.OnMessage('erreur de type');return self.valeur
        v[4]='lin'
        if self.user.checkState():
            v[4]='fix';v[5]=self.listuser
        if self.log.checkState(): v[4]='log'
        if self.auto.checkState(): v[4]='auto'
        #print v
        return v
    
    def onListUser(self,event):
        """ opens a dialog to set some values"""
        if self.user.getValue()==False: return
        data  =[('values','Textlong',[''])]
        if self.valeur!=None:
            if type(self.valeur[5])==type([1]): 
                data =[('values','Textlong',self.valeur[5])]
        #dialg = MyListDialog(self.parent,"",lst1)
        dialg = genericDialog(self.parent,'List of values',data)
        ls = dialg.getValues()[0];#print 'lstuser',self.listuser
        self.listuser = [float(x) for x in ls]
        return
''' ##################"" Instant FIT DIALOG (added in wx 9/4/2017 oa then 20/2/18 for qt) #############""""
'''

class instantFitDialog(QDialog):
    """This dialog allows to choose the type of calculation (head, transport),
    to change the value of the current zone with button
    fix the value of dispersivity for transport
    parent is the instant class in addin
    """
    
    def __init__(self, gui, parent, dic_options):

        self.gui,self.parent,self.dic_options = gui,parent,dic_options
        QDialog.__init__(self,self.gui)
        self.setModal(0)
        self.setWindowTitle("Instant fit")
        self.glWidget = QWidget(self)
        self.gl = QGridLayout(self.glWidget)
        self.gl.setMargin(2)
        self.gl.setSpacing(2)
        txt = QLabel(self.glWidget);txt.setText('type') 
        self.gl.addWidget(txt,0,0,1,1)
        self.chT = QComboBox(self.glWidget)
        lst = ['Head','Tracer']
        self.chT.addItems(lst)
        self.chT.setCurrentIndex(lst.index(dic_options['type']))
        self.gl.addWidget(self.chT,0,1,1,1)
        self.chT.activated['QString'].connect(self.onChoiceT)
        # long dispersivity
        txt = QLabel(self.glWidget);txt.setText('Long. Disp') 
        self.gl.addWidget(txt,1,0,1,1)
        self.aL = QLineEdit(self.glWidget)
        self.aL.returnPressed.connect(self.onAL)
        self.gl.addWidget(self.aL,1,1,1,1)
        # lateral dispersivity
        txt = QLabel(self.glWidget);txt.setText('Lat. Disp') 
        self.gl.addWidget(txt,2,0,1,1)
        self.aT = QLineEdit(self.glWidget)
        self.aT.returnPressed.connect(self.onAT)
        self.gl.addWidget(self.aT,2,1,1,1)
        QMetaObject.connectSlotsByName(self)

    def onChoiceT(self): 
        self.dic_options['type'] = self.chT.currentText()
        self.change()
    def onAL(self): 
        self.dic_options['aL'] = self.aL.text()
        self.change()
    def onAT(self): 
        self.dic_options['aT'] = self.aT.text()
        self.change()
                
    def change(self):
        '''look a the change, sets dic_options and sends it back to parent'''
        self.parent.dic_options = self.dic_options
        self.parent.update('dialog')
        
    def OnExit(self,evt): 
        self.parent.end()
        self.destroy()

#/////////////////////////////////////////:  plot XY    ////////////////::
class plotxy(QDialog):
    
    def __init__(self,gui,x,arry,legd,title,Xtitle,Ytitle,typ='-'):
        QDialog.__init__(self)
        layout = QVBoxLayout(self)
        self.gui,self.data,self.title,self.legd,self.Xtitle,self.x,self.arry = gui,3,title,legd,Xtitle,x,arry
        lab = QLabel(title)
        layout.addWidget(lab)
        glWidget = QWidget(self)
        glWidget.setGeometry(QRect(0, 0, 200, 300))
        self.cnv = qwt.QwtPlot(glWidget);
        x,arry = transpose(array(x,ndmin=2)),array(arry) # x en vertical
        self.lignes,cols = [],['red','blue','green','orange','cyan','black']*5
        # verify size of input vectors
        if len(shape(x))==1:
            x1 = ones((len(x),1))+0.;
        if len(shape(arry))==1:
            arry1 = ones((len(arry),1))+0.; arry1[:,0]=arry; arry=arry1*1.;ny=1
        else :
            [nt,ny] = shape(arry)
        x2 = x*1.; arry2 = arry*1.;
        # creer les lignes
        colors = [Qt.black,Qt.red,Qt.blue,Qt.green,Qt.cyan,Qt.magenta,Qt.darkGreen,Qt.darkRed]
        if ny==1:
            if typ=='-': 
                gobj = qwt.QwtPlotCurve(legd[0])
                gobj.setData(x2,arry2)
            else : #n xy plot
                gobj = qwt.QwtPlotCurve() #legd[1])]
            self.lignes.append(gobj)
        else :
            for i in range(ny):
                if typ=='-': 
                    gobj = qwt.QwtPlotCurve(legd[i])
                    gobj.setPen(QPen(colors[mod(i,8)]))
                    gobj.setData(x2,arry2[:,i])
                else : 
                    gobj = qwt.QwtPlotCurve(legd[i])
                    gobj.setPen(QPen(colors[mod(i,8)]))
                self.lignes.append(gobj)
        self.cnv.insertLegend(qwt.QwtLegend(), Qwt.QwtPlot.RightLegend)
        self.cnv.setAxisTitle(qwt.QwtPlot.xBottom, Xtitle)
        self.cnv.setAxisTitle(qwt.QwtPlot.yLeft, Ytitle)
        # insert curves
        for ligne in self.lignes : ligne.attach(self.cnv)
        layout.addWidget(self.cnv)
        self.cnv.replot()
            
        basWidget = QWidget(self)
        basLayout = QHBoxLayout(basWidget);
        bexport = QPushButton(basWidget);bexport.setText('Export')
        basLayout.addWidget(bexport)
        bexport.clicked.connect(self.onExport)
        bcancel = QPushButton(basWidget);bcancel.setText('Cancel')
        basLayout.addWidget(bcancel)
        bcancel.clicked.connect(self.reject)
        layout.addWidget(basWidget) #,nb,1,2,1)
        
        QMetaObject.connectSlotsByName(self)

    def showDialogAndDisconnect(self):
        self.show()
        self.gui.actionToggleEditing().triggered.disconnect(self.showDialogAndDisconnect)
        
    def getValues(self):
        self.exec_()
        return self.data
    
    def onExport(self,evt):
        #print self.lignes
        #core = self.gui.core
        #fName = core.fileDir+os.sep+core.fileName+self.title+'.txt'
        dlg = myFileDialog('Save')
        fDir,fName = dlg.getsetFile(self.gui,'Save','*.txt')
        if fDir == None: return
        f1 = open(fDir+os.sep+fName+'.txt','w')
        f1.write(self.Xtitle)
        for n in self.legd: f1.write(' '+n)
        f1.write('\n')
        nt,ny = shape(self.arry)
        arr = zeros((nt,ny+1))
        arr[:,0]=self.x;arr[:,1:]=self.arry
        savetxt(f1,arr)
        f1.close()
