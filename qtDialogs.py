# -*- coding: utf-8 -*-

from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from PyQt5.Qt import QFrame
#import qwt
from .config import *
from .geometry import *

from functools import partial
#from scipy import * #OA 2/4/19

def onMessage(gui,text):  QMessageBox.information(gui,"Info",text)
def onQuestion(gui,text):  
    q1 = QMessageBox.question(gui,"Info",text,QMessageBox.Yes,QMessageBox.No)
    if q1 == QMessageBox.Yes: return 'Yes'
    elif q1 == QMessageBox.No: return 'No'

class textDialog(QDialog): # Text dialog for batch, formula, initial chemistry
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
        
class genericDialog(QDialog): # Dialog for addin parameters and options for plot results
    def __init__(self, gui, title, data):
        self.gui,self.data,self.color = gui,data,data[0][2]  # OA 22/11/19 added color
        QDialog.__init__(self)
        self.setWindowTitle(title)
        self.glWidget = QWidget(self)
        y0=0
        for name,typ,value in self.data:
            if typ=='Textlong': y0+=1
        nb = len(self.data)
        self.screenShape = QDesktopWidget().screenGeometry()
        self.glWidget.setGeometry(QRect(5, 5, int(self.screenShape.width()*.15), nb*20+y0*160+30))
        #self.glWidget.setObjectName(_fromUtf8("gridLayoutWidget"))
        self.gl = QGridLayout(self.glWidget)
        self.gl.setContentsMargins(1,1,1,1) # OA 3/10/18
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
                for j,n in enumerate(chlist):
                    self.item[i].addItem("")
                    self.item[i].setItemText(j,n)
                self.item[i].setCurrentIndex(value[1].index(value[0]))
                self.gl.addWidget(self.item[i],i,1,1,1)
            elif typ=='Check':
                self.item[i] = QCheckBox(self.glWidget) #EV 18/02/20
                self.item[i].setChecked(value)
                self.gl.addWidget(self.item[i],i,1,1,1)   
            elif typ=='CheckList':
                self.item[i] = CheckableComboBox()
                chlist = value[1];
                for j,n in enumerate(chlist):
                    self.item[i].addItem(n)
                    self.item[i].setItemChecked(j, False)
                self.gl.addWidget(self.item[i],i,1,1,1)   
            elif typ=='Text':
                self.item[i] = QLineEdit(self.glWidget)
                self.item[i].setText(str(value))
                self.gl.addWidget(self.item[i],i,1,1,1)
            elif typ=='Textlong':
                scrollArea = QScrollArea(self.glWidget)
                #scrollArea.setGeometry(QRect(50, y0, 100, 50))
                scrollArea.setMaximumWidth(int(self.screenShape.width()*.08)) #EV 05/08/19 0.05 -> 0.08
                scrollArea.setFixedHeight(160)
                self.item[i] = QTextEdit(scrollArea)
                #y0 += 30
                if type(value)==type([5,6]): 
                    value = '\n'.join([str(v) for v in value]) # cretaes a text with line returns
                self.item[i].setText(str(value))
                scrollArea.setWidget(self.item[i])                
                self.gl.addWidget(scrollArea,i,1,1,1)
            elif typ=='Color': # OA added 22/11/19
                self.item[i] = QPushButton('Color',self.glWidget)
                self.item[i].clicked.connect(self.onColor)
                self.gl.addWidget(self.item[i],i,1,1,1)
            elif typ=='File': #EV 07/02/20
                self.item[i] = QPushButton('Browse',self.glWidget)
                self.item[i].clicked.connect(self.onFile)
                self.tex= QLineEdit(self.glWidget)
                self.tex.setText(value[3])
                self.gl.addWidget(self.item[i],i,2,1,1)
                self.gl.addWidget(self.tex,i,1,1,1)
                self.value=value
            i+=1
        self.buttonBox = QDialogButtonBox(self)
        self.buttonBox.setOrientation(Qt.Horizontal)
        self.buttonBox.setStandardButtons(QDialogButtonBox.Cancel|QDialogButtonBox.Ok)
        self.glWidget2 = QWidget(self)
        self.glWidget2.setGeometry(QRect(5, nb*20+y0*160+30, int(self.screenShape.width()*.15), 40))
        self.gl2 = QGridLayout(self.glWidget2)
        self.gl2.setContentsMargins(0,0,0,0)
        self.gl2.addWidget(self.buttonBox) #,nb,1,2,1)
        self.state = 'reject'
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
    def onColor(self): # OA added 22/11/19
        self.color = array(QColorDialog.getColor().getRgb())/255.
    
    def onFile(self,evt): #EV 07/02/20
        dlg = myFileDialog()
        title, ext, opt = self.value[0],self.value[1],self.value[2]
        file=self.value[3]
        self.fDir,self.fName=dlg.getsetFile(self.gui,title, ext, opt)
        file = str(self.fDir+self.fName)
        self.tex.setText(file)

    def getValues(self):
        self.exec_()
        nb = len(self.data)
        val = [0]*nb
        for i in range(nb):
            typ = self.data[i][1]
            if typ == 'Choice': val[i] = str(self.item[i].currentText())
            if typ == 'Text': val[i] = str(self.item[i].text())
            if typ == 'Check': val[i] = self.item[i].checkState()
            if typ == 'CheckList':
                val[i]=[]
                print(self.data[i][2][1])
                for j in range(self.item[i].count()):
                    if self.item[i].itemChecked(j): 
                        val[i].append(self.data[i][2][1][j])
            if typ == 'Textlong': 
                v0 = str(self.item[i].document().toPlainText())
                val[i] = v0.split('\n')   
            if typ == 'Color': val[i] = self.color # OA 22/11/19
            if typ == 'File' : val[i] =str(self.tex.text()) #EV 07/02/20
        if self.state =='accept' : return val
        else : return None

class CheckableComboBox(QComboBox):
	def __init__(self):
		super().__init__()
		self._changed = False
		self.view().pressed.connect(self.handleItemPressed)

	def setItemChecked(self, index, checked=False):
		item = self.model().item(index, self.modelColumn()) # QStandardItem object
		if checked: item.setCheckState(Qt.Checked)
		else: item.setCheckState(Qt.Unchecked)
        
	def handleItemPressed(self, index):
		item = self.model().itemFromIndex(index)
		if item.checkState() == Qt.Checked:
			item.setCheckState(Qt.Unchecked)
		else: item.setCheckState(Qt.Checked)
		self._changed = True
        
	def hidePopup(self):
		if not self._changed: super().hidePopup()
		self._changed = False

	def itemChecked(self, index):
		item = self.model().item(index, self.modelColumn())
		return item.checkState() == Qt.Checked
    
class myFileDialog: # Dialog to open or save a file 
    def __init__(self,opt='Open'):
        self.opt = opt
    def getsetFile(self,gui,title,filt,opt=False): #EV 07/02/20
        settings = QSettings("ORTI3D_team","ORTI3D")
        folder = settings.value("last_file")
        if folder=='' : folder = ' '
        dlg = QFileDialog()
        if self.opt=='Open':
            fileName = dlg.getOpenFileName(gui,title,folder,filter=filt)#, QFileDialog.DontUseNativeDialog)
        elif self.opt in ['New','Save']:
            dlg.setFileMode(QFileDialog.AnyFile)
            fileName = dlg.getSaveFileName(gui,title,folder,filter=filt)
        fileName = fileName[0] # OA 1/10 for python 3
        fName = fileName.split('/')[-1]
        fDir = fileName.replace(fName,'')
        if opt!=True : fName = fName.split('.')[0] #EV 07/02/20
        if fDir!='' : settings.setValue("last_file", fDir)# QVariant(QString(fDir)))
        return str(fDir),str(fName)

class myNoteBookCheck(QDialog): # Dialog to choose variable, used for Pest
    def __init__(self, gui,title, dicIn,opt=None):
        QDialog.__init__(self)
        self.setWindowTitle(title)
        self.gui,self.pages,self.layouts,self.dicIn = gui,{},{},dicIn
        self.dicOut = self.dicIn.copy()
        self.layout = QVBoxLayout(self)
        self.setGeometry(QRect(40, 40, 280,530))
        glWidget = QWidget(self)
        nb = QTabWidget(glWidget)
        nb.setGeometry(QRect(5, 5, 250,350))
        self.dwidget = {}
        for n in list(dicIn.keys()):
            if dicIn[n]==None:continue
            if len(dicIn[n])==0: continue
            nbChk = len(dicIn[n])
            pg = QWidget(nb)
            lay = QGridLayout(pg)
            self.layouts[n]= lay
            lay.setContentsMargins(0,0,0,0);lay.setSpacing(0)
            self.dwidget[n] = [0]*nbChk
            #if opt=='sort': dicIn[n] = self.sortList1(dicIn[n]) # OA modif 3/4 # EV 10/1/22
            for i in range(nbChk):
                ic = mod(i,1);il = int(i/1)
                ch = QCheckBox(nb); self.dwidget[n][i] = ch
                ch.setText(str(dicIn[n][i][0])) # OA modif 2/4
                s =(dicIn[n][i][1] in [1,2]) #EV 26/08/19 replaced "True" by True
                ch.setChecked(s) # OA modif 2/4 , EV 26/08/19 
                if n=='Layers': #EV 20/04/20
                    if dicIn[n][nbChk-1][0] == 'All layers':
                        ch.stateChanged.connect(self.onStateChange)
                lay.addWidget(ch,il,ic)
            scroll = QScrollArea()
            scroll.setWidget(pg)
            scroll.setWidgetResizable(True)
            #scroll.setFixedHeight(400)
            self.pages[n] = pg
            nb.addTab(scroll,n)
            pg.show()
        self.layout.addWidget(glWidget)
        self.buttonBox = QDialogButtonBox(self)
        self.buttonBox.setStandardButtons(QDialogButtonBox.Cancel|QDialogButtonBox.Ok)
        self.buttonBox.accepted.connect(self.accept1)
        self.buttonBox.rejected.connect(self.reject1)
        self.layout.addWidget(self.buttonBox)
        QMetaObject.connectSlotsByName(self)
    
    def onStateChange(self, state): #EV 20/04/20
        """ Only for multiplot to change the state of a button in notebook
        Layers whithout doing any action"""
        if state == Qt.Checked:
            if self.sender() == self.dwidget['Layers'][-1]:
                items = self.dwidget['Layers']
                for i in range(len(items)-1):
                    items[i].setChecked(False)
            else : 
                self.dwidget['Layers'][-1].setChecked(False)
        
    def sortList1(self,lst0): # OA added 2/4/19
        '''sort a list by the 1st item of each sublist'''
        l1 = [k[0] for k in lst0]
        indx = argsort(l1)
        indx = indx.astype('int')
        lst0b = array(lst0)
        return lst0b[indx]
        
    def accept1(self): 
        self.close(); self.state = 'accept'
    def reject1(self): 
        self.close(); self.state = 'reject'
    def apply(self): #EV 10/12/18
        self.state = 'accept'
        
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

class myNoteBook(QDialog): # Dialog used for add chemistry and pest parameters
    def __init__(self, gui,title, dicIn):
        QDialog.__init__(self)
        self.setWindowTitle(title)
        self.gui,self.pages,self.dicIn = gui,{},dicIn
        self.dicOut = self.dicIn.copy()
        layout = QVBoxLayout(self)
        self.screenShape = QDesktopWidget().screenGeometry()
        self.setGeometry(QRect(40, 60, int(self.screenShape.width()*.42),int(self.screenShape.height()*.6)))
        glWidget = QWidget(self)
        nb = QTabWidget(glWidget)
        nb.setGeometry(QRect(5, 20, int(self.screenShape.width()*.4),int(self.screenShape.height()*.50)))
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
            self.setRowHeight(il,int(self.screenShape.height()*.03))
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
                #print ('qtdl 269',il,ic,self.item(il,ic),self.type[ic])
                if self.type[ic]=='Text':
                    if self.item(il,ic) != None :
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
                if len(rowText)>self.rowCount(): # OA added if 1/8/19
                    self.setRowCount(len(rowText)) # OA 19/8/19
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
                self.clip.setText(s[:-1]) # OA added [:-1] 1/8/19
            
'''##########################" FOR ZONES  ###############################""
'''

class zoneDialog(QDialog): # Dialog for zone
    def __init__(self, parent, core,model,line, curzones, nb):
        QDialog.__init__(self)
        self.gui, self.core, self.model, self.line,  = parent.gui, core,model,line
        self.nb = nb
        self.glWidget = QWidget(self)
        self.screenShape = QDesktopWidget().screenGeometry()
        self.glWidget.setGeometry(QRect(5, 5, int(self.screenShape.width()*0.4),int(self.screenShape.height()*.5)))
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
        self.state = 'reject' #EV 8/12/21

    def showDialogAndDisconnect(self):
        self.show()
        self.gui.actionToggleEditing().triggered.disconnect(self.showDialogAndDisconnect)
               
    def saveCurrent(self):
        # saving the entered zone
        self.exec_()
        if self.state == 'accept': #EV 22/07/2019
            zp = self.zpanel;
            curdic = self.core.diczone[self.model]
            curzones = curdic.dic[self.line]
            curzones['name'][self.nb] = str(zp.name.text())
            media = zp.media.text()
            if '-' in media : #EV 24/10/18 & 22/07/19
                m1 = media.split('-')
                m2 = list(range(int(m1[0]),int(m1[1])+1))
            else :
                m2 = int(media)
            curzones['media'][self.nb] = m2;#onMessage(self.gui,str(zp.coords.getValues()['data']))
            if self.line != 'dis.1': # OA 20/11/18 this is to create modflow domain for qgis
                curzones['coords'][self.nb]= self.corrCoords(zp.coords.getValues()['data'])
            else : # OA 20/11/18
                curzones['coords'][self.nb]= zp.coords.getValues()['data']
            val0 = ''
            if self.typO: # if several values for one zone
                val0 = self.getOpt()
            if self.typS==1: # for transient values
                v0 = zp.valBox.getValues()['data']
                val = ''
                for b in v0 : 
                    val+=str(b[0])+' '+str(b[1])+'\n'
            else :
                v0 = zp.valBox.document().toPlainText()
                val = str(v0)            
            first = val.split('\n')[0]
            if first.count('.')>2:   #for pht3d zones 1.0.0.0
                val = val.replace('.','')
                val = val.replace(' ','')
                val = val.replace('\n','')
            curzones['value'][self.nb]= val0+val
            return 'OK'
        else : 
            return 'None'
        
    def accept1(self): 
        self.state = 'reject' 
        check = self.checkEntry() #EV 22/07/2019
        if check=='ok': #EV 22/07/2019
            self.close()
            self.state = 'accept'        

    def reject1(self): 
        self.close()
        self.state = 'reject'
    
    def checkEntry(self): #EV 22/07/2019
        zp = self.zpanel;
        media = zp.media.text()
        if not zp.name.text():
            onMessage(self,"Enter a name for the current zone");return
        if '-' in media :
            m1 = media.split('-')
            try :
                m2 = list(range(int(m1[0]),int(m1[1])+1))
            except ValueError: 
                m2 = 'None'
        else :
            try :
                m2 = int(media)
            except ValueError: 
                m2 = 'None'
        if m2=='None':
            onMessage(self,"Enter an integer number for media");return
        if self.typO: # if several values for one zone
            if self.getOpt() == 'None':
                onMessage(self,"Enter correct number in properties");return
        if self.typS==1: # for transient values
            v0 = zp.valBox.getValues()['data']
            for b in v0 : 
                try: float(b[0]);float(b[1])
                except ValueError:
                    onMessage(self,"Enter a number in transient data");return
        else : 
            v0 = zp.valBox.document().toPlainText() # main value
            if self.line not in ('obs.1','ghb.1','drn.1','riv.1'): #EV 29/10/19
                try :float(v0)
                except ValueError:
                    onMessage(self,"Enter a number in Zone value");return 
            if self.line == 'ph.4' : # EV 22/04/20
                nbsol = self.core.dicval['Pht3d']['ph.6'][1]
                if float(v0)>float(nbsol):
                    onMessage(self,"In Zone value, enter a number equal to or " 
                              +'\n'+ "less than the number of pht3d solutions")
                    return
        self.state = 'accept'
        return 'ok'
        
    def corrCoords(self,lcooI):
        '''change coordinates if they are out of the domain'''
        g = self.core.dicaddin['Grid']
        ex = (float(g['x1'])-float(g['x0']))/1e6 # OA 19/12/21
        ey = (float(g['y1'])-float(g['y0']))/1e6
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
            try: float(data[i][0])
            except ValueError: return 'None'
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
        self.coords = myNBpanelGrid(self.core.gui,self,dicCoo,size=(200,int(self.screenShape.height()*.3)))
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
            Parm1 = [x.split() for x in P1.split('\n')];# OA 6/9/19
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
#//////////////////////////// interpolation dialog /////////////
class IntpDialog(QDialog):
    '''This dialog provide different interpolation tools, options and 
    interpolation result'''
    def __init__(self,gui,core,opt):
        self.gui,self.core,self.opt= gui,core,opt
        QDialog.__init__(self,gui) 
        self.state = None
        self.setModal(False)
        self.setWindowTitle('Interpolation')
        screenShape = QtWidgets.QDesktopWidget().screenGeometry()
        #self.setGeometry(QRect(5, 5, screenShape.width()*.65, screenShape.height()*.5))
    ## main horizontal layout
        self.horizontalLayout = QHBoxLayout(self)
        self.horizontalLayout.setContentsMargins(10, 20, 10, 10)
    ## the left panel vertical layout vH1 
        self.verticalLayout = QVBoxLayout()
    ## choise of media 
        #txt = QLabel(self)
        #txt.setText('Media') 
        #self.gl.addWidget(txt,0,0,1,1)
        #self.line = QLineEdit(self)
        #self.line.setText('test')
        #self.gl.addWidget(self.line,0,1,1,1)
    ## getData 
        self.ch,self.mth, self.data = self.getData(self.opt)
    ## Interpolation method layout vH1.1
        txtMhd = QLabel()
        txtMhd.setText('Interpolation method') 
        font = QFont()
        font.setPointSize(9)
        font.setBold(True)
        txtMhd.setFont(font)
        self.verticalLayout.addWidget(txtMhd)
        self.plgroup = QComboBox()
        self.plgroup.addItems(['Kriging','Inverse distance',
                               'Thiessen polygons'])
        self.plgroup.setCurrentIndex(self.mth)
        self.plgroup.activated['QString'].connect(self.onChoiceOption)
        self.verticalLayout.addWidget(self.plgroup)
    ## spacer 1 vH1 
        self.verticalLayout.addSpacerItem(QSpacerItem(
                20, 50, QSizePolicy.MinimumExpanding, QSizePolicy.Minimum))
    ## Interpolation parameters vH1.2
        title2 = QLabel(self)
        title2.setText('Interpolation parameters') 
        font2 = QFont()
        font2.setPointSize(9)
        font2.setBold(True)
        title2.setFont(font2)
        self.verticalLayout.addWidget(title2)
        self.hlayout=QHBoxLayout()
        self.dialg = genericDialog(self.gui,'3D',self.data)
        self.hlayout.addWidget(self.dialg)
        self.dialg.gl2.removeWidget(self.dialg.buttonBox)
        self.dialg.buttonBox.close()
       ## Draw button
        self.drawButton = QPushButton(self)
        self.drawButton.setText('Show result')
        self.dialg.gl2.addWidget(self.drawButton)
        self.drawButton.clicked.connect(self.plotResult)
        self.verticalLayout.addLayout(self.hlayout,3)
    ## Save or recalculate layout vH1.3
        title3 = QLabel(self)
        title3.setText('Save options') 
        font3 = QFont()
        font3.setPointSize(9)
        font3.setBold(True)
        title3.setFont(font3)
        self.verticalLayout.addWidget(title3)
       ## choise save or recalculate
        self.choise = QComboBox(self)
        self.choise.addItems(['Save results as array',
                              'Recalculate each time files are written'])
        self.choise.setCurrentIndex(self.ch)
        self.verticalLayout.addWidget(self.choise) 
    ## spacer 2 vH1 
        self.verticalLayout.addSpacerItem(
                QSpacerItem(20, 50, QSizePolicy.MinimumExpanding,
                            QSizePolicy.Minimum))
       ## save button
        self.saveButton = QDialogButtonBox(self)
        self.saveButton.setOrientation(Qt.Horizontal)
        self.saveButton.setStandardButtons(QDialogButtonBox.Cancel|
                QDialogButtonBox.Save)
        self.verticalLayout.addWidget(self.saveButton, 
                                      alignment=Qt.AlignHCenter)
        self.saveButton.accepted.connect(self.accept1)
        self.saveButton.rejected.connect(self.reject1)
    ## add vertical layout vH1  
        self.horizontalLayout.addLayout(self.verticalLayout,2)
        self.horizontalLayout.addSpacerItem(QSpacerItem(
                20, 20, QSizePolicy.MinimumExpanding, QSizePolicy.Minimum))
    ## the right panel vertical layout
        self.verticalLayout2 = QVBoxLayout()
    ## parameter of variogram
        self.frame = QFrame()
        self.HLayout2 = QHBoxLayout()
        self.txt = QLabel(self)
        self.HLayout2.addWidget(self.txt,1)
    ## variogram matplotlib figure
        self.vario = Figure(tight_layout=True)
        self.cnv2 = FigureCanvas(self.vario)
        self.HLayout2.addWidget(self.cnv2,2)
        self.frame.setLayout(self.HLayout2)
        if self.plgroup.currentIndex() != 0:
            self.frame.hide()
    ## add variogram figure
        self.verticalLayout2.addWidget(self.frame,1)
    ## result matplotlib figure 
        self.figure = Figure(tight_layout=True,figsize=(4, 5), dpi=100) # EV 04/02/20 
        self.cnv = FigureCanvas(self.figure) 
    ## add matplotlib figure
        self.verticalLayout2.addWidget(self.cnv,3)
    ## add vertical layout 2
        self.horizontalLayout.addLayout(self.verticalLayout2,3)
        QMetaObject.connectSlotsByName(self)  #OA 1/6/19
    
    def getData(self,opt):
        if opt :
            val=opt
        else : val=[0,True,True,False,6,'spherical',1e-3,50,0,0]
        self.mth=val[0]
        self.ch=val[-1] 
        if self.mth == 0:
            self.data = [('Automatic parameter','Check',val[1]),
                        ('Log Values','Check',val[2]),
                        ('Plot Variogram','Check',val[3]),
                        ('nlags','Text',val[4]),
                        ('Variogram model','Choice',(val[5],['power','gaussian',
                                    'spherical','exponential'])),
                        ('Sill / Scale*','Text',val[6]),
                        ('Range / Exponent*','Text',val[7]),
                        ('Nugget','Text',val[8])
                        ]
        if self.mth==1:
            self.data = [('Power','Text',val[1])]
        if self.mth==2 :
            self.data = [('Smoothing number','Text',val[1])]
        return self.ch, self.mth, self.data
    
    def onChoiceOption(self):
        intpMtd=int(self.plgroup.currentIndex())
        try : self.dialg 
        except: pass 
        else : self.dialg.setVisible(False) 
        if intpMtd == 0 :
            data = [('Automatic parameter','Check',True),
                            ('Log Values','Check',True),
                            ('Plot Variogram','Check',False),
                            ('nlags','Text',6),
                            ('Variogram model','Choice',
                             ('spherical',['power','gaussian',
                                           'spherical','exponential'])),
                            ('Sill / Scale*','Text',1e-4),
                            ('Range / Exponent*','Text',30),
                            ('Nugget','Text',0)
                            ]
            self.frame.show()
        if intpMtd == 1 :
            data = [('Power','Text',1)]
            self.frame.hide()
        if intpMtd == 2 :
            data = [('Smoothing number','Text',1)]
            self.frame.hide()
        self.dialg = genericDialog(self.gui,'3D',data)
        self.hlayout.addWidget(self.dialg)
        self.dialg.gl2.removeWidget(self.dialg.buttonBox)
        self.dialg.buttonBox.deleteLater()
        del self.dialg.buttonBox 
        self.dialg.gl2.addWidget(self.drawButton)  
        
    def getParms(self):
        nb = len(self.dialg.data)
        parms = [0]*nb
        intpMtd = int(self.plgroup.currentIndex())
        for i in range(nb):
            typ = self.dialg.data[i][1]
            if typ == 'Choice': parms[i]= str(self.dialg.item[i].currentText())
            if typ == 'Text': parms[i] = str(self.dialg.item[i].text())
            if typ == 'Check': 
                if self.dialg.item[i].isChecked() : parms[i] = 1
                else : parms[i] = 0
            if typ == 'Textlong': 
                v0 = str(self.dialg.item[i].document().toPlainText())
                parms[i] = v0.split('\n')   
        parms.insert(0,intpMtd)
        return parms
        
    def plotResult(self):
        self.allOpt=self.getParms(); #print(self.allOpt)
        line = self.gui.currentLine
        media = self.gui.currentMedia 
        model = self.gui.currentModel
     ## Get and plot interpolation result
        value,mess,extent=self.core.runInterp(model,line,media,self.allOpt)
        self._ax=self.figure.add_subplot(1,1,1)
        if self.core.addin.mesh == None or self.core.getValueFromName(model,'MshType')==0 :
            myplot=self._ax.imshow(value,extent=extent)
        else :
            myplot=self._ax.tricontourf(self.core.addin.mesh.trg,value)
        self.figure.colorbar(myplot, ax=self._ax)
        self._ax.figure.canvas.draw()
        self.figure.clf() 
     ## Plot variogram
        if mess :
            if self.allOpt[0] == 0 :
                self._axv=self.vario.add_subplot(1,1,1)
                self.vario.clf()
                self._axv.figure.canvas.draw()
                if self.allOpt[3] == 1 :
                    lags=mess[1]
                    sem =mess[2]
                    var =mess[3]
                    self._axv=self.vario.add_subplot(1,1,1)
                    self._axv.plot(lags, sem, 'r*')
                    self._axv.plot(lags, var, 'k-')
                    self._axv.figure.canvas.draw()
                    self.vario.clf()
                self.txt.setText(mess[0])
            
    def saveResult(self):
        self.exec_()
        parms=self.getParms()
        ch=self.choise.currentIndex()
        parms.append(ch)
        if self.state != 'accept': parms=None
        return parms

    def accept1(self): 
        self.close(); self.state = 'accept' 
        
    def reject1(self): 
        self.close(); self.state = 'reject'

#///////////////////////////// import observation data /////////////
class impObsData(QDialog) :
    def __init__(self,gui,core,option):
        self.gui,self.core= gui,core
        self.option = option
        QDialog.__init__(self,gui) 
        self.setModal(False)
        self.setWindowTitle(self.option+' observation data')
        self.screenShape = QDesktopWidget().screenGeometry()
        self.setGeometry(QRect(40, 60, int(self.screenShape.width()*.42),int(self.screenShape.height()*.6)))
    ## main vertical layout
        self.verticalLayout = QVBoxLayout(self)
        self.verticalLayout.setContentsMargins(10, 20, 10, 10)
    ## label for instruction
        label = str('Set or copy paste here your '+self.option+' observation data')
        self.label = QtWidgets.QLabel(self)
        self.label.setMaximumSize(500, 24)
        self.label.setText(label)
        self.verticalLayout.addWidget(self.label)#, alignment=Qt.AlignHCenter)
    ## grid for paste data
        dicIn=self.getDicObs(self.option)
        self.nbg = myNBpanelGrid(self.gui,self,dicIn)
        #dicOut=self.setDicObs(self.option)
        self.verticalLayout.addWidget(self.nbg)
    ## button ok and cancel
        buttonBox = QDialogButtonBox(self)
        buttonBox.setStandardButtons(QDialogButtonBox.Cancel|QDialogButtonBox.Ok)
        buttonBox.accepted.connect(partial(self.setDicObs,self.option))
        buttonBox.rejected.connect(self.reject1)
        self.verticalLayout.addWidget(buttonBox)
        #QMetaObject.connectSlotsByName(self)
        
    def getDicObs(self,option):
        dicName=str('obs'+option)
        dic = self.core.dicaddin[dicName]
        if dic != {}:
            if option == 'Chemistry':
                dic=self.updateChem(dicName)
            return dic
        else :
            if option == 'Head':
                dicIn = {'cols':['\tWell\t','\tLayer\t','\tTime\t','\tHead\t'],'rows':['1','2','3'],'data':{}}
            if option == 'Tracer':
                dicIn = {'cols':['\tWell\t','\tLayer\t','\tTime\t','\tTracer\t'],'rows':['1','2','3'],'data':{}}
            if option == 'Chemistry':
                lname=self.getChemSol()
                cols=['\tWell\t','\tLayer\t','\tTime\t']+lname
                dicIn = {'cols':cols,'rows':['1','2','3'],'data':{}}   
            self.core.dicaddin[dicName]=dicIn
            return dicIn
    
    def getChemSol(self):
        lname=[]
        zchem=self.core.dicaddin['Chemistry']
        zname=zchem['Chemistry']['Solutions']['rows']
        for i in range(len(zname)) :
            if (zchem['Chemistry']['Solutions']['data'][i][0])!=False :
                lname.append(zname[i])
        return lname
    
    def updateChem(self,dicName):
        lname=self.getChemSol() 
        lname2=self.core.dicaddin[dicName]['cols'][3:] 
        if len(lname)>len(lname2):
            for i in range(len(lname)):
                if lname[i] not in lname2:
                    self.core.dicaddin[dicName]['cols'].insert(i+3,lname[i])
                    [self.core.dicaddin[dicName]['data'][x].insert(i+3,'') 
                    for x in range(len(self.core.dicaddin[dicName]['data']))]
        else :
            for i in range(len(lname2)): #EV 05/03/2019
                if lname2[i] not in lname:
                    self.core.dicaddin[dicName]['cols'].pop(i+3)
                    [self.core.dicaddin[dicName]['data'][x].pop(i+3) 
                    for x in range(len(self.core.dicaddin[dicName]['data']))]  
        return self.core.dicaddin[dicName]
        
    def setDicObs(self,option):
        dic2=self.nbg.getValues()
        dicName=str('obs'+option)
        nrow = len(dic2['data'])
        dic2['rows']=[str(x+1) for x in range(nrow)]
        if dic2 != None:
                self.core.dicaddin[dicName] = dic2
        else : return
        check=self.checkObs(dic2)
        if check == 'No': self.close()

    def reject1(self): 
        self.close()
        
    def checkObs(self,dic):
        #print('test',dic['data'])
        #if any(dic['data']) != False : #11/04/19
        zname=self.core.diczone['Observation'].dic['obs.1']['name']
        zobs=[dic['data'][i][0] for i in range(len(dic['data']))]
        m=''
        for i in range(len(zobs)):
            if zobs[i] not in zname: 
                m+= str(zobs[i])+' is not a model observation well\n'
        if m!='' :
            m2='Warning\n'+m+'\nDo you want modify your data?'
            resp=onQuestion(self.gui,m2)
            return resp
        else : self.close()
        #else : self.close()

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
        self.gl.setContentsMargins(1,1,1,1)
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

        self.butlist = QPushButton('User List',self.glWidget) # OA 1/10/19 added user list
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
        #if self.user.getValue()==False: return # OA removed 1/10/19
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
        self.gl.setContentsMargins(1,1,1,1)
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
from matplotlib.backends.qt_compat import QtCore, QtWidgets
from matplotlib.backends.backend_qt5agg import (
        FigureCanvas, NavigationToolbar2QT as NavigationToolbar)
from matplotlib.figure import Figure
class plotxy(QDialog):
    
    def __init__(self,gui,x,arry,legd,title,Xtitle,Ytitle,typ='-'):
        QDialog.__init__(self,gui) # OA 6/11/18 added gui
        self.setModal(False)
        self.setWindowTitle('Plot of results')
        layout = QVBoxLayout(self)
        self.gui,self.data,self.title,self.legd,self.Xtitle,self.Ytitle,self.x,self.arry = gui,3,title,legd,Xtitle,Ytitle,x,arry
        lab = QLabel(title)
        layout.addWidget(lab)
        glWidget = QWidget(self)
        glWidget.setGeometry(QRect(0, 0, 200, 300))
        self.cnv = FigureCanvas(Figure(figsize=(5, 3)))
        self._ax = self.cnv.figure.add_axes([0.1, 0.15, 0.7, 0.8])#subplots()
        #self._ax = self.cnv.add_axes([0.1, 0.1, 0.6, 0.75])
        x,arry = transpose(array(x,ndmin=2)),array(arry) # x en vertical
        # verify size of input vectors
        if len(shape(x))==1:
            x1 = ones((len(x),1))+0.;
        if len(shape(arry))==1:
            arry1 = ones((len(arry),1))+0.; arry1[:,0]=arry; arry=arry1*1.;ny=1
        else :
            [nt,ny] = shape(arry)
        x2 = x*1.; arry2 = arry*1.;
        # creer les lignes
        self._ax.plot(x2,arry2)
        self._ax.legend(self.legd,bbox_to_anchor=(1.05, 1),loc='upper left', borderaxespad=0.)
        self._ax.ticklabel_format(style='sci',scilimits=(-4,4),axis='both')
        self._ax.set_xlabel(self.Xtitle)
        self._ax.set_ylabel(self.Ytitle)
        layout.addWidget(self.cnv)
             
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
        self.show() # OA 6/11/18 replace exec by show
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

class plotxy_old(QDialog):
    
    def __init__(self,gui,x,arry,legd,title,Xtitle,Ytitle,typ='-'):
        QDialog.__init__(self)
        self.setModal(False)
        self.setWindowTitle(title)
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
