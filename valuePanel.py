import wx
#import wx.grid as wgrid
import os
from .modflowWriter import *
from .geometry import *
from . import core as corebase
from .myDialogs import *

class valuePanel(wx.Panel):
    """ main frame"""
    def __init__(self,parent,core,modName):
        wx.Panel.__init__(self,parent,-1,size=(400,500))
        self.fDir,self.fName= None,None
        self.core,self.modName = core,modName
        self.Mkword = core.dickword[modName] # contains groups and lines
        self.val = core.dicval[modName];# initiate values to 0
        self.zone = core.diczone[modName]
        self.array = core.dicarray[modName]
        self.formula = core.dicformula[modName]
        self.initStoredValues();# to get the conditions from the dictionnary

        grpList = self.core.getUsedModulesList(modName)
        self.chgroups = wx.Choice(self,-1,choices=grpList)
        self.Bind(wx.EVT_CHOICE,self.onChoiceGroup,self.chgroups)
        self.chlines = wx.Choice(self,-1,choices=[''])
        self.Bind(wx.EVT_CHOICE,self.onChoiceLine,self.chlines)
        self.boxkeys = boxKeys(self)
        frameSizer = wx.BoxSizer(wx.VERTICAL)

        frameSizer.Add(self.chgroups, 5) #, wx.EXPAND)
        frameSizer.Add(self.chlines, 5) #, wx.EXPAND)
        frameSizer.Add(self.boxkeys,80, wx.EXPAND)
        frameSizer.SetMinSize((280,440))
        frameSizer.Fit(self)
        #frameSizer.SetSizeHints(self)
        self.SetSizer(frameSizer)
        self.Layout()
            
    def setDictionnaries(self,dicval,diczone,dicarray,dicformula):
        self.val = dicval
        self.zone = diczone
        self.array = dicarray
        self.formula = dicformula
        self.initStoredValues()
        
    def initStoredValues(self):
        """at the model loadings look into the self.val dictionnary
        and gets all value that need to be stored as keywords
        """
        for ll in list(self.Mkword.lines.keys()):
            lk=self.Mkword.lines[ll];
            for i in range(len(lk['kw'])):
                kwd=lk['kw'][i]
                if ll not in list(self.val.keys()): continue
                if len(self.val[ll])>i: val=self.val[ll][i]
                else : val = lk['default'][i]
                if type(val) in [type('r'),type('r')]: continue
                if lk['type'][0][:3]=='arr': continue # do not store array
                kwd=kwd.split('(')[0]
                exec('self.'+kwd+'='+str(val))
    
    def changeStoredValues(self):
        """change the value of the keywords stored in the model"""
        line = self.currentLine
        lname = self.Mkword.lines[line]['kw']
        ltype = self.Mkword.lines[line]['type']
        for i in range(len(lname)):
            kwd=lname[i]
            a=kwd.split('(')
            if len(a)>=1: kwd=a[0]
            if ltype[i] in ['int','float','vecint','vecfloat','layint','layfloat']:
                val = self.val[line][i]
                exec('self.'+kwd+'='+str(val))

    def testConditions(self,lstL):
        """ test if the lines indices given in lstL sastify the condition"""
        lstout=[];
        for l in lstL:
            cond = self.Mkword.lines[l]['cond']
            if self.core.testCondition(self.modName,cond):
                if self.core.gui != None and self.Mkword.lines[l]['type'][0][:3]=='arr':
                    continue  # arrays, when gui is on are not set here
                lstout.append(l)
        return lstout
                       
    def onChoiceGroup(self,evt):
        item = self.FindWindowById(evt.GetId())
        n = item.GetStringSelection()
        groups,lines = self.Mkword.groups,self.Mkword.lines
        if n in groups: # make visible buttons for lines
            self.currentGroup=n
            #self.boxlines.Show(True);
            self.boxkeys.Show(False);
            #self.chlines.title.SetLabel(n)
            lname=[]
            lst1=self.testConditions(groups[n]) #select lines that satisfy the conditions
            for l in lst1:
                lname.append(l+'- '+lines[l]['comm'])
            self.chlines.Clear()
            for n in lname: self.chlines.Append(n)
            #self.boxlines.addButtons(lname,lst1)
            #self.boxlines.bxSizer.Layout()
            
    def onChoiceLine(self,evt):
        """ action when a line choice is clicked : change the interface"""
        lines=self.Mkword.lines
        groups=self.Mkword.groups
        item = self.FindWindowById(evt.GetId());
        n = item.GetStringSelection().split('-')[0]; 
        if n in list(lines.keys()):
            n=str(n);
            self.currentLine=n
            details=None;lname=lines[n]['kw'];nkw=len(lname)
            if 'detail' in lines[n]: details = lines[n]['detail']
            vtype=lines[n]['type']
            self.boxkeys.title.SetLabel(n)
            self.boxkeys.addButtons(lname,self.val[n],details,vtype)
            self.boxkeys.bxSizer.Layout()
            self.boxkeys.Show(True)
        self.Layout()
            
    def OnSetNewVal(self,evt):
        """sets the new values when user click on OK in key box"""
        values=self.boxkeys.getValues()
        self.val[self.currentLine]=[]
        for i in range(len(values)):
            self.val[self.currentLine].append(values[i])        
        lname=[]
        self.changeStoredValues();
        #readapt lines if condition modify them
        lst0=self.Mkword.groups[self.currentGroup]
        lst1=self.testConditions(lst0) #select lines that satisfy the conditions
        for l in lst1:
            lname.append(l+'- '+self.Mkword.lines[l]['comm'])
        self.chlines.Clear()
        for n in lname: self.chlines.Append(n)
        
    def setChoiceList(self,item,chlist):
        item.Clear()
        for ch in chlist: item.Append(ch)

#class boxButtons(wx.Panel):
#    def __init__(self,parent):
#        wx.Panel.__init__(self,parent,-1,size=(380,400))
#        self.parent=parent
#        self.bxSizer=wx.BoxSizer(wx.VERTICAL)
#        self.title= wx.StaticText(self,-1,' ')
#        self.bxSizer.Add(self.title,0)
#        self.butSizer=wx.BoxSizer(wx.VERTICAL)
#        self.bxSizer.Add(self.butSizer,0)
#        self.SetSizer(self.bxSizer)
#        
#    def addButtons(self,ltxt,lname):
#        self.butSizer.DeleteWindows();self.butSizer.Clear()
#        nk=len(ltxt)
#        for k in range(nk):
#            but=wx.Button(self,-1,ltxt[k],name=lname[k])
#            self.butSizer.Add(but,0)
#            self.Bind(wx.EVT_BUTTON,self.parent.OnChange,but)
#        self.bxSizer.Layout()
        
class boxKeys(wx.Panel):
    def __init__(self,parent):
        wx.Panel.__init__(self,parent,-1,size=(350,400))
        self.currentCol=0
        self.parent,self.modName = parent,parent.modName
        self.val = parent.val
        self.bxSizer=wx.BoxSizer(wx.VERTICAL) #|wx.EXPAND)
        tiSizer=wx.BoxSizer(wx.VERTICAL)
        self.title= wx.StaticText(self,-1,' ')
        tiSizer.Add(self.title,0)
        self.grdSizer = wx.FlexGridSizer(8,2,vgap=0,hgap=5)
        self.grd2Sizer = wx.FlexGridSizer(1,4,vgap=0,hgap=3)
        butOK=wx.Button(self,-1,'OK')
        self.Bind(wx.EVT_BUTTON,parent.OnSetNewVal,butOK)
        tunits=wx.StaticText(self,-1,' ')
        self.bxSizer.AddMany([(tiSizer,0),(self.grdSizer,0,wx.EXPAND,0),(self.grd2Sizer,0),
                (butOK,0),(tunits,0)])
        self.SetSizerAndFit(self.bxSizer)
        #self.Bind(wx.EVT_SIZE,self.onSize)
        self.Layout()
        
    def addButtons(self,lname,lval,ldetail,ltype):
        """it adds the new buttons after the user has clicked a choice"""
        self.lValBut=[];self.ltype=ltype;nb=len(lname);
        self.currentLine = self.parent.currentLine
        self.lval=lval
        self.grdSizer.DeleteWindows();self.grdSizer.Clear()
        self.grdSizer.SetCols(2)
        self.grdSizer.SetRows(nb)
        self.grd2Sizer.DeleteWindows();self.grd2Sizer.Clear()
        self.grd2Sizer.SetCols(4)
        self.grd2Sizer.SetRows(2)
        # adapt grid to insert Edit and formula buttons
        #if ltype[0][:3] in ['vec','arr']: self.grdSizer.SetCols(4)
        #else : 
        for i in range(nb):
            txt=lname[i]
            #put the array dimensions in parenthesis
            a=txt.split('(')
            if len(a)>1:
                txt=str(a[0])+'('
                b=a[1][:-1].split(',')
                for s in b:
                    exec('n=self.parent.'+s)
                    txt+=str(n)+','
                txt=txt[:-1]+')'
            #make the choice lists
            if len(lval)>i: curVal = lval[i]
            else : curVal = self.parent.Mkword.lines[self.currentLine]['default'][i]
            if ldetail not in [None,[]]: 
                if type(ldetail[i])==type([5]):
                    self.ltype[i]='choice'
                    txt+=' : '+ldetail[i][0]
                    val=wx.Choice(self,-1,choices=ldetail[i][1:])
                    val.SetSelection(curVal) #1: because 1 st val is title
                else :
                    txt+=' : '+ldetail[i]
                    val=wx.TextCtrl(self, -1, str(curVal))
            else :
                val=wx.TextCtrl(self, -1, str(curVal))
                    
            txtZ = wx.StaticText(self, -1, txt)
            self.grdSizer.Add(txtZ, 0)
            self.lValBut.append(val)
            self.grdSizer.Add(val,0)
        if ltype[i][:3] in ['vec','arr']:
            name = lname[i].split('(')[0]
            butE=wx.Button(self,-1,'E',name=name)
            self.Bind(wx.EVT_BUTTON,self.onEdit,butE)
            butF=wx.Button(self,-1,'F',name=name)
            self.Bind(wx.EVT_BUTTON,self.onFormula,butF)
            butZ=wx.Button(self,-1,'Z',name=name)
            self.Bind(wx.EVT_BUTTON,self.onZone,butZ)
            self.grd2Sizer.AddMany([(butE,0),(butF,0),(butZ,0)])
        self.Layout()
        
    def onSize(self,evt):
        self.Layout()
        
    def onEdit(self,evt):
        pass
    def onFormula(self,evt):
        """opens a dialog to ask for python formula and executes them
        to get the value of the given keyword in the last line
        """
        ll=self.currentLine
        curval=self.val[ll][ik];#print curval,type(curval)
        # find if the present value is a formula
        name='value='
        if curval=='formula': name=self.parent.formula[ll][ik]
        #make the dialog
        dialg = textDialog(self,'input python formula',(340,300),name)
        retour = dialg.ShowModal()
        if retour == wx.ID_OK:
            formula = dialg.getText();
            self.parent.formula[ll][ik]=formula
            self.val[ll][ik]='formula'
            self.lValBut[0].SetValue('formula')
        dialg.Destroy()

    def onZone(self,evt):
        """opens a dialog to ask for zones, save them and calc the array
        to get the value of the given keyword in the last line
        """
        zkeys = ['backg','number','name','coords','media','value','type']
        #make the dialog
        ll=self.currentLine
        if (ll in self.parent.zone)==False: 
            zone = {}
            for k in zkeys : zone[k] = []
        else : zone=self.parent.zone[ll]
        dialg = zoneDialog(self,'zone',(340,400),zone)
        retour = dialg.ShowModal()
        if retour == wx.ID_OK:
            dialg.saveCurrent();
            zone = dialg.getValues()
            self.parent.zone[ll] = zone # saving zone
            self.parent.val[ll]=['zone'] # tell the model that you use zones
            self.parent.array[ll]=zone2grid(self.parent.core,self.modName,ll,zone)# arranging grid from zone
            self.lValBut[0].SetValue('zone')
        dialg.Destroy()
            
    def setValues(self):
        nb=len(self.lval)
        for i in range(nb):
            but=self.lValBut[i];
            val=self.lval[i]
            if self.ltype[i]!='choice': but.SetValue(str(val))
            else: b.SetSelection(val)

    def getValues(self):
        nb=len(self.lValBut)
        self.lval=[]
        for i in range(nb):
            but=self.lValBut[i]
            #val = self.lval[i]
            if self.ltype[i]=='choice':
                self.lval.append(but.GetSelection())
                continue
            if but.GetValue() not in ['formula','zone','array']:
                if self.ltype[i] in ['int','vecint','arrint','layint','layfloat']:
                    self.lval.append(int(but.GetValue()))
                elif self.ltype[i] in ['float','vecfloat','arrfloat']:
                    self.lval.append(float(but.GetValue()))
                else : 
                    self.lval.append(but.GetValue())
            else :  
                self.lval.append(but.GetValue())
        return self.lval


