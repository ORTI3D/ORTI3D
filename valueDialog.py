# -*- coding: utf-8 -*-

import os
from geometry import *
from qtValueDialog import *

from defSPH import *
from configSPH import *


def onMessage(gui,text):  QMessageBox.information(gui,"Info",text)

class valueDialog:
    """this class allows to create the dialog to enter parameter values
    it uses the keyword dictionnaries to create the structure
    and also to test conditions for entering some specific values"""
    def __init__(self,gui,title,core,modName):

# note SPH january 2024:
# there is already a 'gui.core' 'and a 'gui.core.gui'... I'm lost
# so what's new in this new argument core ?!

# what I finally think (after zillions of iterations):
# 1. in the main GUI it is possible to choose the 'Model group'
# 2. this 'Model Group' can be 'Modflow series', 'Modflow USG' or 'Min3p', the other groups are in fact given up and should be removed from the lists
# 3. the 'Model Group' is stored in gui.core.Model['group'] ('gui.core.Model' being a dictionnary)
# 4. other parameters of the model are dimension (which can be 2D or 3D) and type (confined, unconfined, etc.) have no interest here
# 5. self.currentLine (probably created in the 'onChoiceLine' method) contains 'model.?' info, so removing the . and the number gives the model, no need to catch it
# BUT model may be in capital letters or not, depending on the model group: for 'Modflow series' and 'Modflow USG' the model is in capital letters
# for 'Min3p' the model is small letters for 'Flow', 'Transport' and 'Chemistry' families, for the 'Pest' family the models are always in capital letters

        if modifSPH2023_debugLevel>20:
            self.sph_modelGroupForBackup = modName
            if modifSPH2023_useInspect:
                import inspect
                print('debug 240129e:', inspect.getmembers(gui.core))
            print('end debug 240129e')

        if modName == 'Modflow_USG': modName = 'Modflow'
        self.butNb, self.fDir, self.fName = 500, None, None
        self.gui, self.core, self.modName = gui, core, modName
        self.Mkword = core.dickword[modName]  # contains groups and lines
        self.val = core.dicval[modName];  # initiate values to 0
        self.zone = core.diczone[modName]
        self.array = core.dicarray[modName]
        self.formula = core.dicformula[modName]
        self.interp = core.dicinterp[modName]  # EV 20/02/20
        self.initStoredValues();  # to get the conditions from the dictionnary
        # self.currentGroup,self.currentLine = self.gui.currentGroup,self.gui.currentLine

        if modifSPH2023_debugLevel > 10:
            print('debug SPH 240118c:', gui.currentModel)

        self.dialg = qtValueDialog(self, gui, core, modName)

        # some lines will not be shown (they are in the addin)
        self.blind = ['DELR', 'DELC', 'TOP', 'BOTM', 'PERLEN', 'NROW', 'NCOL', 'NLAY', 'NPER', 'WELLS',
                      'NCOMP', 'MCOMP', 'GCOMPN', 'KCOMPN', 'HTOP', 'DZ', 'PRSTY', 'ICBUND', 'SCONC', 'MTRECH',
                      'SP1', 'SP2', 'RC1', 'RC2', 'SWC', 'SDH', 'TLAYCON',
                      'LAYCBD']  # OA 24/5/19 removed al, trpt, trpv, dmcoef... EV 25/09/19 add TLAYCON EV 15/11/2019 add LAYCBD


# print self.Mkword.lines

    def show(self):
        self.dialg.exec_()


    def setDictionnaries(self, dicval, diczone, dicarray, dicformula, dicinterp):  # EV 20/02/20
        self.val = dicval
        self.zone = diczone
        self.array = dicarray
        self.formula = dicformula
        self.interp = dicinterp  # EV 20/02/20
        self.initStoredValues()


    def initStoredValues(self):
        """at the model loadings look into the self.val dictionnary
        and gets all value that need to be stored as keywords
        """
        for ll in list(self.Mkword.lines.keys()):
            # print 'valdlg l47',ll
            lk = self.Mkword.lines[ll];
            for i, kwd in enumerate(lk['kw']):
                if ll not in list(self.val.keys()): continue
                if lk['type'][:3] == 'lay':
                    val = self.val[ll]  # OA added 18/9/19
                else:
                    if i < len(self.val[ll]):
                        val = self.val[ll][i]
                    else: # value of the given keyword does not exist in val
                        val = lk['default'][i]
                        self.val[ll].append(val)
                if type(val) in [type('r'), type('r')]: continue
                if lk['type'][0][:3] == 'arr': continue  # do not store array
                kwd = kwd.split('(')[0]
                setattr(self, kwd, str(val))  # OA 3/10/18


    def changeStoredValues(self):
        """change the value of the keywords stored in the model"""
        l = self.currentLine
        lname = self.Mkword.lines[l]['kw'];
        for i, kwd in enumerate(lname):
            a = kwd.split('(')
            if len(a) >= 1: kwd = a[0]
            val = self.val[l][i]
            if type(val) not in [type('r'), type('r')]:
                setattr(self, kwd, str(val));  # OA 3/10/18


    def testConditions(self, lstL):
        """ test if the lines indices given in lstL satisfy the condition"""
        lstout = [];
        for l in lstL:
            cond = self.Mkword.lines[l]['cond'];  # print 'valudiag 75',l,cond
            ltyp = self.Mkword.lines[l]['type'][0]  # OA added 25/10 + arr condition below
            if (self.core.testCondition(self.modName, cond)) & (ltyp[:3] != 'arr'):
                lstout.append(l);  # print 'valudiag 75',l,cond, 'True'
        # print 'valdilg 80',lstout
        return lstout


    def onChoiceGroup(self, name):
        '''action when a group is chosen, changes the line list'''
        # onMessage(self.gui,str(name))
        groups, lines = self.Mkword.groups, self.Mkword.lines
        if name in groups:  # make visible buttons for lines
            self.currentGroup = name
            lname = []
            lst1 = self.testConditions(groups[name])  # select lines that satisfy the conditions
            for l in lst1:
                lname.append(l + '- ' + lines[l]['comm'])
            self.showBox(self.dialg.boxkeys, False);
            self.changeCombo(self.dialg.chlines, lname)


    def onChoiceLine(self, name):
        """ action when a line choice is clicked : change the interface"""
        lines = self.Mkword.lines
        comm = name.split('-')[1]  # EV 25/09/19
        name = name.split('-')[0]
        if name in list(lines.keys()):
            n = str(name);
            self.currentLine = n
            if 'detail' in lines[n]:
                details = lines[n]['detail']
            else:
                details = [None] * len(self.val[n])
            if len(details) == 0: details = [None] * len(self.val[n])
            # print(n,lines[n]['type'])
            self.changeButtons(name, comm, lines[n]['kw'], self.val[n], details, lines[n]['type']);
            # print('valueD l102',self.val[n])


    def OnSetNewVal(self, evt=''):
        """sets the new values when user click on OK in key box"""

        # note SPH 240122: more convenient to insert the save here rather than in getValues...
        # indeed for the save mecanism it looks equivalent but for load...
        if modifSPH2023_debugLevel > 15:
            if modifSPH2023_inspect:
                import inspect
                print('*************************************')
                print('debug 240119b:', inspect.getmembers(self.gui.core))
                print('*************************************')
            print('debug 240119c', hasattr(self, 'currentGroup'), self.currentGroup)
            print('    ', hasattr(self.gui, 'currentGroup'))
            print('    ', hasattr(self.core, 'currentGroup'))
            print('    ', hasattr(self.gui.core, 'dicaddin'),(self.gui.core.dicaddin['Model'])['group'])

# note the change of self.val was already done here !!! self.val gives finalValues !!! very weird !!!
            #print('debug 240130b:', self.val[self.currentLine])

            # modif SPH january 2024: add signature (empty before)
     # useless ?!
                # if modifSPH2023_enableTrySPH:
                #import copy
                # canceled
                #sph_addSave(calframe=calframe,values=self.values,backupList='dev',mode='pickle')#modifSPH2023_saveList)
                #initialValues = copy.deepcopy(self.val[self.currentLine])
            #values=self.dialg.boxkeys.getValues(sph_groupNameForBackup=self.sph_groupNameForBackup)

        values=self.dialg.boxkeys.getValues()
        #;print ('vdialg, setnew',self.currentLine,values)

        """
        note SPH:
        Modflow series, Flow ==> self.sph_modelGroupForBackup == 'Modflow'
        Modflow USG, Flow ==> self.sph_modelGroupForBackup == 'Modflow_USG'
        Min3p, Flow ==> self.sph_modelGroupForBackup == 'Min3pFlow'
        
        Modflow series, Transport ==> self.sph_modelGroupForBackup == 'Mt3dms'
        Modflow USG, Transport ==> self.sph_modelGroupForBackup == 'MfUSGTrans'
        Min3p, Transport ==> self.sph_modelGroupForBackup == 'Min3pTrans'
        
        Modflow series, Chemistry ==> self.sph_modelGroupForBackup == 'Mt3dms'
        Modflow USG, Chemistry ==> self.sph_modelGroupForBackup == 'MfUSGTrans'
        Min3p, Chemistry ==> self.sph_modelGroupForBackup == 'Min3pChem'
        
        ...
        
        """

        if modifSPH2023:
            # modif 240209
            # sph_addSave(currentLine=self.currentLine,
            #             modelGroup=(self.gui.core.dicaddin['Model'])['group'],  # sph_modelGroupForBackup,
            #             initialValues=self.dialg.boxkeys.previousValues,
            #             # self.val[self.currentLine],warning: self.val is already updated !!!
            #             finalValues=values,
            #             mode=sph_possibleBackupFormat.pickle,
            #             )
            sph_addSave(currentLine=self.currentLine,
                        modelGroup=(self.gui.core.dicaddin['Model'])['group'],#sph_modelGroupForBackup,
                        initialValues=self.dialg.boxkeys.previousValues,#self.val[self.currentLine],warning: self.val is already updated !!!
                        finalValues=values,
                        mode=modifSPH2023_backupFormat,
                        backupList=modifSPH2023_saveList
                        )

        for i in range(len(values)):
            #print (self.val[self.currentLine][i])
            self.val[self.currentLine][i]=values[i]
        names = []
        self.changeStoredValues()
        #readapt lines if condition modify them
        lst0=self.Mkword.groups[self.currentGroup]
        lst1=self.testConditions(lst0) #select lines that satisfy the conditions
        for l in lst1:
            names.append(l+'- '+self.Mkword.lines[l]['comm'])
        self.changeCombo(self.dialg.chlines,names)
        self.dialg.boxkeys.setVisible(False) # OA 10/09/2018
        #readapt unconfined / confined for lpf.2
        if self.currentLine == 'lpf.2':  #EV 25/09/19
            if set(values)=={'0'}:
                self.core.dicaddin['Model']['type']='Confined'
            elif set(values)=={'1'}:
                self.core.dicaddin['Model']['type']='Unconfined'
            else : self.core.dicaddin['Model']['type']='Mix (for 3D only)'

    def showBox(self,box,bool):
        box.setVisible(bool)

    def changeCombo(self,comboName,names):
        comboName.clear()
        for n in names:
            comboName.addItem(n)

    def changeButtons(self,title,comm,names,values,details,types):  #EV 25/09/19
        self.showBox(self.dialg.boxkeys,True)
        self.dialg.boxkeys.addButtons(title,comm,names,values,details,types) #EV 25/09/19
#        if self.gtyp=='qt':
#            self.dialg.boxkeys.title.SetLabel(n)

    def makeButton(self,name,value,detail,typ):
        #print('vdlg 145',name,typ)
        # find the dimension of the array
        txt = name;
        a=txt.split('(')
        if len(a)>1:
            txt=str(a[0])+'('
            b=a[1][:-1].split(',')
            for s in b:
                txt+=str(getattr(self,s))+','  # OA 3/10/18
            txt=txt[:-1]+')'
        #make the choice lists
        curVal = value
        bselect = None
        det1 = '' # OA 20/3/20 this and 2 follow lines
        if detail not in [None,[]]:
            if type(detail) != type([]) : det1 = detail
        if typ == 'choice' :
            txt+=' : '+detail[0]
            bcontent = detail[1:] #1st item is title line
            bselect = curVal
        elif typ == 'textlong' :
            txt+=' : '+det1 # OA 20/3/20 detail to det1
            bcontent = str(curVal)
#        else :
#            typ = 'text'
#            txt+=' : '+det1 # OA 20/3/20 detail to det1
#            bcontent = str(curVal)
        else :  # OA 20/03/20 readded
            if typ[:3] != 'lay':
                bcontent = str(curVal)
                typ = 'text'
            else :
                bcontent = curVal
            txt += ' : '+det1
        return txt,bcontent,bselect,typ # OA 6/11/18 changed order
