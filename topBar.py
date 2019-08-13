from os import sep
from .geometry import *
from .qtDialogs import *
from .config import *
    
class BaseTop:
    def __init__(self,gui,core):
        self.blind,self.gui,self.visu,self.core = {},gui,gui.visu,core
        for k in self.core.modelList: self.blind[k] = []
        self.blind['Mt3dms']=['btn.9','btn.10','uzt.3','uzt.4'],
        self.cfg = Config(core)
        self.gtyp = self.cfg.gtyp
        self.curVar = {}
        
    def modlistFromGroup(self,categ):
        """find the list of models in a category"""
        lshort = ['Sutra','Min3p','Opgeo']
        lshort2 = ['Sutr','Min3','Opge']
        if categ in lshort:
            lmodels = [x for x in self.core.modelList if x[:4]==categ[:4]]
        else :
            if 'series' in categ: # 27/7/19 this line and 3 more changed forUsg transp
                lmodels = ['Modflow','Mt3dms','Pht3d']
            elif 'USG' in categ:
                lmodels = ['Modflow','MfUsgTrans','Pht3d']
        if 'Observation' not in lmodels: 
            lmodels.append('Observation')
        return lmodels
        
    def modelFromLine(self,categ,line):
        """returns the model corresponding to a given line name"""
        lmodels = self.modlistFromGroup(categ)
        for md in lmodels:
            if line in list(self.core.dickword[md].lines.keys()): return md

    def testConditions(self,modName,lstL):
        """ test if the lines indices given in lstL sastify the condition
        also test if the data are like arrays"""
        indexout=[];
        for i,l in enumerate(lstL):
            cond = self.core.dickword[modName].lines[l]['cond']
            atype = self.core.dickword[modName].lines[l]['type'][0]
            if self.core.testCondition(modName,cond) and (atype[:3]=='arr'):
                indexout.append(i)
        return indexout

    def changeVisu(self):
        '''this modifies the visu when a button in topbar is used
        to be modified to consider qgis'''
        line = self.gui.currentLine
        media = self.gui.currentMedia
        self.visu.showVar(line, media)        
        
    def onImportZones(self):
        fdialg = myFileDialog()
        fileDir,fileName = fdialg.getsetFile(self.gui,'choose zone file','*.txt')
        self.core.importZones(fileDir,fileName,self.gui.currentModel,self.gui.currentLine)

    def onFormula(self):
        """opens a dialog to ask for python formula and executes them
        to get the value of the given keyword in the last line
        """
        ll = self.gui.currentLine
        formula = self.core.getFormula(self.gui.currentModel,ll)[0]; #onMessage(self.gui,formula)
        dialg = textDialog(self.gui,'input python formula',(340,300),str(formula))
        formula = dialg.getText()
        if formula != None :
            self.core.dicformula[self.gui.currentModel][ll]=[str(formula)]
            self.core.dictype[self.gui.currentModel][ll]=['formula']
               
    def onInterpolate(self):
        """creates the string for a specific formula"""
        ll = self.gui.currentLine
        model = self.gui.currentModel
        reg = 'Regular'
        if model[:5]=='Opgeo': reg = 'Unstruct'
        s='modName = \''+model+'\'\nline = \''+ll+'\'\nintp = 1\n'
        s+='#choices : interp. Kr, interp. ID; vtype: spher, gauss, gauss3\n'
        s+='opt = \'interp. Kr;vrange=22;vtype=@spher@\'\n'
        s+='value = block'+reg+'(self,modName,line,intp,opt,0)'
        self.core.dicformula[model][ll]=[s]
        self.core.dictype[model][ll]=['formula']
        #self.choiceT.setCurrentIndex(self.typeList.index('formula'))
        self.onFormula()
        
    def getCurVariable(self):
        """used to see the current variable for current medium"""
        mod = self.gui.currentModel
        line = self.gui.currentLine
        opt, iper, plane, section = None,0, 'Z',self.gui.currentMedia; #print 'guisho 275',line,self.curVar
#        if self.curVar.has_key(line): 
#            mat = self.curVar[line]*1
#        else:
        mat = self.core.getValueLong(mod,line,0)*1;#print 'topb 210',shape(mat)
        #self.curVar[line] = mat*1;
        X,Y = getXYmeshSides(self.core,plane,section)
        if self.core.addin.getModelGroup()=='Opgeo':
            return None,None,mat[0][0]
        if self.core.addin.getDim() in ['Radial','Xsection']:
            m2 = mat[-1::-1,0,:]
        else :
            if plane=='Z': m2 = mat[section,:,:] #-1 for different orientation in modflow and real world
            elif plane=='Y': m2 = mat[:,section,:]
            elif plane=='X': m2 = mat[:,:,section]
        return X,Y,m2
            
    def onZoneCreate(self, typeZone, xy):
        """ zone drawn in visu, we get coords
        and open the dialog to fill it
        """
        curzones = self.core.diczone[self.gui.currentModel]
        line = self.gui.currentLine
        curzones.addZone(line)
        iz = curzones.getNbZones(line)-1;
        curzones.setValue(line,'coords',iz,xy)
        curzones.setValue(line,'value',iz,' ')
        curzones.setValue(line,'media',iz,self.gui.currentMedia)
        #onMessage(self.gui,'topbar \n'+str(curzones.dic[line]))
        # create dialog
        dialg = zoneDialog(self, self.core,self.gui.currentModel,line, curzones.dic[line], iz)
        retour = dialg.saveCurrent()
        #self.gui.actions('zoneEnd')
        if retour != 'None':
            lz = curzones.dic[line]
            self.gui.visu.addZone(lz['media'][iz],lz['name'][iz],lz['value'][iz],lz['coords'][iz])
            if self.gtyp == 'qt': self.gui.modifBox.updateChoiceZone(line)
            self.core.dictype[self.gui.currentModel][line]=['zone']
            self.core.makeTtable();#print core.ttable
            #self.modifZones(line)
        else : #cancel
            curzones.delZone(line,iz)
            self.gui.visu.redraw()#line)
        if line=='obs.1':
            onames = self.core.diczone['Observation'].dic['obs.1']['name']
            self.gui.guiShow.setNames('Observation_Zone_L',onames)
        #self.modifZones(line)# adding a zone changes the view of the current variable
