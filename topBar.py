from os import sep
from .geometry import *
from .qtDialogs import *
from .config import *
from .myInterpol import IntpDialog # EV 19/02/20
    
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
        media = self.gui.currentMedia # EV 3/2/20
        formula = self.core.getFormula(self.gui.currentModel,ll,media) # EV 3/2/20
        dialg = textDialog(self.gui,'input python formula',(340,300),str(formula))
        formula = dialg.getText()
        if formula != None :
            self.core.dicformula[self.gui.currentModel][ll][media]=str(formula) # EV 4/2/20 add Media
            self.core.dictype[self.gui.currentModel][ll][media]='formula'# EV 3/2/20 add Media
               
    def onInterpolate(self):
        """creates the string for a specific formula"""
        ll = self.gui.currentLine
        media = self.gui.currentMedia # EV 3/2/20
        nmedia = getNmedia(self.core)
        model = self.gui.currentModel
        if ll in self.core.dicinterp[model]:  # EV 19/02/20
            try : self.core.dicinterp[model][ll][media]
            except IndexError :
                form=len(self.core.dicinterp[model][ll])
                nform=len(form)
                if nform<nmedia : form.extend(['']*(nmedia-form))
                opt=None
            else :
                opt = self.core.dicinterp[model][ll][media]
        else :
            self.core.dicinterp[model][ll] = ['']*nmedia
            opt=None
        m = IntpDialog(self.gui,self.core,opt)
        m.show()
        parms=m.saveResult()
        if parms != None :
            self.core.dicinterp[model][ll][media]=parms
            if parms[-1]== 1 :
                self.core.dictype[model][ll][media]='interpolate'
            else :
                cols,rows,value = self.core.save2array(model,ll,media)
                dlg = myFileDialog('Save')
                fDir,fName = dlg.getsetFile(self.gui,'Save to array','*.gvar') # EV 30/04/20
                if fName == '': file=None
                else :
                    file = str(fDir+fName)
                    f=open(file+'.gvar','w')
                    lcols=' '.join(map(str,cols));lrows=' '.join(map(str,rows))
                    f.write(lcols+'\n')
                    f.write(lrows+'\n')
                    savetxt(f,value)
                    f.close()
                if file :
                    self.core.dicarray[model][ll][media]=fName+'.gvar' #EV 05/05/20
                    self.core.dictype[model][ll][media]='importArray'
    
    def onImportArray(self):
        """Open a dialog to choose an array file and save the name in a dic EV 05/02/20"""
        ll = self.gui.currentLine
        media = self.gui.currentMedia # EV 3/2/20
        model = self.gui.currentModel
        nmedia = getNmedia(self.core)
        if ll in self.core.dicarray[model]: 
            try : self.core.dicarray[model][ll][media]
            except IndexError :
                arr=len(self.core.dicarray[model][ll])
                narr=len(arr)
                if narr<nmedia : arr.extend(['']*(nmedia-narr))
                f=''
            else :
                f = str(self.core.dicarray[model][ll][media])
        else :
            self.core.dicarray[model][ll] = ['']*nmedia
            f=''
        data=[('Media '+str(media),'File',['Choose Array',
               '*.asc;*.dat;*.txt;*.gvar',True,f])] #;*.vtk
        dialg = genericDialog(self.gui,'Choose Array',data)
        retour = dialg.getValues()
        fDir,fName=os.path.split(retour[0]) #print('ret',fName) #EV 05/05/20
        if retour :
            self.core.dicarray[model][ll][media]=fName #EV 05/05/20
        if retour == [''] : 
            self.core.dictype[model][ll][media]='one_value'
        
    def getCurVariable(self,var): #EV 26.11.20
        """used to see the current variable for current medium"""
        mod = self.gui.currentModel
        line = self.gui.currentLine
        opt, iper, plane, section = None,0, 'Z',self.gui.currentMedia;
        #print('topb gcurvar',mod,line,section,var)
#            mat = self.curVar[line]*1
        mdgroup = self.core.addin.getModelGroup() # OA 3/2/21 added
        if line in ['drn.1','riv.1','ghb.1']:
            if self.core.dictype[mod][line][section]=='importArray': #EV 26.11.20
                mat =  self.getTransientArray(line,var)[0];# OA 16/1/21 added 0
            else: 
                if 'USG' in mdgroup: # OA 16/1/21
                    mat=zone2mesh(self.core,mod,line,section,var=var)  # OA 16/1/21
                else :
                    mat=zone2grid(self.core,mod,line,section,opt=var,iper=0)
        else:
            if self.core.addin.getDim()=='3D' and 'USG' in mdgroup: # OA 3/2/21
                mat = zone2mesh(self.core,mod,line,media=section,iper=0)
                if self.core.dictype[mod][line][section]=='importArray':#EV 15/2/21 
                    mat = self.core.getValueLong(mod,line,0)[section]
            else :
                mat = self.core.getValueLong(mod,line,0)[section];#added section
        #self.curVar[line] = mat*1;
        X,Y = getXYmeshSides(self.core,plane,section)
        if 'USG' in mdgroup: # OA 20/11/20
            return None,None,mat  # OA 16/1/21 removed [0]
        #if self.core.addin.getDim() in ['Radial','Xsection']:
            #m2 = mat[-1::-1,0,:]
        #else :
            #if plane=='Z': m2 = mat[section,:,:] #-1 for different orientation in modflow and real world
            #elif plane=='Y': m2 = mat[:,section,:]
            #elif plane=='X': m2 = mat[:,:,section]
        return X,Y,mat #m2
    
    def getTransientArray(self,line,var): #EV 26.11.20
        line=line.split('.')[0]
        if line=='drn': nvar = 2
        if line=='riv': nvar = 3 
        if line=='ghb': nvar = 2
        im=self.gui.currentMedia
        flgU = False # flgU unstructured
        if self.core.addin.mesh == None: xx,yy=getXYmeshCenters(self.core,'Z',0)
        else : m = self.core.addin.mesh.getCenters();xx,yy = m[0],m[1];flgU=True
        ysign,gdx,gdy,arr = self.core.importGridVar(self.core.fileDir,line+str(im)+'.gvar')
        if ysign==-1:
            if gdy != None: gdy = gdy[-1::-1]*1
        intp = False;arr2 = [];grd = self.core.addin.getFullGrid()
        for iv in range(nvar): 
            if ysign==-1 :arr[iv]=arr[iv][::-1] #EV 14/01/2021
            arr2.append(linIntpFromGrid(self.core,grd,arr[iv],xx,yy,intp,gdx,gdy))
        #if ysign==-1 and flgU==False: #EV 11/12/20 #EV 14/01/2021
            #arr3=arr2[var]
            #arr3=[arr3[::-1]]
        #else: arr3=[arr2[var]]
        arr3=[arr2[var]]
        return arr3
        
    def onZoneCreate(self, typeZone, xy):
        """ zone drawn in visu, we get coords
        and open the dialog to fill it
        """
        curzones = self.core.diczone[self.gui.currentModel]
        line = self.gui.currentLine
        curzones.addZone(line)
        iz = curzones.getNbZones(line)-1;
        #xy = [(nice(a),nice(b)) for a,b in xy] # OA 26/8/19 to have short numbers
        #xy = [tuple([nice(x) for x in xy[i]]) for i in range(len(xy))]
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
            self.core.dictype[self.gui.currentModel][line][self.gui.currentMedia]='zone' #EV 13/02/20
            self.core.makeTtable();#print core.ttable
            #self.modifZones(line)
        else : #cancel
            curzones.delZone(line,iz)
            self.gui.visu.redraw(line)  # OA 28/1/21
        #if line=='obs.1': #EV 06/03/20
            #onames = self.core.diczone['Observation'].dic['obs.1']['name']
            #self.gui.guiShow.setNames('Observation_Zone_L',onames)
        #self.modifZones(line)# adding a zone changes the view of the current variable
