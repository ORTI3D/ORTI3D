# -*- coding: utf-8 -*-

import os, sys
from config import *
from importExport import *
from qtPyConsole import *
from qtDialogs import *
import zipfile as zp
import requests # OA 1/10
import shutil

class Menus:
    def __init__(self, gui, core):
        self.gui,self.core = gui,core
        self.cfg = Config(core)
        self.dialogs = self.cfg.dialogs
        self.gtyp = self.cfg.gtyp

    def OnNew(self,evt=None):
        """creates a new file"""
        #self.askSave(evt)
        dlg = self.dialogs.myFileDialog('New')
        self.core.fileDir,self.core.fileName = dlg.getsetFile(self.gui,'New Model',"*.iqpht;*.orti")
        if self.core.fileDir == '': return
        self.core.addin.initAddin()
        self.core.initAll()
        if self.gtyp =='qt':
            self.gui.visu.setVisu(self.core)
            self.gui.updateTitle()

    def OnOpen(self,evt=None):
        """opens a file"""
        if self.core.fileDir!=None:
            self.askSave(evt)
        dlg = self.dialogs.myFileDialog('Open')
        fDir,fName =dlg.getsetFile(self.gui,'Open','*.iqpht;*.orti');
        if fName == '': return
        self.onOpen1(fDir,fName) # OA 16/12/18 to separate if file dir is already known
        if self.core.mfUnstruct: self.gui.onGridMesh('Mesh') # OA 18/2/20 removed from core, set here
        self.gui.onSetMediaNb(getNmedia(self.core),getNlayers(self.core))  # OA 14/3/21
        
    def onOpen1(self,fDir,fName): # OA 16/12/18 to be able to use only that part
        self.core.openModel(fDir,fName)
        a = self.core.makeTtable()
        self.core.addin.set3D()
        mtype = self.core.dicaddin['Model']['group']
        self.gui.guiShow.openModel() # OA 14/10
        if self.gtyp =='qt':
            if 'obs.1' in list(self.core.diczone['Observation'].dic.keys()):
                onames = self.core.diczone['Observation'].dic['obs.1']['name']
                #self.gui.guiShow.setNames('Observation_Zone_L',onames) #EV 19/12/18
            self.gui.visu.setVisu(self.core)
            self.gui.updateTitle()
        if self.gtyp=='qgis':
            self.gui.visu.removeOrtiLayers()
            self.gui.visu.zonesCore2qgs()
        self.gui.varBox.chooseCategory(mtype)
        self.gui.visu.initDomain()
        tl2 = self.core.ttable['wtimes'] #self.core.getTlist2()
        self.gui.guiShow.setNames('Model_Tstep_L',tl2)
        self.dialogs.onMessage(self.gui,'File opened')
            
    def OnSave(self,evt=None):
        if self.core.fileDir!=None:
            if self.gtyp =='qt':
                self.gui.updateTitle()
            if self.gtyp=='qgis':
                self.gui.visu.zonesQgs2core()
            self.core.saveModel(str(self.core.fileDir),str(self.core.fileName))
            self.dialogs.onMessage(self.gui,'File saved')
        else :
            self.OnSaveAs(evt)
            
    def OnSaveAs(self,evt=None):
        dlg = self.dialogs.myFileDialog('Save')
        fDir,fName = dlg.getsetFile(self.gui,'Save as','*.iqpht;*.orti')
        if fName == '': return
        self.core.saveModel(str(fDir),str(fName))
        if self.gtyp =='qt': #EV 27/11/18
            self.gui.updateTitle()
        if self.gtyp=='qgis':
            self.gui.visu.zonesQgs2core()
       
    def OnImportVersion1(self,evt=None):
        dlg = self.dialogs.myFileDialog()
        fDir,fName =dlg.getsetFile(self.gui,'Open',"*.ipht");#print fDir,fName
        importer = impFile(self.gui,self.core)
        importer.impVersion1(fDir,fName)
    '''    
    def OnImportModflowAscii(self,evt=None):
        dlg = self.dialogs.myFileDialog()
        fDir,fName =dlg.getsetFile(self.gui,'Open',"*.nam");#print fDir,fName
        importer = impAsciiModflow(self.core,fDir,fName)
        importer.readAll()
    '''    
    def OnImport3DgeomDis(self,evt=None):
        dlg = self.dialogs.myFileDialog()
        fDir,fName =dlg.getsetFile(self.gui,'Open',"*.disrw");#print fDir,fName
        importer = impFile(self.gui,self.core)
        bool = importer.imp3DgeomDis(self.core,fDir+os.sep+fName)
        if bool: self.dialogs.onMessage(self.gui,'File imported')

    def askSave(self,evt=None):
        if self.gtyp=='qgis': #EV 27/11/18
            if self.core.fileDir==None:return
        message = self.dialogs.onQuestion(self.gui,"Do you want to save the Orti file?") # OA modif 22/10/18
        if message == 'Yes':
            self.OnSave(evt)
        else : return
        #message.Destroy()
                
    def OnImportData(self,evt=None):
        """import external data to be used for representation"""  
        dlg = self.dialogs.myFileDialog()
        fDir,fName =dlg.getsetFile(self.gui,'Open data file',"*.txt")
        if fDir == None: return
        else : 
            self.data=self.core.importData(fDir,fName)
            self.dialogs.onMessage(self.gui,'Data imported')
        
    #def OnImportSolutions(self,evt=None):                               #EV 14/11/19
    #    """import a text file to store solutions"""
    #    dlg = self.dialogs.myFileDialog()
    #    fDir,fName =dlg.getsetFile(self.gui,'Open solutions',"*.txt")
    #    if fDir == None: return
    #    else : 
    #        self.core.importSolutions(fDir,fName)
    #        self.dialogs.onMessage(self.gui,'Solutions imported')
    #        #self.gui.OnMessage("Fichier donnees importe")
            
    def OnImportUserSpecies(self,evt=None):
        dlg = self.dialogs.myFileDialog()
        fDir,fName =dlg.getsetFile(self.gui,'Open solutions',"*.txt;*.out")
        if fDir == None: return
        elif fName == 'selected' :
            d = self.core.addin.pht3d.readSelectOut(fDir)
            #self.gui.guiShow.userSpecies = d
            self.gui.guiShow.setNames('Chemistry_User_L',list(d.keys()))
        else : 
            f1= open(fDir+os.sep+fName+'.txt')
            dicSp = {}
            for l in f1:
                if '=' in l: 
                    a,b=l.split('=');dicSp[a]=b
            self.gui.guiShow.setUserSpecies(dicSp)    
            nameBox = 'Chemistry_User_L'
            self.gui.guiShow.setNames(nameBox,list(dicSp.keys()))
        self.dialogs.onMessage(self.gui,'User species imported')
        
    def OnImportPostfixSpecies(self,evt=None): # oa added 30/6/19 (some modifs, added core.)
        if 'postfix.phrq' in os.listdir(self.core.fileDir):
            #print('ok')
            if 'selected.out' in os.listdir(self.core.fileDir):
                d = self.core.addin.pht3d.readSelectOut(self.core.fileDir)
                if self.gui != None:
                    self.gui.guiShow.setUserSpecies(d);#print 'core203',d
                    self.gui.guiShow.setNames('Chemistry_User_L',list(d.keys()))
                self.dialogs.onMessage(self.gui,'Postfix species imported')
            
    def OnExportParm(self,evt=None): # added 28/3/17 oa
        model,line,media = self.gui.currentModel,self.gui.currentLine,self.gui.currentMedia
        if not line : #EV 11/12/19
            self.dialogs.onMessage(self.gui,'Select a parameter to export')
        else :
            name = line.replace('.','')
            fname = self.core.fileDir+os.sep+ name
            data = self.core.getValueLong(model,line,0)[media2layers(self.core,media)[0]] # returns just the 1st layer of the considered medium
            savetxt(fname+'.txt',data)  
            self.dialogs.onMessage(self.gui,'File '+name+' saved')

    def OnExportResu(self,evt=None): # modif 28/3/17 oa for correct name
        if not self.gui.guiShow.curName : #EV 11/12/19
             self.dialogs.onMessage(self.gui,'Select a result to export')
        else :
            if self.gui.guiShow.curSpecies != True :
                name = self.gui.guiShow.curName+self.gui.guiShow.curSpecies;# OA modified 10/5/17
            else : 
                name = self.gui.guiShow.curName
            fname = self.core.fileDir+os.sep+ name
            data = self.gui.guiShow.data[-1]
            savetxt(fname+'.txt',data)  
            self.dialogs.onMessage(self.gui,'File '+name+' saved')
        
    def OnExportParmVtk(self,evt=None): # added 28/3/17 oa
        model,line = self.gui.currentModel,self.gui.currentLine
        if not line : #EV 11/12/19
            self.dialogs.onMessage(self.gui,'Select a parameter to export')
            return
        dlg = self.dialogs.myFileDialog('Save')
        fDir,fName = dlg.getsetFile(self.gui,'Save vtk',"*.vtk")
        if fName : 
            name = line.replace('.','')
            data = self.core.getValueLong(model,line,0)
            s = writeVTKstruct(self.core,model,data)
            f1=open(fDir+os.sep+fName+'.vtk','w');f1.write(s);f1.close()
            self.dialogs.onMessage(self.gui,'File '+fName+' saved')
        
    def OnExportResuVtk(self,evt=None): # modif 28/3/17 oa for correct name
        model = self.gui.currentModel
        data = self.gui.guiShow.arr3#;print('menu 151',shape(data))
        if not shape(data) : #EV 11/12/19
            self.dialogs.onMessage(self.gui,'Select a result to export')
            return
        dlg = self.dialogs.myFileDialog('Save')
        fDir,fName = dlg.getsetFile(self.gui,'Save vtk',"*.vtk")
        if fName :
            s = writeVTKstruct(self.core,model,data)
            f1=open(fDir+os.sep+fName+'.vtk','w');f1.write(s);f1.close()
            self.dialogs.onMessage(self.gui,'file '+fName+' saved')

    def OnExportVectorVtk(self,evt=None): # modif 28/3/17 oa for correct name
        model = self.gui.currentModel;
        data = self.gui.guiShow.vect3;
        if not shape(data) : #EV 11/12/19
            self.dialogs.onMessage(self.gui,'Select a result to export')
            return
        dlg = self.dialogs.myFileDialog('Save')
        fDir,fName = dlg.getsetFile(self.gui,'Save vtk',"*.vtk")
        if fName :
            s = writeVTKstruct(self.core,model,data,'vectors')
            f1=open(fDir+os.sep+fName+'.vtk','w');f1.write(s);f1.close()
            self.dialogs.onMessage(self.gui,'file '+fName+' saved')

    def OnExportResuVtkAll(self,evt=None): # modif 28/3/17 oa for correct name
        model,g = self.gui.currentModel,self.gui.guiShow
        name = g.curName
        dlg = self.dialogs.myFileDialog('Save')
        fDir,fName = dlg.getsetFile(self.gui,'Save vtk',"*.vtk")
        tlist = self.core.getTlist2()
        if fName :
            for it in range(len(tlist)):
                data = self.gui.guiShow.getArray3D(g.curGroup,name,it,g.curSpecies)
                s = writeVTKstruct(self.core,model,data)
                f1=open(fDir+os.sep+fName+str(it)+'.vtk','w')
                f1.write(s);f1.close()
            s = '{"file-series-version" : "1.0",\n "files" : \n [ \n'
            for it in range(len(tlist)):
                t = tlist[it]
                s += '{ "name" : "'+name+str(it)+'.vtk", "time" : '+str(t)+'},\n'
            s += '] \n }'
            f1=open(fDir+os.sep+fName+'.series','w');f1.write(s);f1.close()
        self.dialogs.onMessage(self.gui,'file '+fName+' saved')
            
    def OnHelpI(self,evt=None): #,lang):
        """calling help file"""
        os.startfile(self.gui.mainDir+os.sep+'doc'+os.sep+"interfaceHelp.chm")
    def OnHelpM(self,evt=None): #,lang):
        """calling help file"""
        os.startfile(self.gui.mainDir+os.sep+'doc'+os.sep+"modelsHelp.chm")
            
    def OnDownloadLast(self,evt=None):
        self.onDownload('master')
        
    def OnDownloadDev(self,evt=None):
        self.onDownload('develop')
        
    def OnDownloadLocal(self,evt=None): # Remove from orti3dGui.py EV 10/10/18
        dlg = self.dialogs.myFileDialog()
        fDir,fName =dlg.getsetFile(self.gui,'Open','*.zip');#print fDir,fName
        self.onDownload(fDir+os.sep+fName,'local')
        
    def OnHelpPy(self,evt=None):
        """starting a python console"""
        console=pyConsole(self.core,self.gui)
        console.exec_()

    def onDownload(self,fname,opt='web'):
        maindir=self.gui.mainDir
        dirdoc=maindir+os.sep+'doc'
        dirutil=maindir+os.sep+'utils'
        dirbin=maindir+os.sep+'bin'
        dirdoc=os.path.normpath(dirdoc)
        dirutil=os.path.normpath(dirutil)
        dirbin=os.path.normpath(dirbin)
        if self.cfg.typInstall=='python': 
            dirlib=maindir+os.sep+'ilibq'
        else : 
            dirlib=maindir+os.sep+'lib'+os.sep+'ilibq' # EV 23/10/18 new exe version
        dirlib=os.path.normpath(dirlib)
        lfu=os.listdir(dirutil)
        if 'newlib.zip' in lfu:
            os.system('copy '+dirutil+os.sep+'newlib.zip '+dirutil+os.sep+'oldlib.zip')
        f2=dirutil+os.sep+'newlib.zip'
        if opt=='web':
            htname = 'https://www.github.com/ORTI3D/ORTI3D/archive/'+fname+'.zip'
            r = requests.get(htname)
            with open(f2,"wb") as code: # OA 1/10
                code.write(r.content) # OA 1/10
        else :
            f2 = fname+'.zip'
        znew=zp.ZipFile(f2,'r')
        if self.cfg.typInstall=='python': #the python version
            znew.extractall(dirutil)
            #os.system('xcopy /Y '+dirutil+os.sep+'ORTI3D-'+fname+os.sep+'ilibq '+dirlib)   
            os.system('xcopy /Y /E '+dirutil+os.sep+'ORTI3D-'+fname+' '+dirlib) #EV 26/11/18
            print ('ok')
            for n in os.listdir(dirlib):
                if ('.chm' in n) or ('.pdf' in n): 
                    os.system('move '+dirlib+os.sep+n+' '+dirdoc)
                if ('.png' in n) or ('.dbs' in n): #EV 30/09/19 add .png
                    os.system('move '+dirlib+os.sep+n+' '+dirutil)
                if ('.exe' in n) : #EV 30/09/19 
                    os.system('move '+dirlib+os.sep+n+' '+dirbin)
        else : # the windows version
            znew.extractall(dirutil) # EV 23/10/18 new exe version
            os.system('xcopy /Y /E '+dirutil+os.sep+'ORTI3D-'+fname+' '+dirlib) #EV 26/11/18
            for n in os.listdir(dirlib): #EV 30/09/19 
                if ('.chm' in n) or ('.pdf' in n) : #EV 30/09/19 add .png
                    os.system('move '+dirlib+os.sep+n+' '+dirdoc)
                if ('.png' in n) or ('.dbs' in n): #EV 30/09/19 add .png
                    os.system('move '+dirlib+os.sep+n+' '+dirutil)
                if ('.exe' in n) : #EV 30/09/19 
                    os.system('move '+dirlib+os.sep+n+' '+dirbin)
            #zlib=zp.ZipFile(dirlib,'r'); #print 'menu dwnload 157', dirlib
            #zlib.extractall(maindir+os.sep+'temp')
            #zlib.close()
            #shutil.rmtree(maindir+os.sep+'temp'+os.sep+'ilibq')
            #znew.extractall(maindir+os.sep+'ilib1')
            #os.system('xcopy /Y '+maindir+os.sep+'ilib1'+os.sep+'ORTI3D-'+fname+os.sep+'iliblast '+maindir+os.sep+'temp'+os.sep+'ilibq')            
            #self.zip_folder(maindir+os.sep+'temp',maindir+os.sep+'zout.zip')
            #os.chdir(maindir)
            #os.system('del '+dirlib)
            #os.system('rename zout.zip library.zip')
            #shutil.rmtree('temp')
            #shutil.rmtree('ilib1')
        znew.close()
        if self.gui.gtyp!='qgis':
            self.dialogs.onMessage(self.gui,'lib changed, ORTi3D will stop, then restart it')
            sys.exit()
        else : 
            self.dialogs.onMessage(self.gui,'lib changed, please reloading qORTi3D plugin')
        
    def zip_folder(self,folder_path, output_path):
        parent_folder = os.path.dirname(folder_path)
        # Retrieve the paths of the folder contents.
        contents = os.walk(folder_path)
        zip_file = zp.ZipFile(output_path, 'w', zp.ZIP_DEFLATED)
        for root, folders, files in contents:
            for folder_name in folders:
                absolute_path = os.path.join(root, folder_name)
                relative_path = absolute_path.replace(parent_folder + '\\', '')
                relative_path = relative_path.replace('temp\\','')
                zip_file.write(absolute_path, relative_path)
            for file_name in files:
                absolute_path = os.path.join(root, file_name)
                relative_path = absolute_path.replace(parent_folder + '\\', '')
                relative_path= relative_path.replace('temp\\','')
                zip_file.write(absolute_path, relative_path)
        
    def OnBackVersion(self,evt=None):
        dirutil=self.gui.mainDir+os.sep+'utils'
        lf=os.listdir(dirutil)
        if 'oldlib.zip' not in lf: self.dialogs.onMessage('sorry no old lib')
        zin=zp.ZipFile(dirutil+os.sep+'oldlib.zip','r')
        zin.extractall(self.gui.mainDir+os.sep+'ilibq')
        self.dialogs.onMessage(self.gui,'lib changed, iPht3D will stop, then restart')
        self.gui.Destroy()
