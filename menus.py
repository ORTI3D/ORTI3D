# -*- coding: utf-8 -*-

import os
from .config import *
from .importExport import *
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
        self.core.openModel(fDir,fName)
        a = self.core.makeTtable()
        listSpec = self.core.addin.chem.getListSpecies() # just the names
        self.core.addin.set3D()
        mtype = self.core.dicaddin['Model']['group']
        self.gui.guiShow.init()
        if self.gtyp =='qt':
            if 'obs.1' in list(self.core.diczone['Observation'].dic.keys()):
                onames = self.core.diczone['Observation'].dic['obs.1']['name']
                self.gui.guiShow.setNames('Observation_Zone_L',onames)
            self.gui.visu.setVisu(self.core)
            self.gui.updateTitle()
        if self.gtyp=='qgis':
            self.gui.visu.zonesCore2qgs()
        self.gui.varBox.chooseCategory(mtype)
        self.gui.visu.initDomain()
        tl2 = self.core.getTlist2()
        self.gui.guiShow.setNames('Model_Tstep_L',tl2)
        self.gui.guiShow.setChemSpecies(listSpec)
        self.dialogs.onMessage(self.gui,'file opened')
            
    def OnSave(self,evt=None):
        if self.core.fileDir!=None:
            if self.gtyp =='qt':
                self.gui.updateTitle()
            if self.gtyp=='qgis':
                self.gui.visu.zonesQgs2core()
            self.core.saveModel()
            self.dialogs.onMessage(self.gui,'file saved')
        else :
            self.OnSaveAs(evt)
            
    def OnSaveAs(self,evt=None):
        dlg = self.dialogs.myFileDialog('Save')
        fDir,fName = dlg.getsetFile(self.gui,'Save as','*.iqpht;*.orti')
        if fName == '': return
        self.core.saveModel(str(fDir),str(fName))
       
    def OnImportVersion1(self,evt=None):
        dlg = self.dialogs.myFileDialog()
        fDir,fName =dlg.getsetFile(self.gui,'Open',"*.ipht");#print fDir,fName
        importer = impFile(self.gui,self.core)
        importer.impVersion1(fDir,fName)
        
    def OnImportModflowAscii(self,evt=None):
        dlg = self.dialogs.myFileDialog()
        fDir,fName =dlg.getsetFile(self.gui,'Open',"*.nam");#print fDir,fName
        importer = impAsciiModflow(self.core,fDir,fName)
        importer.readAll()
        
    def askSave(self,evt=None):
        message = self.dialogs.onQuestion(self.gui,"Do you want to save the project?")
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
            self.core.importData(fDir,fName)
            self.dialogs.onMessage(self.gui,'Data imported')
        
    def OnImportSolutions(self,evt=None):
        """import a text file to store solutions"""
        dlg = self.dialogs.myFileDialog()
        fDir,fName =dlg.getsetFile(self.gui,'Open solutions',"*.txt")
        if fDir == None: return
        else : 
            self.core.importSolutions(fDir,fName)
            self.dialogs.onMessage(self.gui,'Solutions imported')
            #self.gui.OnMessage("Fichier donnees importe")
            
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
        self.dialogs.onMessage(self.gui,'User Species imported')
            
    def OnExportParm(self,evt=None): # added 28/3/17 oa
        model,line,media = self.gui.currentModel,self.gui.currentLine,self.gui.currentMedia
        name = line.replace('.','')
        fname = self.core.fileDir+os.sep+ name
        data = self.core.getValueLong(model,line,0)[media2layers(self.core,media)[0]] # returns just the 1st layer of the considered medium
        savetxt(fname+'.txt',data)  
        self.dialogs.onMessage(self.gui,'file '+name+' saved')

    def OnExportResu(self,evt=None): # modif 28/3/17 oa for correct name
        if self.gui.guiShow.curSpecies != True :
            name = self.gui.guiShow.curName+self.gui.guiShow.curSpecies;# OA modified 10/5/17
        else : name = self.gui.guiShow.curName
        fname = self.core.fileDir+os.sep+ name
        data = self.gui.guiShow.data[-1]
        savetxt(fname+'.txt',data)  
        self.dialogs.onMessage(self.gui,'file '+name+' saved')
        
    def OnExportParmVtk(self,evt=None): # added 28/3/17 oa
        dlg = self.dialogs.myFileDialog('Save')
        fDir,fName = dlg.getsetFile(self.gui,'Save vtk',"*.vtk")
        model,line,media = self.gui.currentModel,self.gui.currentLine,self.gui.currentMedia
        name = line.replace('.','')
        data = self.core.getValueLong(model,line,0)
        s = writeVTKstruct(self.core,data)
        f1=open(fDir+os.sep+fName,'w');f1.write(s);f1.close()
        self.dialogs.onMessage(self.gui,'file '+fName+' saved')
        
    def OnExportResuVtk(self,evt=None): # modif 28/3/17 oa for correct name
        dlg = self.dialogs.myFileDialog('Save')
        fDir,fName = dlg.getsetFile(self.gui,'Save vtk',"*.vtk")
        data = self.gui.guiShow.arr3;print('menu 151',shape(data))
        s = writeVTKstruct(self.core,data)
        f1=open(fDir+os.sep+fName,'w');f1.write(s);f1.close()
        self.dialogs.onMessage(self.gui,'file '+fName+' saved')
            
    def OnHelp(self,evt=None): #,lang):
        """calling help file"""
        os.startfile(self.gui.mainDir+os.sep+'doc'+os.sep+"iPht3dDoc_En.chm")
            
    def OnDownloadLast(self,evt=None):
        self.onDownload('master')
        
    def OnDownloadDev(self,evt=None):
        self.onDownload('develop')
        
    def OnDownloadLocal(self,evt=None): # Remove from orti3dGui.py EV 10/10/18
        dlg = self.dialogs.myFileDialog()
        fDir,fName =dlg.getsetFile(self.gui,'Open','*.zip');#print fDir,fName
        self.onDownload(fDir+os.sep+fName,'local')
        
    def onDownload(self,fname,opt='web'):
        maindir=self.gui.mainDir
        dirdoc=maindir+os.sep+'doc'
        dirutil=maindir+os.sep+'utils'
        if self.cfg.typInstall=='python': 
            dirlib=maindir+os.sep+'ilibq'
        else : 
            dirlib=maindir+os.sep+'library.zip'
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
            os.system('xcopy /Y '+dirutil+os.sep+'ORTI3D-'+fname+os.sep+'ilibq '+dirlib)            
            for n in os.listdir(dirlib):
                if ('.chm' in n) or ('.pdf' in n): 
                    os.system('move '+dirlib+os.sep+n+' '+dirdoc)
                if ('.gif' in n) or ('.dbs' in n): 
                    os.system('move '+dirlib+os.sep+n+' '+dirutil)
        else : # the windows version
            zlib=zp.ZipFile(dirlib,'r'); #print 'menu dwnload 157', dirlib
            zlib.extractall(maindir+os.sep+'temp')
            zlib.close()
            shutil.rmtree(maindir+os.sep+'temp'+os.sep+'ilibq')
            znew.extractall(maindir+os.sep+'ilib1')
            os.system('xcopy /Y '+maindir+os.sep+'ilib1'+os.sep+'ORTI3D-'+fname+os.sep+'iliblast '+maindir+os.sep+'temp'+os.sep+'ilibq')            
            self.zip_folder(maindir+os.sep+'temp',maindir+os.sep+'zout.zip')
            os.chdir(maindir)
            os.system('del '+dirlib)
            os.system('rename zout.zip library.zip')
            shutil.rmtree('temp')
            shutil.rmtree('ilib1')
        znew.close()
        self.dialogs.onMessage(self.gui,'lib changed, iPht3D will stop, then restart it')
        self.gui.Destroy()
        
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
