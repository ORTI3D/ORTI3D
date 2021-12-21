# -*- coding: utf-8 -*-
"""
Created on Sun Feb 16 09:47:03 2014

@author: olive
"""
from .valueDialog import valueDialog
from .config import *
import os

class BaseParms:
    def __init__(self,gui,core):
        cfg = Config(core)
        self.gtyp = cfg.gtyp
#        if self.gtyp=='wx': 
#            self.gui,self.visu = gui,visu
#        elif self.gtyp=='qt': 
#            self.gui,self.visu = gui.gui,gui.visu
        self.dialogs = cfg.dialogs
        self.core,self.addin,self.gui,self.visu = core,core.addin,gui,gui.visu
        
        self.groups={'1.Model':['Map'],
            '2.Flow':['Parameters','Write','Run'],
            '3.Transport':['Parameters','Write','Run'],
            '4.Chemistry':['Parameters','Write','Run'],
            '5.Pest':['Parameters','Write','Run']
            }
        self.tipNames={}
        self.dicaction={
            'Mo_Map':'self.openMap()',
            'Fl_Parameters':'self.getParameters(\'Modflow\')',
            'Fl_Write': 'self.dialogs.onMessage(self.gui,self.writeModel(\'Modflow\'))',
            'Fl_Run': 'self.runModel(\'Modflow\')',
            'Tr_Parameters':'self.getParameters(\'Mt3dms\')',
            'Tr_Write': 'self.dialogs.onMessage(self.gui,self.writeModel(\'Mt3dms\'))',
            'Tr_Run': 'self.runModel(\'Mt3dms\')',
            'Ch_Parameters':'self.getParameters(\'Pht3d\')',
            'Ch_Write': 'self.dialogs.onMessage(self.gui,self.writeModel(\'Pht3d\'))',
            'Ch_Run': 'self.runModel(\'Pht3d\')',
            'Pe_Parameters':'self.getParameters(\'Pest\')',
            'Pe_Write': 'self.dialogs.onMessage(self.gui,self.writeModel(\'Pest\'))',
            'Pe_Run': 'self.runModel(\'Pest\')',
            }
    
    def action(self,name):
        action=self.dicaction[str(name)]
        #print (self.core.dicaddin['Model']['group'],action)
        mgroup = self.core.dicaddin['Model']['group']
        if 'USG' in mgroup: # only the transport is different in USG
            action=action.replace('Mt3dms','MfUsgTrans')            
            action=action.replace('Modflow','Modflow_USG')         # OA 17/2/20
        if mgroup =='Min3p':
            action=action.replace('Modflow','Min3pFlow')
            action=action.replace('Mt3dms','Min3pTrans')
            action=action.replace('Pht3d','Min3pChem')
        if mgroup =='Opgeo':
            action=action.replace('Modflow','OpgeoFlow')
            action=action.replace('Mt3dms','OpgeoTrans')
            action=action.replace('Pht3d','OpgeoChem')
        if mgroup =='Sutra':
            action=action.replace('Modflow','Sutra')
            action=action.replace('Mt3dms','Sutra')
        exec(action)
        
    def openMap(self):
        dlg = self.dialogs.myFileDialog('Open')
        fDir,fName = dlg.getsetFile(self.gui,'Choose map',"*.png")
        if fDir == '': return
        if fDir != None:          
            self.gui.map = str(fDir)+os.sep+str(fName)+'.png' 
            self.gui.visu.createMap()
            self.gui.guiShow.dlgShow.onTickBox('Model','Map','B',True)
        #else : return

    def getParameters(self,modName):
        """this method opens a valuedialog in a frame and enters all parameters"""
        #print modName
        self.core.updateDicts()
        dlg = valueDialog(self.gui,'parameters',self.core,modName)
        retour = dlg.show()
        
    def writeModel(self,modName):
        #print modName
        if self.gtyp=='qgis':
            self.gui.visu.zonesQgs2core()
        messg = self.core.writeModel(modName)
        a = self.core.makeTtable()
        tl2 = self.core.getTlist2()
        listSpec = self.core.addin.chem.getListSpecies() # just the names
        self.gui.guiShow.init() # OA 22/8/19 added
        self.gui.guiShow.setChemSpecies(listSpec)
        #self.gui.guiShow.resetDicContour() # the show panel nows that everything has to be reread # OA modifs 1/10/19
        #self.gui.guiShow.setNames('Model_Tstep_L',tl2) #EV 16/12/21
        return messg
        
    def runModel(self,modName):
        message = self.core.runModel(modName);#print message
        self.dialogs.onMessage(self.gui,message)
        self.gui.guiShow.redraw() #EV 16/12/21
        #if modName == 'Pht3d': self.core.usePostfix() #OA 3/10/19 removed, usePostFix does not exist anymore?
