# -*- coding: utf-8 -*-
"""
Created on Thu Jan  4 09:19:11 2024

@author: olivi
"""
from .config import *

class plugins:
    
    def __init__(self,core):
        self.core = core
        self.pl_list=['Radon','immobile oil'] #'sorptionAW'
        for pl_name in self.pl_list:
            if pl_name not in core.dicplugins.keys():
                core.dicplugins[pl_name]={'active':False,'data':None}
    
    def setGui(self,gui):
        self.gui= gui
        cfg = Config(self.core)
        if gui != None:
            self.dialogs = cfg.dialogs
        
    def action(self):
        '''this is the action when the plugin is called (generally a dialog)'''
        pl_name = self.gui.menuBar().sender().text()[1:];
        if pl_name=='sorptionAW':
            act=self.core.dicplugins[pl_name]['active']
            data = [('Use Air-Water sorption?','Check',act)] 
            dialg = self.dialogs.genericDialog(self.gui,'AWsorption',data)
            retour = dialg.getValues()
            if retour != None:
                self.core.dicplugins[pl_name]['active']=True
                dk={'kw':'OISAW','detail':['Air-water sorption','no','active'],'type':'choice','default':0}
                self.core.dickword['OpenTrans'].addKeyword('rct.1',dk)
                self.core.addin.setMtSpeciesList(['Koc','RC1','RC2','AW_a','AW_b'])
                
        if pl_name=='immobile':
            act=self.core.dicplugins[pl_name]['active']
            lspc=self.core.dicplugins[pl_name]['species']
            data = [('Use immobile species?','Check',act),
                    ('Species','Textlong','\n'.join(lspc))] 
            dialg = self.dialogs.genericDialog(self.gui,'Immobile compnts',data)
            retour = dialg.getValues()
            if retour != None:
                self.core.dicplugins[pl_name]['active']=retour[0]
                self.core.dicplugins[pl_name]['species']=split(retour[1])
            
    
    def writer(self,pl_name):
        '''this is called when opfoam writes its files'''
        if pl_name=='sorptionAW':
            s='\n activateSorptionAW 1; \n'
            self.core.opfWriter.addTransportPropeties(s)
            
        if pl_name=='immobile':
            s='\n'.join(self.core.dicplugins[pl_name]['species'])
            f1=open(self.core.fileDir+'constant\\options\\'+immobile,'w')
            f1.write(s);f1.close()
        