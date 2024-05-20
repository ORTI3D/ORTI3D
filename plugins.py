# -*- coding: utf-8 -*-
"""
Created on Thu Jan  4 09:19:11 2024

@author: olivi
"""
from config import *

class plugins:
    
    def __init__(self,core):
        self.core = core
        self.pl_list=['Coupling','Foam']#,'Immobile oil'] #'sorptionAW'
        for pl_name in self.pl_list:
            if pl_name not in core.dicplugins.keys():
                core.dicplugins[pl_name]={'active':False,'data':[]}
    
    def setGui(self,gui):
        self.gui= gui
        cfg = Config(self.core)
        if gui != None:
            self.dialogs = cfg.dialogs
                    
    def action(self):
        '''this is the action when the plugin is called (generally a dialog)'''
        pl_name = self.gui.menuBar().sender().text()[1:];
                
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

        if pl_name=='Radon':
            act=self.core.dicplugins[pl_name]['active']
            d = self.core.dicplugins[pl_name]['data']
            data = [('Use gas decay?','Check',act),
                    ('gas decay rate (s-1)','Text',d)] 
            dialg = self.dialogs.genericDialog(self.gui,'Radon',data)
            retour = dialg.getValues()
            if retour != None:
                self.core.dicplugins[pl_name]['active']=retour[0]
                self.core.dicplugins[pl_name]['data']=float(retour[1])
                
        if pl_name=='Coupling':
            act=0
            if pl_name in self.core.dicplugins.keys():
                act=self.core.dicplugins[pl_name]['active']
            else : 
                self.core.dicplugins[pl_name]={'active':0,'data':[]}
            #first dialog
            quest = [('Use coupling?','Check',act)] 
            dialg = self.dialogs.genericDialog(self.gui,'Coupling',quest)
            retour = dialg.getValues()
            if retour != None:
                self.core.dicplugins[pl_name]['active']=retour[0]
            dic = {'hEqn':{}}
            dic['hEqn']['rows']=['','','']
            dic['hEqn']['cols']=['Chk','Ykey','Xvar','Xref','typ','nparms','a0','a1','a2']
            dic['hEqn']['data']=[[False,'muw','T',25,'linear',1,-1e-3],
                    [False,'muw','C',0,'linear',1,1e-2],
                    [False,'K','eps',0,'kozeny',0,0]]
            #print("dic plug ",dic)
            dics = self.core.dicplugins[pl_name]['data']
            if len(dics)>0:
                for k in dic.keys():
                    if k in dics.keys(): 
                        for i,l in enumerate(dics[k]['data']):
                            dic[k]['data'][i]=l
            dialg = self.dialogs.myNoteBook(self.gui,"coupling",dic)
            retour = dialg.getValues()
            if retour != None:
                self.core.dicplugins[pl_name]['data']=retour  
                
        if pl_name=='Foam':
        # foam injection
            if pl_name in self.core.dicplugins.keys():
                act=self.core.dicplugins[pl_name]['active']
                d = self.core.dicplugins[pl_name]['data']
            else : 
                d=[1e4,1000,0.35,1e-4,1e4];act=0
                self.core.dicplugins[pl_name]={'active':0,'data':d}
            data = [('Use Foam?','Check',act),
                    ('Fmmob','Text',d[0]),('epdry','Text',d[1]),
                    ('Fmdry','Text',d[2]),('Cref','Text',d[3]),
                    ('Fc','Text',d[4])] 
            dialg = self.dialogs.genericDialog(self.gui,'Foam',data)
            retour = dialg.getValues()
            if retour != None:
                self.core.dicplugins[pl_name]['active']=retour[0]
                self.core.dicplugins[pl_name]['data']=retour[1:]
            #print(self.core.dicplugins)
            
    def writer(self,pl_name):
        '''this is called when opfoam writes its files'''
        if pl_name=='sorptionAW':
            s='\n activateSorptionAW 1; \n'
            self.core.opfWriter.addTransportPropeties(s)
            
        if pl_name=='immobile':
            s='\n'.join(self.core.dicplugins[pl_name]['species'])
            f1=open(self.core.fileDir+'constant\\options\\immobile','w')
            f1.write(s);f1.close()

        if pl_name=='Radon':
            # get gas names and write the file
            gspc=self.core.addin.chem.getDictSpecies()['g'];print('writing Rn')
            s=str(len(gspc))+' 1\n'
            for g in gspc:
                s+=g
                if g=='Rn(g)': s+=' '+str(self.core.dicplugins['Radon']['data'])+'\n'
                else : s+=' 0\n'
            f1=open(self.core.fileDir+'constant\\options\\gasDecay','w')
            f1.write(s);f1.close()
            
        if pl_name=='Coupling':
            if pl_name not in self.core.dicplugins.keys(): return
            if self.core.dicplugins[pl_name]['active']==0: return
            lDic = self.core.dicplugins[pl_name]['data']
            s = ''
            for k in lDic.keys():
                dic=lDic[k];print(dic)
                for d in dic['data']:
                    if d[0]: #d[0] is to check state
                        s +=' '.join(str(a) for a in d[1:])+'\n'
            f1=open(self.core.fileDir+'constant\\options\\coupling','w')
            f1.write(s);f1.close()
            
        if pl_name=='Foam':
            s='\n'.join(self.core.dicplugins[pl_name]['data'])
            f1=open(self.core.fileDir+'constant\\options\\foam','w')
            f1.write(s);f1.close()
       