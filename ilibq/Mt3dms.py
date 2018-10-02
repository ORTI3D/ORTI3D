import os

class MT3DMS:
    """A class that gathers the actions for Mt3d, except writing that is left to
    mtphtwriter"""
    def __init__(self, core):
        self.core = core
        self.initBase()
        
    def initBase(self):
        self.Base= {}
        self.Base['Species']={'Solutions':None}
        #self.Base['DualPoro']=[('Immob_porosity',0.),('Transf_coeff',0.)]
        self.Base['Immobile']={'Solutions':None}
        self.nsol = 4
        
    def resetBase(self): 
        self.Base.update(self.core.dicaddin['Mt3dData'])
    def setSolutions(self,solu):
        self.Base['Species']['Solutions']=solu;#print solu
        self.nsol=len(solu['cols'])+1
    def setBaseFromNames(self,namelist):
        dic={}
        dic['rows']=namelist
        
    def setImportedSolutions(self,dicData):
        rows,dataIn=dicData['rows'],dicData['data']
        nl,nc=shape(dataIn);data=[]
        for i in range(nl):
            a=[True];a.extend(list(dataIn[i,:]));data.append(a)
        cols=['C','Backgrd']
        cols.extend(['solu'+str(i+1) for i in range(len(data[0])-2)]);
        dic={'cols':cols,'data':data,'rows':rows};#print dic
        #print len(cols),shape(data),data
        self.setSolutions(dic)
        
    def updateChemistry(self):
        """update the Chemistry dictionnary from the database
        elements, rates and phases, Data is a list of lists, one list for
        one colon, the list are ordered according to the rows names"""
        # for chemistry
        old=self.Base['Chemistry'].copy();dic={};#print old['Solutions']
        kw=['Solutions','Rates','Phases','Exchange','Kinetic_Minerals','Surface']
        dbkw=['SOLUTION_MASTER_SPECIES','RATES','PHASES','EXCHANGE_SPECIES',\
              'RATES','SURFACE_MASTER_SPECIES']
        for n in kw: 
            dic[n]={'rows':[],'cols':[],'data':[],'mmol':[]}
        # update all dictionnaries
        lc1=['C','Backgrd'];
        nsolu=self.nsolu
        if old['Solutions']!=None: 
            nsolu=max(nsolu,len(old['Solutions']['cols'])-1)
        for i in range(nsolu-1): 
            lc1.append('solu'+str(i+1))
        lcolk=['C','IM']
        lcolk.extend(['parm'+str(i+1) for i in range(self.npk)])
        lcolk.append('Formula');
        lc2=['C','Backg_SI','Backg_Moles','Ass1_SI','Ass1_Moles','Ass2_SI','Ass2_Moles',\
            'Ass3_SI','Ass3_Moles','Ass4_SI','Ass4_Moles']
        lc3=['C','Backgrd','Assembl1','Assembl2','Assembl3','Assembl4','Assembl5']
        lc4=['C','Site_back','Sites1','Sites2','Sites3','Sites4','Specif_area','mass','name','switch']
        lckp=lcolk[:-1]
        lckp.remove('IM')
        lcols=[lc1,lcolk,lc2,lc3,lckp,lc4]
        le1=self.tempDbase['PHASES']
        l0=self.tempDbase['RATES']
        l2=[]
        for r in l0:
            if r not in le1: l2.append(r)
        lexclude=[None,le1,None,None,l2,None]
        for i in range(len(kw)):
            dic=self.updateDict(dic,old,kw[i],dbkw[i],lcols[i],lexclude[i])
        self.Base['Chemistry'] = dic;#print dic
        
        # for immobile : dual domain poro
        if 'Immobile' in self.Base: old=self.Base['Immobile'].copy()
        else : old={'Solutions':None,'Phases':None,'Exchange':None,'Surface':None}
        dic={};
        # add poro and exch to solutions
        self.tempDbase['SOLUTION_MASTER_SPECIES']['Imm_poro']=''
        self.tempDbase['SOLUTION_MASTER_SPECIES']['Transf_coeff']=''
        kw=['Solutions','Phases','Exchange','Surface']
        for n in kw:
            dic[n]={'rows':[],'cols':[],'data':[],'mmol':[]}
        # update all dictionnaries
        dbkw=['SOLUTION_MASTER_SPECIES','PHASES','EXCHANGE_SPECIES','SURFACE_MASTER_SPECIES']
        lcols=[lc1,lc2,lc3,lc4]
        for i in range(len(kw)):
            dic=self.updateDict(dic,old,kw[i],dbkw[i],lcols[i],None)
        self.Base['Immobile'] = dic; #print 'Immob',dic
        
    def getChemDict(self,keyword):
        return self.Base['Chemistry'][keyword]
    
    def setChemDict(self,keyword,dict):
        self.Base['Chemistry'][keyword]=dict
        
        
    def getListSpecies(self):
        dicE, listE = self.getDictSpecies(),[]
        short=['k','i','kim','p','e','kp','s']
        for s in short: listE.extend(dicE[s])
        return listE
