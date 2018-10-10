import os
from .config import *
from copy import deepcopy

class Min3pChem:
    """A class that gathers thhe actions for min3p chemistry, except writing that is left to
    min3pwriter"""
    def __init__(self, core):
        self.core = core
        self.initBase()
        
    def initBase(self):
        self.Base= {}
        bs = {'rows':[],'cols':[],'data':[]}
        self.Base['MChemistry']={'comp':bs.copy(),'complex':bs.copy(),'mineral':bs.copy(),
                'sorption':bs.copy(),'exchange':bs.copy(),'redox':bs.copy(),'gases':bs.copy()}
        self.temp={}
        self.temp['Dbase']={}
        self.temp['DBfile']=''
        self.listTemps,self.listUserSpecies=[],[' ']
        self.Eactuel,self.Uactuel=0,'mmol/L'
        
    def resetBase(self): self.Base.update(self.core.dicaddin['MChemistry'])
    def getEactuel(self): return self.Eactuel
    def setUserSpecies(self,lEsp): self.listUserSpecies=lEsp
        
    def getMmol(self,esp):
        lesp1=self.Base['MChemistry']['Solutions']['rows']
        lesp2=self.Base['MChemistry']['Rates']['rows']
        if esp in lesp1:
            iesp=lesp1.index(esp)
            mmol=self.Base['MChemistry']['Solutions']['mmol'][iesp]
        if esp in lesp2:
            iesp=lesp2.index(esp)
            mmol=self.Base['MChemistry']['Rates']['mmol'][iesp]
        else : mmol=0.
        return float(mmol)
        
    def getListSpecies(self,group=None):
        """returns the list of species, if nothing specified returns all, if
        group specified returns only the species of the group"""
        l = []
        Bchem = self.Base['MChemistry']
        if group==None:
            for k in list(Bchem.keys()):
                for i,r in enumerate(Bchem[k]['rows']):
                    if Bchem[k]['data'][i][0]: l.append(r)
        else :
            for i,r in enumerate(Bchem[group]['rows']):
                if Bchem[group]['data'][i][0]: l.append(r)
        return l
        
    def importDB(self,oldChem):
        """to import a min3p database and returns a dictionnary in ipht3d format"""
        keyw = ['comp','complex','gases','sorption','exchange','redox','mineral']
        dicDB = {}
        for k in keyw: dicDB[k]=[]
        utilDir = self.core.gui.mainDir+os.sep+'utils'
        # read comp file
        f1=open(utilDir+os.sep+'comp.dbs')
        for l in f1: 
            if l[:3]=='end': break
            a = l.split()
            k = a[0].replace('\'','')
            dicDB['comp'].append((k,k in oldChem['comp']['rows']))
        f1.close()
        # read gases sorption and complex file
        for n in ['gases','complex','sorption']:
            f1=open(utilDir+os.sep+n+'.dbs')
            while 1>0: 
                l = f1.readline() # reads first line
                if l == None or l[:3]=='end' or l in ['','\n']: break
                a = l.split();
                k = a[0].replace('\'','')
                if '-x' in k:
                    dicDB['exchange'].append((k,k in oldChem[n]['rows']))
                else :
                    dicDB[n].append((k,k in oldChem[n]['rows']))
                l = f1.readline() # 2nd line is not used
            f1.close()
        # read redox and mineral phases
        for n in ['redox','mineral']:
            f1=open(utilDir+os.sep+n+'.dbs')
            for l in f1: 
                if l[0]== '!': continue  # remove commented lines
                l1 = l.replace('\'','').split();#print l1
                kwd = ['transport','irreversible','surface']
                if len(l1)!=1: continue
                if l1[0] in kwd: continue
                k = l1[0]
                dicDB[n].append((k,k in oldChem[n]['rows']))
            f1.close()
        #sorting database
        for n in keyw:
            str_list,sec_list = list(zip(*dicDB[n]))
            temp = sorted(zip(str_list, sec_list), key=lambda x: len(x[0]))
            str_list, sec_list = list(map(list, list(zip(*temp))))  
            dicDB[n] = list(zip(str_list, sec_list))
        #print dicDB
        self.temp['Dbase'] = dicDB

    def updateChemistry(self):
        """when clicking on species names, it creates new species with all data
        and also keep info of the old data already entered
        """
        baseOut = {}
        Db = self.temp['Dbase'].copy()
        lival = {'comp':[True,1e-12,1e-12,1e-12,1e-12,'free'],
                'mineral':[True,0.,0.,0.,0.,'.true.','constant',1e-9,1e-12,0.], 
                'exchange':[True,0.,1.87],
                'sorption':[True,'=soh',1.,1.,1.,10.,6],
                'linear sorption':[False,0,0,0,0],
                'gases':[True,1e-7,2e-5,4.5,300,1.],
                'complex':[True],
                'redox':[True,1e-9],
        }
        linames = {'comp':['C','backgrd','solu1','solu2','solu3','type'],
                'mineral':['C','phibackg','phim1','phim2','phim3','minequil','update_type','phimin','keff','sat'], 
                'exchange':['C','cec','rhobulk'],
                'sorption':['C','surface name','mass0','mass1','mass2','area','density'],
                'linear sorption':['C','Kd_bak','Kd1','Kd2','Kd3'],
                'gases':['C','Diff_coeff','WKvisco','LJ_sigma','LJ_eK','mass'],'complex':['C'],
                'redox':['C','scaling']
        }
        #print Db
        lins = 'linear sorption'
        if lins not in list(Db.keys()): 
            Db[lins]=deepcopy(Db['comp'])
            nrows = len(Db['comp'])
            self.Base['MChemistry'][lins]={'rows':Db['comp']}
        for n in list(Db.keys()):
            base = {'rows':[],'cols':[],'data':[],'text':[]}
            bIn = deepcopy(self.Base['MChemistry'][n])
            rowtrue = [name for (name,boo) in Db[n] if boo];#print rowtrue
            if n=='comp': rowtrue.extend(['ph','pe'])
            base['cols'] = linames[n]
            for name in rowtrue:
                base['rows'].append(name)
                if name in bIn['rows']:
                    line = bIn['data'][bIn['rows'].index(name)]
                    base['data'].append(line)
                else :
                    val = lival[n]*1
                    if name in ['ph','pe']: val[-1]=name
                    base['data'].append(val)
                if n=='gases' and name in['o2(g)','n2(g)','co2(g)','ch4(g)']:
                   base['data'][-1][1] = 'LJ'
                   base['data'][-1][2:] = self.getLJparms(name)
            baseOut[n] = base.copy()
        self.Base['MChemistry'] = deepcopy(baseOut)
        
    def readKinetics(self,lspecies):
        """the format of kinetic data : {'formula':'','type':['equilibrium','irreversible']
        'hyperbolic T^a':'','inhibition T^a':'','inhibition phi^m':''}
        """
        #reads the kinetics of the given species in redox 
        utilDir = self.core.gui.mainDir+os.sep+'utils'
        f1=open(utilDir+os.sep+'redox.dbs')
        l = f1.readline()
        self.Base['MChemistry']['redox']['text']=['']*len(lspecies)
        while 'end of database' not in l: 
            l = f1.readline()
            sp = ' '
            if len(l)>2: sp = l.replace('\'','').split()[0];#print sp
            if sp in lspecies:
                form = f1.readline()
                typ = f1.readline()
                l = f1.readline()
                while len(l.rstrip())== 0: 
                    l = f1.readline();# read following blank lines
                hypT = l
                while len(l.rstrip())>0: # hyperbolic
                    l = f1.readline(); hypT += l
                l = f1.readline()
                s, inhT = '',l
                if l[1:4]=='inh': # inhibiton by solutes
                    while len(l.rstrip())>0:
                        l = f1.readline();inhT += l
                l = f1.readline()
                s,inhP = '',l
                if l[1:4]=='inh': # inhibition by minerals
                    while len(l.rstrip())>0:
                        l = f1.readline();inhP += l
                txt = form+typ+hypT+inhT+inhP
                indx = lspecies.index(sp)
                self.Base['MChemistry']['redox']['text'][indx]=txt
        #print self.Base['MChemistry']['redox']['text']
        f1.close()
        
    def readMinerals(self,lspecies):
        """the format of minerals data
        """
        #reads the kinetics of the given species in redox 
        utilDir = self.core.gui.mainDir+os.sep+'utils'
        f1=open(utilDir+os.sep+'mineral.dbs')
        l = f1.readline()
        self.Base['MChemistry']['mineral']['text']=['']*len(lspecies)
        while 'end of database' not in l: 
            l = f1.readline()
            sp = ' '
            if len(l)>2: sp = l.replace('\'','').split()[0];#print 'mineral',sp
            if sp in lspecies:
                typ1 = f1.readline() # classical : 'surface'
                gfw = f1.readline()  # weigth and density
                form = f1.readline()# stochiometry
                typ2 = f1.readline() # type : reversible irreversible...
                l = f1.readline()
                while len(l.rstrip())== 0: 
                    l = f1.readline();# read following blank lines
                s = ''
                while l[1:4] in ['fra','hyp','inh','log']: # read different keyw
                    s += l
                    while len(l.rstrip())>0:
                        l = f1.readline();s += l # read constants
                    while len(l.rstrip())== 0: 
                        l = f1.readline();# read following blank line 
                txt = typ1+gfw+form+typ2+s
                indx = lspecies.index(sp)
                self.Base['MChemistry']['mineral']['text'][indx]=txt
        f1.close()
                
    def readOtherDb(self,name,lspecies):
        """read the other databases files comp, complex, gases to retain
        only the species of interest"""
        utilDir = self.core.gui.mainDir+os.sep+'utils'
        f1=open(utilDir+os.sep+name+'.dbs')
        l = f1.readline()
        while l[:3] != 'end':
            sp = l.split()[0]
            if sp in lspecies : 
                txt = l
                if name in ['gases','complex']:
                    f1.readline() # to read the second line
                self.Base['MChemistry'][name]['text'].append(txt)
            l = f1.readline()
        f1.close()

    def writeDb(self,name):
        """returns the correct string of kinetics to be printed in the database"""
        dicKin = self.Base['MChemistry'][name];#print name
        s = ''
        for i,ek in enumerate(dicKin['rows']):
            if name in ['redox','mineral']:
                s += '!\n\n\''+ek +'\''+ '\n' # name of the species                
            if ek in ['ph','pe']: continue
            if 'text' in dicKin: s += dicKin['text'][i]
        if name=='mineral': return s + '!\n! end of database\n!\n\'end\''
        else : return s + '\n!end'
        
    def getLJparms(self,name):
        '''wilke visco, LJ sigma and LJ e/K'''
        d = {'n2(g)':[1.923e-5,3.667,99.8,28.01],'o2(g)':[1.952e-5,3.433,113,32],
             'co2(g)':[1.442e-5,3.996,190,41.01],'ch4(g)':[2.056e-5,3.78,154,16.04]}
        return d[name]

        