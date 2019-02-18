import os
from .config import *
from copy import deepcopy
from .geometry import *

class Min3p:
    """A class that gathers thhe actions for min3p chemistry, except writing that is left to
    min3pwriter"""
    def __init__(self, core):
        self.core = core
        self.initBase()
        self.nodes = None
        
    def initBase(self):
        self.Base= {}
        bs = {'rows':[],'cols':[],'data':[]}
        self.Base['MChemistry']={'comp':bs.copy(),'complex':bs.copy(),'mineral':bs.copy(),
                'sorption':bs.copy(),'exchange':bs.copy(),'redox':bs.copy(),'gases':bs.copy()}
        self.temp={}
        self.temp = {'Dbase':{},'Dbadd':{},'DBfile':''}
        self.listTemps,self.listUserSpecies=[],[' ']
        self.Eactuel,self.Uactuel=0,'mmol/L'
        
    def resetBase(self): self.Base.update(self.core.dicaddin['MChemistry'])
    def getEactuel(self): return self.Eactuel
    def setUserSpecies(self,lEsp): self.listUserSpecies=lEsp
    
    def buildMesh(self,nlay=1,opt='build'):
        '''the objective is to build the mesh and to store most of tis components
        here. With the options read, the mesh is not rebuilt from gmsh, but read from vtk
        reading can seem a little weird, but this is the fastest'''
        self.meshtype = 'grid'
        if self.core.getValueFromName('Min3pFlow','P_Uns')==0: return
        self.core.addin.mesh = self
        self.meshtype='mesh'
        if opt=='build':
            #print 'opt build'
            dct = self.core.diczone['Min3pFlow'].dic
            dicD = dct['spat.6'] # where the domain boundaries are + line .. + points (only pts and line will be used)
            lname,lcoords,lvalue =[],[],[]
            dicM = {'name':[],'coords':[],'value':[]}
            self.polynames = lname
            s = gmeshString(self.core,dicD,dicM)
            os.chdir(self.core.fileDir)
            f1 = open('gm_in.txt','w');f1.write(s);f1.close()
            bindir = self.core.baseDir+os.sep+'bin'+os.sep
            os.system(bindir+'gmsh gm_in.txt -2 -optimize_netgen -o gm_out.msh')
            f1 = open('gm_out.msh','r');s = f1.read();f1.close()
            nodes,elements = readGmshOut(s);#print shape(nodes)
            self.nodes,self.elements,self.nnod,self.nel = nodes,elements,len(nodes),len(elements)
            s= '# vtk DataFile Version 3.0 \n2D scalar \nASCII \n'
            s += 'DATASET UNSTRUCTURED_GRID \nPOINTS '+str(self.nnod)+' float \n'
            if self.core.addin.getDim() in ['Radial','Xsection']: 
                nd = around(nodes[:,[1,3,2]],3)
            else : nd = around(nodes[:,1:],3)
            s += arr2string1(nd)
            s += '\n\nCELLS '+str(self.nel)+' '+str(self.nel*4)+'\n'
            elt2 = elements[:,-4:] # the last four columns
            elt2[:,0] = 3 # this is the nb of points in polygon, here 3 for triangle
            s += self.arr2string2(elt2);#print 'nb pts elts', self.nnod,self.nel
            s += '\nCELL_TYPES '+str(self.nel)+'\n'+'\n'.join(['5']*self.nel)
            f1 = open(self.core.fileDir+os.sep+self.core.fileName+'.vtk','w');
            s = f1.write(s);f1.close()        
        elif opt=='read':
            #print 'opt read'
            f1 = open(self.core.fileDir+os.sep+self.core.fileName+'.vtk','r')
            a = f1.read().replace('\n\n','\n');a =a.replace('\n\n','\n')
            f1.close()
            b = a.split('POINTS')
            c = b[1].split('CELLS');c1=c[0].split('\n')
            nodes = array([a.split() for a in c1[1:-1]]).astype('float')
            nn = shape(nodes)[0]
            nodes = c_[arange(nn),nodes]
            if self.core.addin.getDim() in ['Radial','Xsection']: 
                nodes[:,2] = nodes[:,3]
                nodes[:,3] *= 0
            self.nodes,self.nnod = nodes,nn
            d = c[1].split('CELL_TYPES');d1 = d[0].split('\n')
            elts = array([a.split() for a in d1[1:-1]]).astype('int')
            self.elements, self.nel = elts, shape(elts)[0]
        self = createTriangul(self)
        self.elzones = self.elements[:,1]
        #self.nodzones = self.makeNodeZones(self.elzones)
        
    def getCenters(self):    
        if self.core.getValueFromName('Min3pFlow','P_Uns')==0:
            return getXYmeshCenters(self.core,'Z',0)
        else :
            return self.nodes[:,1],self.nodes[:,2]
            
    def getNumber(self,typ):
        return self.nnod
        
    def arr2string2(self,arr):
        s=''
        nr,nc = shape(arr)
        for i in range(nr):
            s += str(int(arr[i,0]))+' '
            for j in range(1,nc):
                s += str(arr[i,j])+' '
            s += '\n'
        return s
        
    def makeNodeZones(self,elzones):
        '''transforms zones on elements to zones on the nodes corresponding to the elements
        the zone that has the highest number wins!!'''
        nbzones=max(elzones);#print shape(elzones)
        nodzones = zeros(self.nnod)
        for i in range(nbzones+1):
            a0 = unique(ravel(self.trg.triangles[elzones==i,:]));#print a0
            nodzones[a0] = i
        return nodzones
    
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
                    #print 'min3p l 145',k,i,r
                    if Bchem[k]['data'][i][0]: l.append(r)
        else :
            for i,r in enumerate(Bchem[group]['rows']):
                if Bchem[group]['data'][i][0]: l.append(r)
        return l
        
    def importDB(self,oldChem):
        """to import a min3p database and returns a dictionnary in ipht3d format"""
        keyw = ['comp','complex','gases','sorption','exchange','redox','mineral']
        dicDB = {}
        dicAdd = {'sorp2':{},'exch2':{}}
        for k in keyw: dicDB[k]=[]
        localDir = self.core.fileDir
        utilDir = self.core.gui.mainDir+os.sep+'utils'
        # read comp file
        f1=open(self.getCorrectFile(localDir,utilDir,'comp.dbs'))
        for l in f1: 
            if l[:3]=='end': break
            a = l.split()
            k = a[0].replace('\'','')
            dicDB['comp'].append((k,k in oldChem['comp']['rows']))
        f1.close()
        # read gases and complex file
        for n in ['gases','complex']:
            f1=open(self.getCorrectFile(localDir,utilDir,n+'.dbs'))
            while 1>0: 
                l = f1.readline() # reads first line
                if l == None or l[:3]=='end' or l in ['','\n']: break
                a = l.split();
                k = a[0].replace('\'','')
                l = f1.readline()  # OA all folloiwng ofr g and compl largely modified
                if n=='complex':
                    dicDB[n].append([k,l.split()[1::2]])
                else :
                    dicDB[n].append([k,k in oldChem[n]['rows']])
            f1.close()
        # read sorption file
        f1=open(self.getCorrectFile(localDir,utilDir,'sorption.dbs'))
        while 1>0: 
            n='sorption'
            l = f1.readline() # reads first line
            if l == None or l[:3]=='end' or l in ['','\n']: break
            a = l.split();
            k = a[0].replace('\'','');#print k
            l = f1.readline() 
            if '-x' in k: n='exchange'
            dicDB[n].append([k,k in oldChem[n]['rows']])
            dicAdd[n[:4]+'2'][k]=l.split()[1::2]
        f1.close()
        # read redox and mineral phases
        for n in ['redox','mineral']:
            f1=open(self.getCorrectFile(localDir,utilDir,n+'.dbs'))
            for l in f1: 
                if l[0]== '!': continue  # remove commented lines
                l1 = l.replace('\'','').split();#print l1
                kwd = ['transport','irreversible','surface','reversible']
                if len(l1)!=1: continue
                if l1[0] in kwd: continue
                k = l1[0]
                dicDB[n].append((k,k in oldChem[n]['rows']))
            f1.close()
        #sorting database
        for n in keyw:
            #print n,dicDB[n]
            str_list,sec_list = list(zip(*dicDB[n]))
            temp = sorted(zip(str_list, sec_list), key=lambda x: len(x[0]))
            str_list, sec_list = list(map(list, list(zip(*temp))))  
            dicDB[n] = list(zip(str_list, sec_list))
        #print dicDB
        self.temp['Dbase'] = dicDB
        self.temp['Dbadd'] = dicAdd
        
    def getCorrectFile(self,localDir,utilDir,name):
        if name in os.listdir(localDir): return localDir+os.sep+name
        else : return utilDir+os.sep+name

    def updateChemistry(self):
        """when clicking on species names, it creates new species with all data
        and also keep info of the old data already entered
        """
        baseOut = {}
        Db = self.temp['Dbase'].copy();#print Db
        Dbadd = self.temp['Dbadd'].copy()
        lival = {'comp':[True,False,1e-12,1e-12,1e-12,1e-12,'free'],
                'mineral':[True,0.,0.,0.,0.,'.true.','constant',1e-9,1e-12,0.], 
                'exchange':[True,2.0,0.,0.,0.,0.],
                'sorption':[True,1.,1.,1.,10.,6],
                'linear sorption':[False,0,0,0,0],
                'gases':[True,0.,1e-7,2e-5,4.5,300,1.],
                'complex':[True],
                'redox':[True,1e-9],
        }
        linames = {'comp':['C','IM','backgrd','solu1','solu2','solu3','type'],
                'mineral':['C','phimbackg','phim1','phim2','phim3','minequil','update_type','phimin',
                           'keff','sat','rad_init','rad_surf','rad_min'], 
                'exchange':['C','rhob','cecBack','cec1','cec2','cec3'], # OA 24/5 modif to have several exchangers
                'sorption':['C','massBack','mass1','mass2','area','density'],
                'linear sorption':['C','Kd_bak','Kd1','Kd2','Kd3'],
                'gases':['C','1stO dec. rate','Diff_coeff','WKvisco','LJ_sigma','LJ_eK','mass'],
                'complex':['C'],
                'redox':['C','scaling']
        }
        #print Db
        compnames = [name for (name,boo) in Db['comp'] if boo]  # OA 24/5 added to read complexes
        compnames.extend(['h2o','h+','oh-'])
        lins = 'linear sorption'
        if lins not in list(Db.keys()): 
            Db[lins]=deepcopy(Db['comp'])
            self.Base['MChemistry'][lins]={'rows':Db['comp']}
        #print Dbadd
        for n in list(Db.keys()):
            base = {'rows':[],'cols':[],'data':[],'text':[]}
            bIn = deepcopy(self.Base['MChemistry'][n])
            rowtrue = [name for (name,boo) in Db[n] if boo];#print rowtrue
            if n == 'exchange':
                lst1 = [name for (name,boo) in Db[n] if boo]
                rowtrue = ['-x'] # adds '-x' at the begnnning
                rowtrue.extend(lst1);#print rowtrue
                base['text']=lst1*1
            if n == 'sorption':
                lst1 = [name for (name,boo) in Db[n] if boo];#print lst1
                rowtrue = []
                for k in lst1: 
                    basespec = Dbadd['sorp2'][k][1]
                    if basespec not in rowtrue: rowtrue.append(basespec)
                rowtrue.extend(lst1);#print rowtrue
                base['text']=lst1*1
            if n=='complex': 
                rowtrue = [name for (name,cpl) in Db[n] if self.useComplex(compnames,cpl)]
            base['cols'] = linames[n]
            for name in rowtrue:
                base['rows'].append(name)
                if name in bIn['rows']:
                    line = bIn['data'][bIn['rows'].index(name)]
                    base['data'].append(line)
                else :
                    val = lival[n]*1
                    if name in base['text']: val[0]=False
                    base['data'].append(val)
                if n=='gases' and name in['o2(g)','n2(g)','co2(g)','ch4(g)']:
                   base['data'][-1][1] = 'LJ'
                   base['data'][-1][2:] = self.getLJparms(name)
            baseOut[n] = base.copy()
        #print baseOut['complex']
        self.Base['MChemistry'] = deepcopy(baseOut)
        
    def getBase(self): 
        #print 'min3p 250',self.core.getValueFromName('Min3pTrans','Diff_choice')
        b1 = self.Base['MChemistry'].copy()
        if self.core.getValueFromName('Min3pTrans','Diff_choice')==0:
            #used to modify the view according to the diffusion type
            b1['gases']['cols']=['C','1stO dec. rate']
            data = b1['gases']['data']
            b1['gases']['data'] = [d[:2] for d in data]
        # remove secondary species for exchange and sorption
        b1['exchange']['rows'] = ['-x']
        lst0 = b1['sorption']['rows']*1
        for n in lst0 : 
            if n in b1['sorption']['text']: lst0.remove(n)
        b1['sorption']['rows'] = lst0
        return b1
    
    def useComplex(self,complist,complexIn): # OA 24/5
        '''finds if all compounds in a complex are in the list of the current compounds'''
        test = True
        for n in complexIn : 
            if n not in complist: test = False; break
        return test
            
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

        
