import os
from .config import * # oa 8/6/17
from .geometry import * #OA 9/6/17

class PHT3D:
    """A class that gathers thhe actions for Pht3d, except writing that is left to
    mtphtwriter"""
    def __init__(self, core):
        self.core = core
        self.initBase()
        
    def initBase(self):
        self.Base= {}
        self.Base['Chemistry']={'Solutions':None,'Rates':None,'Phases':None,
                'Exchange':None,'Surface':None,'Kinetic_Minerals':None,'Gases':None}
        self.Base['DualPoro']=[('Immob_porosity',0.),('Transf_coeff',0.)]
        self.Base['Immobile']={'Solutions':None,'Phases':None,'Exchange':None,'Surface':None}
        self.temp={}
        self.temp['Dbase']={}
        self.temp['DBfile']=''
        self.temp['reac']=[0,0,0]  #isoth, ireac, isolver
        self.temp['period']=[0]  # liste de duree et conc de chque periode
        self.temp['run']=0  # Pht3d n'a pas tourne
        self.nsolu=4
        self.listTemps,self.listUserSpecies=[],[' ']
        self.Eactuel,self.Uactuel=0,'mmol/L'
        
    def resetBase(self): self.Base.update(self.core.dicaddin['Chemistry'])
    def getEactuel(self): return self.Eactuel
    def setUserSpecies(self,lEsp): self.listUserSpecies=lEsp
    def setSolutions(self,solu):
        self.Base['Chemistry']['Solutions']=solu;#print solu
        self.nsol=len(solu['cols'])+1
        
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


    def calcNbParm(self):
        """calculates the nb of kinetics parameters, for case the database
        was not imported at that moment"""
        mxp=0;parmk={}
        rates=self.Base['Chemistry']['Rates'];
        npk=len(rates['cols'])-2
        if rates==None: return parmk
        for n in rates['rows']:
            iek=rates['rows'].index(n);mxp=0
            for ip in range(1,npk): # compter le nb de parametres
                if rates['data'][iek][ip+1]!=0:mxp=ip
            parmk[n]=mxp
            if mxp>npk: npk=mxp
        kmin=self.Base['Chemistry']['Kinetic_Minerals']
        if kmin==None: return parmk
        for n in kmin['rows']:
            iek=kmin['rows'].index(n);mxp=0
            for ip in range(1,len(kmin['cols'])): # compter le nb de parametres
                if kmin['data'][iek][ip]!=0:mxp=ip #EV 13/11 replace >0 by !=0
            parmk[n]=mxp
            if mxp>npk: npk=mxp
        return parmk
        
    def updateChemistry(self):
        """update the Chemistry dictionnary with the new database
        elements, rates and phases, Data is a list of lists, one list for
        one colon, the list are ordered according to the rows names"""
        # for chemistry
        old=self.Base['Chemistry'].copy();dic={};#print old['Solutions']
        kw=['Solutions','Rates','Phases','Exchange','Kinetic_Minerals','Surface','Gases']
        dbkw=['SOLUTION_MASTER_SPECIES','RATES','PHASES','EXCHANGE_SPECIES',\
              'RATES','SURFACE_MASTER_SPECIES','GASES']
        nsolu = int(self.core.getValueFromName('Pht3d','NB_SOLU'))+1
        nphase = int(self.core.getValueFromName('Pht3d','NB_PHASE'))
        nexch = int(self.core.getValueFromName('Pht3d','NB_EXCH'))
        nsurf = int(self.core.getValueFromName('Pht3d','NB_SURF'))
        for n in kw: 
            dic[n]={'rows':[],'cols':[],'data':[],'mmol':[]}
        # creating all lists for columns names
        lc1=['C','Backgrd'];
        if old['Solutions']!=None: 
            nsolu=max(nsolu,len(old['Solutions']['cols'])-1)
        for i in range(nsolu-1): 
            lc1.append('solu'+str(i+1))
        lk=['C','IM']
        lk.extend(['parm'+str(i+1) for i in range(self.npk)])
        lk.append('Formula');
        lph=['C','Backg_SI','Backg_Moles']
        for i in range(nphase):
            lph.extend(['Ass'+str(i+1)+'_SI','Ass'+str(i+1)+'_Moles'])
        lex=['C','Backgrd']
        for i in range(nexch): lex.append('Assembl'+str(i+1))
        lsu=['C','Site_back']
        for i in range(nsurf): lsu.append('Sites'+str(i+1))
        lsu.extend(['Specif_area','mass','name','switch'])
        lkp=lk[:-1]
        lkp.remove('IM')
        lcols=[lc1,lk,lph,lex,lkp,lsu,lph]
        # to exclude the kinetic minerals from minerals
        le1=self.tempDbase['PHASES']
        l0=self.tempDbase['RATES']
        l2=[]
        for r in l0:
            if r not in le1: l2.append(r)
        lexclude=[None,le1,None,None,l2,None,None]
        for i in range(len(kw)):
            dic=self.updateDict(dic,old,kw[i],dbkw[i],lcols[i],lexclude[i])
        self.Base['Chemistry'] = dic;#print dic
        
        ##################### for immobile : dual domain poro
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
        lcols=[lc1,lph,lex,lsu]
        for i in range(len(kw)):
            dic=self.updateDict(dic,old,kw[i],dbkw[i],lcols[i],None)
        self.Base['Immobile'] = dic; #print 'Immob',dic
        
    def getChemDict(self,keyword):
        return self.Base['Chemistry'][keyword]
    
    def setChemDict(self,keyword,dict):
        self.Base['Chemistry'][keyword]=dict
        
    def updateDict(self,dic,old,kw,dbkw,lcols,lexclude=None):
        """change the dictionnaries when a new database is imported, to keep
        the values already entered in the dictionnaries"""
        if kw not in list(old.keys()):
            old[kw]=None; return dic
        lrows=list(self.tempDbase[dbkw].keys());lrows.sort();lrows2=[]
        if lexclude!=None:
            for r in lrows:
                if r not in lexclude: lrows2.append(r)
            lrows=lrows2*1
        if len(lrows)!=0:
            dic[kw]['rows']=lrows;dic[kw]['cols']=lcols
            for i in range(len(lrows)):
                if kw=='Rates':
                    a=[False,False];a.extend([0.]*(len(lcols)-2))
                elif kw == 'Surface':
                    a=[False];a.extend([0.]*(len(lcols)-1));a[-2]='';a[-1]=''
                else :
                    a=[False];a.extend([0.]*(len(lcols)-1))
                dic[kw]['data'].append(a)
            if old[kw]!=None:
                for esp in lrows: # get data from current phases
                    if esp not in old[kw]['rows']: continue
                    iold=old[kw]['rows'].index(esp);
                    inew=dic[kw]['rows'].index(esp)
                    for col in lcols:
                        icnew=lcols.index(col);oldcol=old[kw]['cols']
                        if col not in oldcol: continue
                        icold=oldcol.index(col)
                        dic[kw]['data'][inew][icnew]=old[kw]['data'][iold][icold]
        return dic
    
    def getMmol(self,esp):
        lesp1=self.Base['Chemistry']['Solutions']['rows']
        lesp2=self.Base['Chemistry']['Rates']['rows']
        if esp in lesp1:
            iesp=lesp1.index(esp)
            mmol=self.Base['Chemistry']['Solutions']['mmol'][iesp]
        if esp in lesp2:
            iesp=lesp2.index(esp)
            mmol=self.Base['Chemistry']['Rates']['mmol'][iesp]
        else : mmol=0.
        return float(mmol)
            
    def getDictSpecies(self):
        chem=self.Base['Chemistry'];dicE={}
        if chem['Solutions']==None:
            return {'k':[],'i':[],'kim':[],'p':[],'g':[],'e':[],'kp':[],'s':[]}
        #for previous versions
        for kw in ['Kinetic_Minerals','Surface']:
            if kw not in list(chem.keys()):chem[kw]=None
        #rates
        dicE['k']=[];dicE['kim']=[]
        if chem['Rates']!=None:
            rates=chem['Rates'];
            for ir in range(len(rates['rows'])):
                if rates['data'][ir][0]:
                    if rates['data'][ir][1]: # if true immobile kinetics
                        dicE['kim'].append(rates['rows'][ir])
                    else :
                        dicE['k'].append(rates['rows'][ir])
        # solutions
        dicE['i']=[]
        solu=chem['Solutions'];rows=solu['rows']
        self.nsol=len(solu['data'][0])-2;#print self.nsol
        nrow=len(rows);anySol=[False]*self.nsol;#print anySol
        for ir in range(nrow):
            esp = rows[ir]
            if solu['data'][ir][0]:
                if esp not in dicE['kim']:
                    dicE['i'].append(rows[ir])
                for s in range(self.nsol):
                    if float(solu['data'][ir][2+s])!=0.: anySol[s]=True
        #remove species already present in kinetic
        for n in dicE['k']:
            if n in dicE['i']: dicE['i'].remove(n)
        mcomp=len(dicE['i'])+len(dicE['k'])
        # phases,exchange,kinet minerals,surface
        kw=['Kinetic_Minerals','Phases','Exchange','Surface','Gases']
        short=['kp','p','e','s','g']
        for ik in range(len(kw)):
            dicE[short[ik]]=[]
            if kw[ik] not in list(chem.keys()): chem[kw[ik]]=None
            if chem[kw[ik]]!=None:
                d=chem[kw[ik]];
                for ir in range(len(d['rows'])):
                    if (short[ik]=='p') and (d['rows'][ir] in dicE['kp']):
                        continue  # case of phases defined in kinetic minerals
                    if (d['data'][ir][0]): dicE[short[ik]].append(d['rows'][ir]) 
        gcomp = len(dicE['g'])
        # other dissolved species (complexes...)
##        lists=[]
##        if chem.has_key('Species'):
##            for n in chem['Species']['rows']:lists.append(n);
        ncomp=mcomp+len(dicE['kim'])+gcomp
        for k in ['p','e','s','kp']:ncomp+=len(dicE[k])
        dicE['ncomp'],dicE['mcomp'],dicE['gcomp'] = ncomp,mcomp-2,gcomp
        return dicE
        
    def getListSpecies(self):
        dicE, listE = self.getDictSpecies(),[]
        short=['k','i','kim','g','p','e','s','kp']
        for s in short: listE.extend(dicE[s])
        return listE
        
    def importDB(self,fname):
        """to import a phreeqc database and returns a dictionnary in ipht3d format"""
        dicDB={'GASES':None}
        keyw=['SOLUTION_MASTER_SPECIES','SOLUTION_SPECIES','PHASES','RATES',
              'EXCHANGE_MASTER_SPECIES','EXCHANGE_SPECIES','SURFACE_MASTER_SPECIES',
              'SURFACE_SPECIES','GASES','END']
        #klow=['log_k','delta_h','-gamma']
        for k in keyw: dicDB[k]={}
        f1=open(fname)
        curkw='';bufr='';cursp='';parmk={};npk=0;master={};curComplx=''
        for li in f1:
            li2=li.strip() # remove starting blanks
            if li2=='':continue
            if li2[0]=='#':continue
            li2=li2.split('#')[0]; # remove the comments at the end of the line
            if li2 in keyw:
                curkw=li2;continue
            if curkw=='SOLUTION_MASTER_SPECIES': #master species
                a=li.split()
                master[a[0]]=a[1:]
                    
            if curkw=='SOLUTION_SPECIES': # complexes
                li3=li2.split();
                lstk = ['-delta_h','-delta_H','-analytic','-gamma']
                if len(li2.split('='))>1: 
                    b = li2.split('=')
                    if b[0].strip() == b[1].strip() : continue # for the first useless reactions H=H
                    bufr = li2
                    curComplx = li2.split('=')[1].split()[0];#print curComplx # the complex is the species after the = sign
                elif li3[0]=='log_k': # get logk and stores data
                    logk = li3[1]
                    dicDB['SOLUTION_SPECIES'][curComplx] = [bufr,logk] # store name of the phase
                    bufr = ''
                elif li3[0] in lstk: a=0

            if curkw=='PHASES': # minerals and gases
                li3 = li2.split();
                lstk = ['-analytic','delta_h','-Vm','-P_c','-T_c','-Omega']
                if len(li2.split('='))>1: bufr = li2
                elif li3[0]=='log_k': # get logk and stores data
                    logk = li3[1]
                    if '(g)' in curPhase: #gases
                        dicDB['GASES'][curPhase]= [bufr,logk]
                    else : #minerals
                        dicDB['PHASES'][curPhase] = [bufr,logk] # store name of the phase
                    bufr = ''
                elif li3[0] in lstk: a=0
                else: curPhase = li3[0] #the line that contains the name
                        
            if curkw=='RATES': # kinetics
                li3=li2.strip()
                if (li3 in list(master.keys())) or (li3 in dicDB['PHASES']):
                    cursp=li3
                else: bufr=bufr+li
                if li2[:4]=='-end': # store data
                    dicDB['RATES'][cursp]=bufr;mxp=1
                    for i in range(1,20): # counts nb of parametres
                        if bufr.count('parm('+str(i)+')')>0: mxp=i
                    bufr=''
                    parmk[cursp]=mxp
                    if mxp>npk: npk=mxp
    
            if curkw=='EXCHANGE_SPECIES': # exchange species
                a=li.replace('\n','');a=a.split('=');
                if len(a)>1:dicDB[curkw][a[-1]]=0
    
            if curkw=='SURFACE_MASTER_SPECIES': #surface master species
                a=li.split()
                dicDB[curkw][a[0]]=0
                    
        elts=list(master.keys());
        # remove base species if redox are present
        for esp in elts:
            if len(esp.split('('))>1:
                debut=esp.split('(')[0];
                if debut in master:
                    master.pop(debut);
        for n in ['E','O(-2)','Alkalinity','H(0)','H(1)']:#OA 9/6/19  alk and h in list
            if n in master: master.pop(n)
        master['pH']='';master['pe']='';
        dicDB[keyw[0]]=master;
        f1.close(); 
        return dicDB,npk
        
    def readSelectOut(self,fDir):
        '''reads a selected output file made by pht3d with postfix option 
        with user punch (don't let phreeqc write the classical selected output) and get it as arrays 
        (added 8/6/17)
        a difficulty : for large models not all cells are present...
        the selected output file must contain 1st time, 2nd cell_no, then the variables'''
        # get the values
        fname = fDir+os.sep+'selected.out'
        f1=open(fname,'r');a=f1.readline();f1.close()
        lnames = a.split()[2:] #; print(lnames)
        mat = loadtxt(fname,skiprows=1)
        nr,nc = shape(mat)
        #print('mat',shape(mat))
        grd = self.core.addin.getFullGrid()
        nx,ny = grd['nx'],grd['ny'];#print nx,ny
        nlay = getNlayers(self.core)
        ncell = nx*ny*nlay
        times = unique(mat[:,0])
        nt = len(times)
        values = {}
        for n in lnames : values[n] = zeros((nt,nlay,ny,nx))
        for it,t in enumerate(times):
            for i,n in enumerate(lnames):
                a = zeros(ncell)
                indx = mat[mat[:,0]==t,1].astype('int')
                a[indx-1]=mat[mat[:,0]==t,i+2]  ### EV indx start to 1 
                b = reshape(a,(nlay,ny,nx))
                #print('n',b[:,-1::-1,:])
                values[n][it]= b[:,-1::-1,:] #always this inversion...
                #print('val',values)
        return values
    
