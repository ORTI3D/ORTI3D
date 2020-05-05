# -*- coding: cp1252 -*-
from array import array as arr2
import os,time
import numpy as np
from .mtUsgKeywords import Mtu
from .geometry import *
from .modflowWriter import * # OA 6/5/19

class mtUsgWriter:

    def __init__(self, core,fDir, fName):
        self.core,self.Mkey = core,Mtu()
        self.fDir,self.fName = fDir,fName
        self.fullPath = fDir+os.sep+fName;#print self.fullPath
        self.mfloW = modflowWriter(core,fDir,fName) # OA 6/5/19

    def writeMtphtFiles(self,listEsp,opt,parmk=None):
        self.mgroup = self.core.dicaddin['Model']['group']
        self.listEsp = listEsp
        usgTrans = {'active':True}
        self.ttable = self.core.makeTtable();#print 'mtpht ttable',self.ttable
        self.dim = self.core.addin.getDim()
        self.core.updateDicts()
        tlist = array(self.ttable['tlist'])
        self.per = tlist[1:]-tlist[:-1]
        self.nper = len(self.per)#;print('writempht l.26',self.nper)
        usgTrans['nam'] = self.writeNamString(opt)
        if amax(self.core.getValueLong('MfUsgTrans','crch.1',0))>0: 
            usgTrans['rch'] = self.writeRchString()
        if 'cwell.1' in self.core.diczone['MfUsgTrans'].dic.keys():
            usgTrans['wel'] = self.writeWelValues(opt)            
        if 'bct.20' in self.core.diczone['MfUsgTrans'].dic.keys():
            usgTrans['chd'] = self.writeChdValues(opt)            
        self.mfloW.writeModflowFiles(self.core,usgTrans=usgTrans)
        self.writeBCT(opt)
        if opt=='Pht3d':
            self.writePhFile(self.core,listEsp,parmk)
            #self.writePhreeqc(self.core,listEsp);
        if 'pcb.2' in list(self.core.diczone['MfUsgTrans'].dic.keys()):
            if self.mgroup == 'Modflow USG_rect': self.writePcbRectFile(self.core)
            else :self.writePcbFile(self.core)
        return 
        
    def writeNamString(self,opt):
        s='BCT  41 '+self.fName+'.bct\n'
        if 'pcb.2' in list(self.core.diczone['MfUsgTrans'].dic.keys()):
            s += 'PCB  42 '+self.fName+'.pcb\n'
        if opt=='Pht3d':
            s += ' PHT  64    Pht3d_ph.dat\n'
        s += 'DATA(BINARY) 101  '+self.fName+'.conc\n'
        s += 'DATA(BINARY) 102  '+self.fName+'.cbb\n'
        return s
        
    #*********************** BCT file writer ****************
    def writeBCT(self,opt):
        """to write modflow usg BCT transport file.
        reads the keyword file and prints all keywords by types : param (0D)
        vector (1D) array (2D). types are found by (dim1,dim2).."""
        lexceptions, s =['bct.5'],''
        llist=self.Mkey.groups['BCT'];#print n1,name
        if opt=='Pht3d':
            self.core.setValueFromName('MfUsgTrans','MCOMP',self.listEsp['mcomp']+3)
            nimcomp = self.listEsp['ncomp']-self.listEsp['mcomp']
            self.core.setValueFromName('MfUsgTrans','NIMCOMP',nimcomp)
        else :
            self.core.setValueFromName('MfUsgTrans','MCOMP',1)
            self.core.setValueFromName('MfUsgTrans','NIMCOMP',0)
        for ll in llist:
            cond=self.Mkey.lines[ll]['cond'];#print('mtw 50',ll)
            if self.testCondition(cond)==False : continue
            kwlist=self.Mkey.lines[ll]['kw']
            ktyp=self.Mkey.lines[ll]['type']
            lval=self.core.dicval['MfUsgTrans'][ll];
            if ll in lexceptions:
                s += self.writeExceptions(ll);continue
            if ktyp[0] in ['vecint','vecfloat','arrint','arrfloat']:
                s += self.writeArray(opt,ll,ktyp[0]) + '\n'
            elif ktyp[0]=='title': # case of a title line
                s += '#'+str(ktyp[0]).rjust(10)+'\n'
            else : # classical keywords
                for ik in range(len(kwlist)):
                    if ik<len(lval): 
                        if type(lval[ik])==type(5):s += ' %3i ' %lval[ik] # OA 25/4/20
                        elif type(lval[ik])==type(5.0):s += ' %9.3e ' %lval[ik]
                        else :s +=  lval[ik]
                    else : s += '0'.rjust(10)
                s += '\n'
        f1=open(self.fDir+os.sep+self.fName +'.bct','w')
        f1.write(s);f1.close()

    def testCondition(self,cond):
        """ test if the condition is satisfied"""
        return self.core.testCondition('MfUsgTrans',cond)
        
    def writeExceptions(self,line):
        if line == 'bct.5': 
            if self.mgroup == 'Modflow USG_rect': s=''
            else : s = 'CONSTANT    '+str(self.core.dicval['MfUsgTrans'][line][0])+'\n'
            return s
        
    def writeArray(self,opt,line,ktyp):
        """writes arrays, need specific treatment for btn concentrations if pht3d
        and also for react modules of MfUsgTrans"""
        if (opt=='Pht3d') and (line == 'bct.20'):
            s = ''
            self.Conc, self.Names = self.getConcInit('main','ph.3',iper=0);
            nspec = len(self.Names)
            for i in range(nspec):
                s += self.formatBlockMt(self.Conc[i],self.Names[i]) +'\n'
        else: # normal print of array
            arr = self.core.getValueLong('MfUsgTrans',line,0);
            s = self.formatBlockMt(arr,ktyp)
        return s 

    def getConcInit(self,typ,line,iper=0):
        """returns the concentrations arrays to be printed for initial conditions
        initial chemistry has been removed
        radial also (made with the grid in usg?)
        """
        listC=[];names=[] # this will be the list of conc arrays and names of species
        dictE = self.core.addin.pht3d.getDictSpecies()
        Chem = self.core.addin.pht3d.Base['Chemistry']
        pht = self.core.getValueLong('Pht3d',line,0)
        dim = self.core.addin.getDim()
        grd = self.core.addin.getFullGrid()
        dx,ny = array(grd['dx']),int(grd['ny'])
        if typ=='rech': pht=pht[0] # only the 1st layer for recharge
        dInd={'Solutions':pht/1000.,
              'Phases':mod(pht,1000)/100,'Gases':mod(pht,1000)/100,
              'Exchange':mod(pht,100)/10,
              'Surface':mod(pht,10)}
        shortn=['k','i','kim','g','p','e','s','kp']
        longn=['Solutions','Solutions','Solutions','Gases','Phases','Exchange','Surface',
               'Phases']
        # classical case
        listC=[pht*0,pht*0,pht*0];names=['a1','a2','a3'] # unknown species!!
        for i,kw in enumerate(shortn):
            m1 = dInd[longn[i]].astype('int');
            chm = Chem[longn[i]]
            data,rows = chm['data'],chm['rows']
            for ie,e in enumerate(dictE[kw]):
                ncol = len(data[0])
                names.append(e)
                inde = rows.index(e) # finds the index of the species e
                # set background value
                rcol = list(range(2,ncol)) # the range of columns to be read
                if longn[i] in ['Phases','Gases']: 
                    m0 = pht*0.+float(data[inde][2])
                    rcol = list(range(3,ncol)) # OA 13/12/19 removed SI for asemble
                else : 
                    m0 = pht*0.+float(data[inde][1])
                if longn[i]=='Surface': rcol=list(range(2,ncol-3))
                for c in rcol:
                    if longn[i] in ['Phases','Gases']: 
                        m0[m1==(c-2)] = float(data[inde][c]) # OA 13/12/19
                    else : 
                        m0[m1==(c-1)] = float(data[inde][c])
                listC.append(m0)
        return listC,names
        
    #************************ write wells, rech, evt... ******************
    def writeWelValues(self,opt=None):
        '''a table of value that will be used in modflow'''
        if opt == 'Pht3d': 
            return self.getTransTable(opt,'ph.4','wel.1')
        else : 
            return self.getTransTable(opt,'cwell.1','wel.1')

    def writeChdValues(self,opt=None):
        if opt == 'Pht3d': 
            return self.getTransTable(opt,'ph.4','bas.5')
        else :
            return self.getTransTable(opt,'bct.20','bas.5')
               
    def getTransTable(self,opt,tline,mfline):
        '''returns a transient table containing the species, it is ordered
        in the same way as modflow (thanks to connectZones)'''
        #if tline not in self.core.ttable.keys(): return None
        ttable = self.core.ttable[tline];nt,nzo = shape(ttable)
        lsolu = unique(ttable) 
        if opt == 'Pht3d' :
            spec0,spec = self.phval2conc(len(lsolu))
        else :
            spec0,spec = ' 0.000e+00',lsolu
        lout,nzmf = self.connectZones(tline,mfline) # the index correct for modflow
        ttabl1 = np.full((nt,nzmf),spec0)
        if len(lout)>0:
            for i,solu in enumerate(lsolu):
                it,izo = where(ttable==solu)
                izo2 = array(lout)[izo]
                ttabl1[it,izo2] = spec[i] # lout[izo] allow to place the ph zone in correc tplace for modflow
        return ttabl1
        
    def phval2conc(self,nsol):
        '''returns the composition of solutions for 'k','i','kim'''
        schem = ['']*nsol;
        listE = self.core.addin.pht3d.getDictSpecies()
        chm = self.core.addin.pht3d.Base['Chemistry']['Solutions']
        data,rows = chm['data'],chm['rows']
        ncol=len(data[0]);#rcol=list(range(2,ncol)) # the range of columns to be read
        for il in range(nsol): # fills the matrix with the solutions
            schem[il] = ' %9.3e %9.3e %9.3e' %(1e-15,1e-15,1e-15)
        for kw in['k','i','kim']:
            if len(listE[kw])==0 : continue
            for e in listE[kw]:
                if e in ['pH','pe']: continue
                inde = rows.index(e) # finds the index of the species e
                for il in range(nsol): # fills the matrix with the back and solutions
                    schem[il] += ' %9.3e' %float(data[inde][il+1])
        schem0 = ' 0.000e+00'*len(schem[0].split())
        return schem0,schem
        
    def connectZones(self,tline,mfline):
        '''connect the modflow zones to the pht3d zones, returns a list of
        indices that correspond to the order of the ones existing in modflow
        the first number is the ph.4 zone number corresponding to the 1st modlfow zone
        '''
        if tline[:2]=='ph': 
            dicz = self.core.diczone['Pht3d'].dic['ph.4']
            litrans = zptsIndices(self.core,dicz)
        else : 
            dicz = self.core.diczone['MfUsgTrans'].dic[tline] # OA 3/5/20
            litrans = zptsIndices(self.core,dicz)
        dicz = self.core.diczone['Modflow'].dic[mfline]
        limod = zptsIndices(self.core,dicz)
        lout = []
        for i0,ipts in enumerate(limod):
            for ipt2 in litrans: # the zone in pht
                a = 0
                for coo in ipt2 :
                    if coo in ipts: a += 1 # nb of points recognized in ipts
                if a == len(ipt2):
                    lout.append(i0)
        return lout,len(limod)
                    
    def writeRchString(self,opt=None):
        '''this list of strings will be used by modflow, one string for each
        stress period'''
        ls = []
        for iper in range(self.nper): 
            if (opt=='Pht3d'):
                s = ''
                self.Conc, self.Names = self.getConcRch('ph.3',iper=0);
                nspec = len(self.Names)
                for i in range(nspec):
                    s += self.formatBlockMusg(self.Conc[i],self.Names[i])
                ls.append(s)
            else :
                m = block(self.core,'MfUsgTrans','crch.1',False,None,iper);
                ls.append(self.writeVecModflow(m[0],'arrfloat'))
        return ls

    def getConcRch(self,typ,line,iper=0):
        """returns the concentrations arrays to be printed for recharge
        order : 'k','i','kim','g',
        """
        listC=[];names=[] # this will be the list of conc arrays and names of species
        listE=self.core.addin.pht3d.getDictSpecies()
        chm = self.core.addin.pht3d.Base['Chemistry']['Solutions']
        data,rows = chm['data'],chm['rows']
        ncol=len(data[0]);rcol=list(range(2,ncol)) # the range of columns to be read
        pht = block(self.core,'Pht3d',line,False,None,iper).astype('int') # solution number
        for kw in['k','i','kim']:
            for e in listE[kw]:
                names.append(e)
                inde = rows.index(e) # finds the index of the species e
                m0 = pht*0.+float(data[inde][1]) # fills with background
                for c in rcol: # fills the matrix with the solutions
                    m0[pht==(c-1)]=float(data[inde][c])
                listC.append(m0)
        # gases
        gases = self.core.addin.pht3d.Base['Chemistry']['Gases']
        data,rows = gases['data'],gases['rows']
        if len(gases['rows'])>0:
            ncol = len(data[0]);rcol=list(range(2,ncol,2))
            for e in listE['g']:
                names.append(e)
                inde = rows.index(e) # finds the index of the species e
                m0 = pht*0.+float(data[inde][2]) # fills with background
                for c in rcol: # fills the matrix with the solutions
                    m0[pht==(c/2)]=float(data[inde][c])
                listC.append(m0)
        # # other species, set to 0 in recharge    
        # for kw in ['p','e','s','kp']:
        #     for e in listE[kw]:
        #         names.append(e)
        #         listC.append(m0*0.)
        return listC,names
            
    #************************ file for boundary conditions ******************
    def writePcbFile(self,core):
        # finds if the zones are transient and their values
        nzones,line = 0,'pcb.2'
        dicz = core.diczone['MfUsgTrans'].dic[line]
        if line in self.ttable:
            clist = self.ttable[line];
            nzones = len(dicz['name'])
        # find the nodes of the zones
        nbe = core.addin.mesh.getNumber('elements');
        pindx,nbelts = zeros(nbe),0 # this will be the index, 0 for the background, then poly number
        for iz in range(nzones):
            media = dicz['media'][iz]
            idx = zmesh(core,dicz,media,iz)
            try : len(idx) 
            except TypeError : continue # the zone media is not the right one
            pindx[idx] = iz+1
            nbelts += sum(idx==1)
        #write the values for every period
        s1 = str(nbelts*self.nper)+' 0\n' # MXPCB IPCBCB
        for ip in range(self.nper):
            s1 += str(nbelts).rjust(10)+'  0\n' # ITMP MP
            for iz in range(nzones):
                nodelist = where(pindx==iz+1)[0]
                for n in nodelist : 
                    s1 += str(n+1).rjust(10)+ '        1 ' # Node(Cell) iSpecies_No
                    s1 += str(clist[ip,iz]).rjust(10)+'\n' # Conc
        f1=open(self.fullPath+'.pcb','w')            
        f1.write(s1);f1.close()
        
    #********************************* write Pht3d file ********************
    def writePhFile(self,core,listE,parmk):
        # order of species 'k','i','kim','g','p','e','s','kp'
        # writes the first two lines
        s = ' 2  25  1   0   0\n0 \n' # !!!!! to be modified
        # dicval = core.dicval['Pht3d']
        # for v in dicval['ph.1'] : 
        #     s += str(v).rjust(10)
        # s += '\n'+str(dicval['ph.2'][0]).rjust(10)+'\n'
        # determine nb of kinetic species and make a list
        nInorg=len(listE['i']); #+len(self.lists);
        nKmob = len(listE['k'])
        nKimob = len(listE['kim']);
        nMinx= len(listE['p'])
        nExch=len(listE['e'])
        nGas=len(listE['g'])
        nSurf=len(listE['s'])
        s += '%3i\n' %(nInorg)# PH3 nb inorg compounds
        # PH4 nb minerals and gases
        nMin2 = nMinx+nGas
        s += '%3i\n' %(nMin2) #, nGas)
        s += '%3i %3i\n' %(nExch, 0)# PH5 nb ech ions, 0 je sais pas pourquoi
        # PH6 surface complexation 
        s += '%3i\n' %nSurf
        # PH7 Record: NR_MOB_KIN NR_MIN_KIN NR_SURF_KIN NR_IMOB_KIN
        # nb of kinetic species mobiles, minerales, surfaces et substeps pour plus tard
        nKsurf,nKmin = 0,len(listE['kp'])
        s += '%3i %3i %3i %3i\n' %(nKmob, nKmin, nKsurf, nKimob)
        # PH8 : NR_OUTP_SPEC (complexes) PR_ALKALINITY_FLAG (futur)
        #s += '%9i %9i\n' %(0,0)
        ## nb units are useless because pht3d considers that it is always days
        Chem = core.addin.pht3d.Base['Chemistry']
        if 'Rates' in Chem:
            rates=Chem['Rates'];
            for nom in listE['k']:
                iek = rates['rows'].index(nom);
                s += nom+'%5i \n '%(parmk[nom])  #param k
                for ip in range(parmk[nom]):
                    s += str(rates['data'][iek][ip+2])+' \n'
                if type(rates['data'][iek][-1])==type(0.):
                    self.core.gui.onMessage(nom+' missing formula')
                s += '-formula '+rates['data'][iek][-1] +'\n' # formula
        for n in listE['i']:
            add='';
            #if optsu.strip() == n: add=' charge'
            s += n.replace('(+','(')+add+'\n' # que reac isntant et AE, pH pe
        # kinetic immobile
        if 'Rates' in Chem:
            rates=Chem['Rates'];
            for nom in listE['kim']:
                iek = rates['rows'].index(nom);
                s += nom+'%5i \n '%(parmk[nom])  #param k
                for ip in range(parmk[nom]):
                    s += str(rates['data'][iek][ip+2])+' \n'
                if type(rates['data'][iek][-1])==type(0.):
                    self.core.gui.onMessage(nom+' missing formula')
                s += '-formula '+rates['data'][iek][-1] +'\n' # formula
        for p in listE['g']:
            ip=Chem['Gases']['rows'].index(p)
            s += p+'  '+str(Chem['Gases']['data'][ip][1])+' \n' # phase name + SI backgrd
        for p in listE['p']:
            ip=Chem['Phases']['rows'].index(p)
            s += p+'  '+str(Chem['Phases']['data'][ip][1])+' \n' # phase name + SI backgrd
        # exchanger
        for n in listE['e']: s += n+' -1 \n' 
        if 'Surface' in Chem: # surface
            su=Chem['Surface']
            for esp in listE['s']:
                icol=su['cols'].index('Specif_area')
                st=esp;ies=su['rows'].index(esp)
                st+=' '+str(su['data'][ies][icol])
                st+=' '+str(su['data'][ies][icol+1])
                if su['data'][ies][icol+2]!='': st+=' '+su['data'][ies][icol+2]  # name
                if su['data'][ies][icol+3]!='': st+=' '+su['data'][ies][icol+3]  # name
                s += st+' \n' 
            if len(listE['s'])>0 : 
                su_opt = self.core.getValueFromName('Pht3d','SU_OPT')
                if su_opt in ['no_edl','diffuse_layer','Donnan','cd_music']:
                    s+= '-'+su_opt+'\n'
        if 'Kinetic_Minerals' in Chem:
            kp=Chem['Kinetic_Minerals']
            for nom in listE['kp']:
                iek = kp['rows'].index(nom);
                s += nom+'%5i \n '%(parmk[nom])  #param k
                for ip in range(parmk[nom]):
                    s += '%9.5e \n' %float(kp['data'][iek][ip+1])
        f1 = open(self.fDir+os.sep+'Pht3d_ph.dat','w')
        f1.write(s)
        f1.close()
        

    #  '''''''''''''''' fonction writeBlockMusg '''''''''''''''''''''''''''
    #------------------------- fonction  writevect, writemat -------------------
    def formatVecMt(self, v,ktyp):
        #print shape(v),amin(v),amax(v)
        l=len(v);ln=3;s=''
        if ktyp[3:]=='int': typ='I' #OA 1/8/17
        else : typ='G'
        if amin(v)==amax(v):
            if typ=='I': s += 'CONSTANT     %9i  ' %amin(v)
            else : s += 'CONSTANT     %9.5e  ' %amin(v)
            return s
        if typ=='I': fmt='1    ('+str(l)+'I'+str(ln)
        else : fmt='0    ('+str(l)+'G12.4'           
        s += 'INTERNAL     '+fmt+')     3 \n'
        
        if typ=='I': fmt='%'+str(ln)+'i'
        else : fmt='%+11.4e '            

        for i in range(l):
            s += fmt %v[i]
        return s

    def formatMatMt(self, m, ktyp):
        #print 'mfw',shape(m),m
        if len(shape(m))==1: return self.formatVecMt(m,ktyp)
        [l,c] = shape(m);ln=3
        if ktyp[3:]=='int': typ='I' #OA 1/8/17
        else : typ='G'
        s = ''
        if amin(amin(m))==amax(amax(m)):
            if typ=='I': s += 'CONSTANT     %9i  ' %(amin(amin(m)))
            else : s += 'CONSTANT     %9.5e  ' %(amin(amin(m)))
            return s
        if typ=='I':
            fmt='1    ('+str(c)+'I'+str(ln)
        else :
            fmt='0    ('+str(c)+'G12.4' #+str(ln)            
        s += 'INTERNAL     '+fmt+')     3  \n'      
        if typ=='I':
            fmt='%'+str(ln)+'i'
        else :
            fmt='%+11.4e ' #'+str(ln)+'e '            
        for i in range(l-1,-1,-1): # to write the rows from top to bottom
            for j in range(c):
                s+=fmt %(m[i][j])
            s+='\n'
        return s[:-1]

    def formatBlockMt(self,m,ktyp):
        #print shape(m),m
        s = ''
        if len(shape(m))==3:
            nlay,a,b=shape(m);
            for l in range(nlay):
                s += self.formatMatMt(m[l],ktyp)
                if l<nlay-1: s += '\n'
        elif self.core.addin.mesh != None: # unstructured case write nlay vectors
            nlay,a=shape(m);
            for l in range(nlay):
                s += self.formatVecMt(m[l],ktyp)
                if l<nlay-1: s += '\n'       
        else : 
            s = self.formatMatMt(m,ktyp)
        return s
              
class mtUsgReader:

    """ this is the reder of UCN files """
    def __init__(self, fDir, fName):
        self.fDir = fDir
        self.fName = fName

    def readUCN(self,core,opt,iper,iesp,specname=''): 
        """ read .conc file, here opt, iesp, specname are not used
        in free flow Thksat from flo file must be added (not done)"""  
        if opt=='Pht3d': 
            return self.readACN(core,opt,iper,iesp)
        if core.mfUnstruct and core.getValueFromName('Modflow','MshType')>0:
            nlay,ncell = getNlayers(core),core.addin.mfU.getNumber('elements') # only 1 layer up to now
            cnc=zeros((nlay,ncell));#print('mfw 491', shape(cnc))
        else :
            nlay,ncol,nrow = self.getGeom(core)
            ncell = ncol*nrow
            cnc=zeros((nlay,nrow,ncol))
        try : f1 = open(self.fDir+os.sep+self.fName+'.conc','rb')
        except IOError: return None
        blok=44+ncell*4; # v210 60
        for il in range(nlay):
            f1.seek(iper*nlay*blok+blok*il+44)
            data = arr2('f')
            data.fromfile(f1,ncell)
            if core.mfUnstruct and core.getValueFromName('Modflow','MshType')>0: 
                cnc[il] = data
            else : 
                cnc[il] = reshape(data,(nrow,ncol)) #
                if core.mfUnstruct == False : cnc[il] = cnc[il][::-1] #OA 23/4/2
        f1.close()  
        return cnc
    
    def readACN(self,core,opt,iper,iesp,specname=''): 
        if opt != 'Pht3d': return
        if iesp<9: f1 = self.fDir+os.sep+'PHT3D00'+str(iesp+1)+'.ACN'
        else : f1 = self.fDir+os.sep+'PHT3D0'+str(iesp+1)+'.ACN'
        m = loadtxt(f1);
        if core.mfUnstruct and core.getValueFromName('Modflow','MshType')>0:
            nlay,ncell = getNlayers(core),core.addin.mfU.getNumber('elements') # only 1 layer up to now
            ncellt = nlay*ncell
            cnc=reshape(m[iper*ncellt:(iper+1)*ncellt],(nlay,ncell))
        else :
            nlay,ncol,nrow = self.getGeom(core)
            ncell = nlay*nrow*ncol
            cnc=reshape(m[iper*ncell:(iper+1)*ncell],(nlay,nrow,ncol))
            cnc = cnc[:,::-1,:]
        return cnc
    
    def getGeom(self,core):
        grd = core.addin.getFullGrid()
        ncol, nrow = grd['nx'], grd['ny']
        nlay=getNlayers(core);#print iper, nlay,ncol,nrow
        if core.addin.getDim() in ['Xsection','Radial']:
            nlay=nrow;nrow=1
        return nlay,ncol,nrow
