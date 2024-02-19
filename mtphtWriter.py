# -*- coding: cp1252 -*-
from array import array as arr2
import os,time
#from phtDbase import *
from mtPhtKeywords import Mt
from geometry import *
from modflowWriter import * # OA 6/5/19

class mtphtWriter:

    def __init__(self, core,fDir, fName):
        self.core = core
        self.fDir,self.fName,self.Mkey = fDir,fName,Mt()
        self.fullPath = fDir+os.sep+fName;#print self.fullPath
        self.link = {'uzt.3':('Modflow','uzf.7'),'uzt.4':('Modflow','uzf.7')}
        self.mfloW = modflowWriter(core,fDir,fName) # OA 6/5/19

    def writeMtphtFiles(self,listEsp,opt,parmk=None):
        self.ttable = self.core.makeTtable();#print 'mtpht ttable',self.ttable
        self.radfact = 1.
        self.dim = self.core.addin.getDim()
        self.core.updateDicts()
        self.rct = 0
        if 'RCT' in  self.core.getUsedModulesList('Mt3dms'): self.rct = 1
        self.writeNamFile(opt)
        tlist = array(self.ttable['tlist'])
        self.per = tlist[1:]-tlist[:-1]
        self.nper = len(self.per)#;print('writempht l.26',self.nper)
        mcomp,ncomp,gcomp = listEsp['mcomp'],listEsp['ncomp'],listEsp['gcomp']
        self.nesp = ncomp
        nkim = len(listEsp['kim'])
        # things need to be adjusted
        self.core.setValueFromName('Mt3dms','NCOMP',ncomp)
        self.core.setValueFromName('Mt3dms','MCOMP',mcomp)
        self.core.setValueFromName('Mt3dms','GCOMPN',gcomp)
        self.core.setValueFromName('Mt3dms','KCOMPN',nkim)
        self.core.setValueFromName('Mt3dms','NPER',self.nper)
        self.writeFiles(opt)
        if opt=='Pht3d':
            self.writePhFile(self.core,listEsp,parmk)
            self.writePhreeqc(self.core,listEsp);
        self.writeSsmFile(self.core,opt)
        return True #EV 04/01/22
        
    def writeNamFile(self,opt):
        f1=open(self.fDir+os.sep+opt +'.nam','w');#print 'mtpw l41',opt
        f1.write(' List  7  '+opt+'.out\n')
        if 'VDF' in  self.core.getUsedModulesList('Mt3dms') and opt=='Mt3dms':
            f2=open(self.fullPath+'.nam','r');f2.readline()
            for ll in f2: f1.write(ll)
            f2.close()
            f1.write(' vdf  47 '+opt+'.vdf\n')
            if 'VSC' in  self.core.getUsedModulesList('Mt3dms'):
                f1.write(' vsc  48 '+opt+'.vsc\n')
        else : 
            f1.write(' FTL  66 '+self.fName+'.flo\n')
        f1.write(' btn  41 '+opt+'.btn\n adv  42    '+opt+'.adv\n')
        f1.write(' dsp  43 '+opt+'.dsp\n ssm  44    '+opt+'.ssm\n')
        f1.write(' gcg  45 '+opt+'.gcg\n')
        if opt=='Pht3d':
            f1.write(' PHC  64    Pht3d_ph.dat\n')
        if self.rct>0:
            f1.write(' RCT  46 '+opt+'.rct\n')
        self.uzt=False
        if 'UZT' in self.core.getUsedModulesList('Mt3dms'):
            self.uzt = True
            f1.write(' UZT  9 '+opt+'.uzt\n')
        f1.close()

    #*********************** generic file writer ****************
    def writeFiles(self,opt):
        """to write all modflow file.
        reads the keyword file and prints all keywords by types : param (0D)
        vector (1D) array (2D). types are found by (dim1,dim2).."""
        lexceptions=['adv.1','dsp.1','dsp.3','dsp.4','dsp.5','btn.7','btn.8',
            'btn.12','btn.21','btn.22','rct.1']
        lexceptions.extend(['uzt.'+str(a) for a in range(2,11)])
        for grp in self.core.getUsedModulesList('Mt3dms'):
            if grp == 'RCT' and self.rct==0: continue
            if grp == 'SSMs': continue
            f1=open(self.fDir+os.sep+opt +'.'+ grp.lower(),'w')
            llist=self.Mkey.groups[grp];#print n1,name
            for ll in llist:
                cond=self.Mkey.lines[ll]['cond'];#print 'mtw 72',ll
                if self.testCondition(cond)==False : continue
                kwlist=self.Mkey.lines[ll]['kw']
                ktyp=self.Mkey.lines[ll]['type']
                lval=self.core.dicval['Mt3dms'][ll];#print 'mtw 77',self.core.dicval['Mt3dms'],lval,kwlist,ktyp
                if (ll in lexceptions) or (ktyp[0][:3]=='lay'):
                    self.writeExceptions(ll,kwlist,ktyp[0],f1,opt)
                    continue
                if ktyp[0] in ['vecint','vecfloat','arrint','arrfloat']:
                    self.writeArray(f1,opt,ll,0)
                elif ktyp[0]=='title': # case of a title line
                    f1.write('#'+str(ktyp[0]).rjust(10)+'\n')
                else : # classical keywords
                    for ik in range(len(kwlist)):
                        if ik<len(lval): f1.write(str(lval[ik]).rjust(10))
                        else : f1.write('0'.rjust(10))
                    f1.write('\n')
            f1.close()
            #print grp+' written'

    def testCondition(self,cond):
        """ test if the condition is satisfied"""
        return self.core.testCondition('Mt3dms',cond)
        
    def writeArray(self,f1,opt,line,ik):
        """writes arrays, need specific treatment for btn concentrations if pht3d
        and also for react modules of mt3dms"""
        grd = self.core.addin.getFullGrid()
        dx,ny = array(grd['dx']),int(grd['ny'])
        s = ''
        #print 'mt113',opt,line
        if (opt=='Pht3d') and (line in ['btn.13','rct.2c']):
            self.Conc, self.Names = self.getConcInit('main',line,iper=0);#print self.Conc
            nspec = len(self.Names)
            for i in range(nspec):
                s += self.formatBlockMt3d(self.Conc[i],self.Names[i])
        elif opt=='Mt3dms' and line=='btn.13': # to remove negative values
            initChem = self.core.dicaddin['InitialChemistry'];#print initChem
            tstep = 0
            #if 'tstep' in initChem['formula']:
                #tstep=int(initChem['formula'].split('tstep :')[1].split('\n')[0])
            # for RESTART seek the initial conc in UCNs
            if 'importUCN' in initChem['formula']:
                a = initChem['tstep']
                tstep = int(a)
                arr = self.core.transReader.readUCN(self.core,'MT3D',tstep,0)                    
            else : 
                arr = self.correctBtn()  # OA 15/4/20
            arr[arr<0.] = mean(arr)
            s = self.formatBlockMt3d(arr,line)
        elif line in ['rct.3','rct.4','rct.5','rct.6']:
            arr,p = self.getRctParms(opt,line)
            if self.dim =='Radial' and line=='rct.4': # for the transfer coefficient in dual poro medium
                for l in range(ny): 
                    arr[l] = arr[l]*(cumsum(dx)-dx/2.)*6.28
            for isp in range(len(p)):
                s += self.formatBlockMt3d(arr*p[isp],str(isp+1)+line)
        else: # normal print of array
            arr = self.core.getValueLong('Mt3dms',line,ik);#print 'mtw',line,ik,value
            s = self.formatBlockMt3d(arr,line)
        f1.write(s)
        
    def correctBtn(self): # OA added 15/4/20
        '''this take conc to be written in btn only if there is a -1 in btn.12
        (or in modflow bas.3 not sure)'''
        mtbc = self.core.getValueLong('Mt3dms','btn.12',0)
        arr = self.core.getValueLong('Mt3dms','btn.13',0)
        #mfbc = self.core.getValueLong('Modflow','bas.3',0)
        return arr*(mtbc!=-1) #+(mfbc==-1)) # OA 19/12/21 OA changed to !=
        
    def writeExceptions(self,line,kwlist,ktyp,f1,opt):
        """to write some things mt3d wants in a specific format"""
        nx,ny,a,b = getXYvects(self.core)
        nlay=getNlayers(self.core)
        if line == 'btn.7':
            delr = self.core.dicval['Modflow']['dis.4']
            f1.write(self.formatVecMt3d(delr,'delr'))

        if line == 'btn.8': # row height needs to be inverted
            delc = self.core.dicval['Modflow']['dis.5']
            if self.dim in ['2D','3D'] : f1.write(self.formatVecMt3d(delc[-1::-1],'delc'))
            elif self.dim=='Radial': f1.write('       0      1   \n')
            elif self.dim=='Xsection': 
                front = self.core.getValueFromName('Modflow','TOP')
                end = self.core.getValueFromName('Modflow','BOTM')
                f1.write('      0     '+str(front-end)+'\n')

        if line =='btn.12':
            arr = self.core.getValueLong('Mt3dms',line,0)
            s = self.formatBlockMt3d(abs(arr),line) # don't write -1
            f1.write(s)
            
        if line=='btn.21': # periods characteristics 4 values per period
            lval = self.core.dicval['Mt3dms']['btn.21'] # contains period sze, end time, 3 things to print
            lval2 = self.core.dicval['Mt3dms']['btn.22']
            s=''
            for ip in range(len(self.per)):
                s+=str(self.per[ip])[:10].rjust(10)+str(lval[1]).rjust(10)+\
                    str(lval[2]).rjust(10)  # OA 3/11/18
                if len(lval)>3: s+=str(lval[3]).rjust(10)+'\n'
                else : s+='\n'
                for i in range(4): s+=str(lval2[i]).rjust(10)
                s+='\n'
            f1.write(s)
            
        if line =='adv.1':
            lval = self.core.dicval['Mt3dms'][line]
            f1.write(str(lval[0]-1).rjust(10))
            for v in lval[1:]: f1.write(str(v).rjust(10))
            f1.write('\n')
            
        if line == 'dsp.1':
            self.lstDff = self.core.getValueFromName('Mt3dms','IMDIFF')
            if self.lstDff !='#':
                f1.write('$ multidiffusion \n')
            
        if line in ['dsp.3','dsp.4','dsp.5']:
            value = self.core.getValueLong('Mt3dms',line,0);#print 'mtw',line,ik,value
            s = self.formatBlockMt3d(value[0],line)
            if line=='dsp.5' and self.lstDff !='#':
                dff = self.getXiDiffusion(opt,self.lstDff)
                for i in range(len(dff)):
                    s = self.formatBlockMt3d(dff[i],line)
                    for l in range(nlay): f1.write(s)
            else : 
                f1.write(s)
                    
        if line == 'rct.1':
            # for reaction (ik=1) : 0th order(2) needs to be written as 100!
            lval = self.core.dicval['Mt3dms']['rct.1'];#print 'mptphw l169',lval
            s = ''
            for ik in range(len(kwlist)):
                val = lval[ik]*1
                if ik==1 and val==2: val=100 
                s += str(val).rjust(10)
            f1.write(s+'\n')
            
        if line == 'uzt.2': #only the top layer shall be written
            arr = self.core.getValueLong('Mt3dms','uzt.2',0);
            f1.write(self.formatBlockMt3d(arr[0],line))
        if line == 'uzt.3': # water content data have to come from modflow uzf
            arr = self.core.getValueLong('Modflow','uzf.7',0);
            f1.write(self.formatBlockMt3d(arr,line))
        if line == 'uzt.4': # sat thickness data have to come from modflow bas.5
            arr = self.core.getValueLong('Modflow','bas.5',0);
            f1.write(self.formatBlockMt3d(arr,line))
        if line in ['uzt.6','uzt.8','uzt.10']:
            #arr = self.core.getValueLong('Mt3dms',line,0)
            val = self.core.dicval['Mt3dms'][line][0] # OA 3/10/19
            if opt == 'Mt3dms':
                #print('writempht l.220',self.nper)
                s0 = '1 \n         0 '+str(val)+'\n'+'\n'.join(['-1']*self.nper)+'\n' # OA 3/10/19
            else : # pht3d get concentrations
                phline={'uzt.6':'ph.5','uzt.8':'ph.7','uzt.10':'ph.8'} # correspondance btw pt et pht
                s0 = ''
                for ip in range(self.nper) :
                    s0+='1 \n'
                    rch ,names = self.getConcRch('main',phline[line],ip)
                    for i in range(len(rch)):
                        s0 += self.formatBlockMt3d(rch[i][0],names[i])
            f1.write(s0)
                    
        if ktyp[:3] == 'lay':
            ilay=getNlayersPerMedia(self.core)
            lval = self.core.dicval['Mt3dms'][line] # in 3D, a list of values EV 26/09/19
            if line=='btn.6' : 
                lval = array(self.core.dicval['Modflow']['lpf.2']).astype('int') # EV 26/09/19 for type of layer confined or convertible
            s=''
            if self.dim == '2D' : s=str(lval[0])
            elif self.dim == '3D': # EV 26/09/19
                if len(lval)==len(ilay):
                    lval1 = [[lval[x]]*ilay[x] for x in range(len(ilay))]
                    lval= [item for sublist in lval1 for item in sublist]
                    s=' '+str(lval[0]).rjust(2) # OA 11/5/2
                    for i in range(1,len(lval)): # OA 11/5/2
                        if mod(i,40)==0: s+='\n' # OA 11/5/20 added no more than 40 values
                        s+=' '+str(int(lval[i])).rjust(2)
                else : 
                    s=' '+str(lval[0]).rjust(2)
                    for i in range(1,nlay): # OA 11/5/20
                        if mod(i,40)==0: s+='\n'
                        s+=' '+str(lval[0]).rjust(2)
            else : # radial and xsection
                s=str(lval[0]).rjust(2)
                for i in range(1,ny):
                    if mod(i,40)==0: s+='\n'
                    s+=str(lval[0]).rjust(2) 
            f1.write(s+'\n')     

    def getConcInit(self,typ,line,iper=0):
        """returns the concentrations arrays to be printed for initial conditions"""
        listC=[];names=[] # this will be the list of conc arrays and names of species
        dictE = self.core.addin.pht3d.getDictSpecies()
        Chem = self.core.addin.pht3d.Base['Chemistry']
        if line=='btn.13': pht = self.core.getValueLong('Pht3d','ph.3',0) #OA 3/7/22
        elif line=='rct.2c': pht = self.core.getValueLong('Mt3dms',line,0)*1000 #OA 3/7/22
        dim = self.core.addin.getDim()
        grd = self.core.addin.getFullGrid()
        dx,ny = array(grd['dx']),int(grd['ny'])
        if typ=='rech': pht=pht[0] # only the 1st layer
        dInd={'Solutions':pht/1000.,
              'Phases':mod(pht,1000)/100,'Gases':mod(pht,1000)/100,
              'Exchange':mod(pht,100)/10,
              'Surface':mod(pht,10)}
        shortn=['k','i','kim','g','p','e','s','kp']
        longn=['Solutions','Solutions','Solutions','Gases','Phases','Exchange','Surface',
               'Phases']
        initChem = self.core.dicaddin['InitialChemistry'];#print initChem
        # for RESTART seek the initial conc in UCNs
        if 'importUCN' in initChem['formula']:
            if 'All' in initChem['name']:
                a = initChem['tstep']
                tstep = int(a) 
                ie = 0; #print tstep,ie
                for kw in shortn:
                    for e in dictE[kw]:
                        names.append(e)
                        value = self.core.transReader.readUCN(self.core,'PHT3D',tstep,ie)
                        listC.append(value)#[:,::-1,:]) EV 14/05/19
                        ie += 1
                return listC,names
        # for RESTART from acsii file ############# TEMPORARY DEVELOPMENT ##############
        names_imp=[];listC_imp=[] #EV 14/02/19
        if initChem['name']!='': #EV 05/02/19
            names_imp=[];listC_imp=[] #EV 14/02/19
            ntxt=initChem['name']
            ntxt=ntxt.replace('core','self.core')
            value=self.core.formExec(ntxt)
            names_imp=value #EV 14/02/19
            vtxt= initChem['formula']
            vtxt=vtxt.replace('core','self.core')
            value=self.core.formExec(vtxt)
            listC_imp=value #EV 14/02/19
            #return listC,names
        # classical case
        listC=[];names=[]
        for i,kw in enumerate(shortn):
            m1=dInd[longn[i]].astype('int');#print 'mtph l 238',i,kw,longn[i],m1
            chm=Chem[longn[i]]
            data,rows = chm['data'],chm['rows']
            for ie,e in enumerate(dictE[kw]):
                if e in names_imp: #EV 14/02/19
                    ind=names_imp.index(e)
                    names.append(e)
                    listC.append(listC_imp[ind])
                else:
                    ncol=len(data[0])
                    names.append(e)
                    inde=rows.index(e) # finds the index of the species e
                    # set background value
                    rcol=list(range(2,ncol)) # the range of columns to be read
                    if longn[i] in ['Phases','Gases']: 
                        m0=pht*0.+float(data[inde][2])
                        rcol=list(range(3,ncol)) # OA 13/12/19 removed SI for asemble
                    else : 
                        m0=pht*0.+float(data[inde][1])
                    if longn[i]=='Surface': rcol=list(range(2,ncol-3))
                    for c in rcol:
                        if longn[i] in ['Phases','Gases']: 
                            m0[m1==(c-2)]=float(data[inde][c]) # OA 13/12/19
                        else : 
                            m0[m1==(c-1)]=float(data[inde][c])
                    if dim =='Radial' and longn[i] in ['Phases','Gases','Exchange']:
                        for l in range(ny): 
                            m0[l] = m0[l]*(cumsum(dx)-dx/2.)*6.28
                    if dim =='Radial' and longn[i]=='Surface':
                        #scol=chm['cols'].index('Specif_area') # comment 3/7/22
                        #if data[inde][scol+2]=='': # comment 3/7/22
                        for l in range(ny): 
                            m0[l] = m0[l]*(cumsum(dx)-dx/2.)*6.28
                    listC.append(m0)
        return listC,names
        
    def getConcRch(self,typ,line,iper=0):
        """returns the concentrations arrays to be printed for recharge"""
        listC=[];names=[] # this will be the list of conc arrays and names of species
        listE=self.core.addin.pht3d.getDictSpecies()
        chm = self.core.addin.pht3d.Base['Chemistry']['Solutions']
        data,rows = chm['data'],chm['rows']
        ncol=len(data[0]);rcol=list(range(2,ncol)) # the range of columns to be read
        #pht = block(self.core,'Pht3d',line,[3],None,iper).astype('int') # solution number #EV 04/02/20
        pht = self.core.getValueLong('Pht3d',line,0,iper).astype('int') #EV 04/02/20
        for kw in['k','i','kim']:
            for e in listE[kw]:
                names.append(e)
                inde=rows.index(e) # finds the index of the species e
                m0=pht*0.+float(data[inde][1]) # fills with background
                for c in rcol: # fills the matrix with the solutions
                    m0[pht==(c-1)]=float(data[inde][c])
                listC.append(m0)
        # gases
        gases = self.core.addin.pht3d.Base['Chemistry']['Gases']
        data,rows = gases['data'],gases['rows']
        if len(gases['rows'])>0:
            ncol=len(data[0]);rcol=list(range(2,ncol,2))
            for e in listE['g']:
                names.append(e)
                inde=rows.index(e) # finds the index of the species e
                m0=pht*0.+float(data[inde][2]) # fills with background
                for c in rcol: # fills the matrix with the solutions
                    m0[pht==(c/2)]=float(data[inde][c])
                listC.append(m0)
        # other species, set to 0 in recharge    
        for kw in ['p','e','s','kp']:
            for e in listE[kw]:
                names.append(e)
                listC.append(m0*0.)
        return listC,names
        
    def getRctParms(self,opt,line):
        """returns the parms for rct modules for several species. two cases occur
        species are only in mt3dms or they come from pht3d"""
        flag = self.core.addin.getMtSpecies()['flag']
        kw = self.core.dickword['Mt3dms'].lines[line]['kw'][0].split('(')[0]
        react = self.core.addin.getMtReact()
        givenE = react['rows']
        ngiv = len(givenE)
        arr = self.core.getValueLong('Mt3dms',line,0)
        icol = react['cols'].index(kw)
        if flag =='Pht3d': #writes values for each pht3d species
            listE = self.core.addin.pht3d.getListSpecies();#print 'mtpht l 262',listE
            nspec = len(listE)
            p = zeros(nspec)
            for i in range(ngiv): #find given species in listE and set value
                p[listE.index(givenE[i])] = float(react['data'][i][icol])
        elif flag =='Mt3dms':
            p = [float(react['data'][i][icol]) for i in range(ngiv)]            
        #print 'rct parms',arr,p
        return arr,p
        
    def getXiDiffusion(self,opt,listIn):
        """to transform a list of xidiff for some species in a formatted list containing all
        species for pht3d"""
        if '\n' in listIn: 
            d = listIn.split('\n')
            defined = [a.split() for a in d]
        else : 
            defined = [listIn.split()]
        if opt =='Pht3d':
            listE = self.core.addin.pht3d.getListSpecies()
            nspec = len(listE)
            dff = ones(nspec)*self.core.getValueFromName('Mt3dms','DMCOEF')
            for i in range(len(defined)):
                dff[listE.index(defined[i][0])] = defined[i][1]
            #print dff
            return dff
        
    #************************ file for transient data *************************************
    def writeSsmFile(self,core,opt):
        f1=open(self.fDir+os.sep+opt+'.ssm','w')
        nx,ny,xvect,yvect = getXYvects(core);nlay=getNlayers(core)
        BC = core.getValueLong('Modflow','bas.3',0);
        BCmt = core.getValueLong('Mt3dms','btn.12',0)
        SSMspec = core.getValueLong('Mt3dms','ssms.1',0) # specific conditions given in ssm
        # for ghb,riv and drn l 423 to 436 added OA 6/5/19
        dicSSM = {};lSSM = ['wel.1','drn.1','riv.1','ghb.1']
        dicz = self.core.diczone['Modflow']
        for n in lSSM :
            dicSSM[n] = BC*0  # OA 8/7/22 nx,ny to BC*0
            for iz in range(dicz.getNbZones(n)): 
                ilay,irow,icol,zvect = self.mfloW.xyzone2Mflow(core,n,iz)
                dicSSM[n][ilay,irow,icol] = 1
            dicSSM[n] = dicSSM[n][:,::-1,:]
        nwells,nBC = sum(ravel(dicSSM['wel.1'])),sum(ravel(BC)==-1)
        mxpts = (nBC+nwells)*self.nper;#print type(BC),nBC,mxpts,self.nper
        if opt=='Pht3d' : self.createConcStrings()
        # writes first line with flags
        flg = {}
        for n in ['WEL','DRN','RCH','EVT','RIV','GHB']:
            if n in core.getUsedModulesList('Modflow'): 
                f1.write(' T')
                flg['i'+n] =True # modfi OA 3/10
            else : 
                f1.write(' F')
                flg['i'+n] = False # modfi OA 3/10
        f1.write('\n')
        
        if opt == 'Mt3dms' :  line = 'btn.13' # line for concentrations
        elif opt =='Pht3d' : line = 'ph.4'  # line for concentrations for ssm
        nzones = 0
        if line in self.ttable:
            clist = self.ttable[line];#print "mfw transient",line,clist
            zones = core.diczone[opt].dic[line]
            nzones = len(zones['name'])
        lpts,buff = [],''
        
        typ = list(range(nzones));npts=0
        for iz in range(nzones):
            lpts.append([])
            xy = zones['coords'][iz];s='';
            x,y = list(zip(*xy))
            z=x*1 # dummy value for zone2index
            icol,irow,a = zone2index(core,x,y,z);
            irow,icol = where(fillZone(nx,ny,icol,irow,a)>0)
            n0 = len(icol)
            imed = core.diczone[opt].getMediaList(line,iz)
            ilay = media2layers(core,imed)
            #print 'ssm l 318',iz,ilay,irow,icol
            dm = core.addin.getDim()
            if dm in ['3D','2D']:
                icol, irow = list(icol)*len(ilay),list(irow)*len(ilay)
                ilay = list(ilay)*n0;ilay.sort()
            if dm in ['Xsection','Radial']:
                ilay=[ny-x-1 for x in irow];irow=[0]*len(irow);ir2=irow*1 #[ny-x-1 for  x in irow]
            else : 
                ir2= [ny-x-1 for x in irow] # change to modflow coordinates
            #print 'ssm l 327',iz,ilay,ir2,icol
            # define type (icol,ilay are in ipht3d coords, ir2 in modflow coords) and add to string
            npts+=len(irow)
            for i in range(len(irow)):
                typ[iz] = SSMspec[ilay[i],irow[i],icol[i]] # set to the specified ssm value
                if typ[iz]==-1: # nothing was specified in ssm, then search for other sources
                    if BC[ilay[i],irow[i],icol[i]]==-1: typ[iz] = 1
                    if dicSSM['wel.1'][ilay[i],irow[i],icol[i]]!=0 : typ[iz] = 2
                    if dicSSM['drn.1'][ilay[i],irow[i],icol[i]]!=0 : typ[iz] = 3 #added OA 6/5/19
                    if dicSSM['riv.1'][ilay[i],irow[i],icol[i]]!=0 : typ[iz] = 4 #added OA 6/5/19
                    if dicSSM['ghb.1'][ilay[i],irow[i],icol[i]]!=0 : typ[iz] = 5 #added OA 6/5/19
                #print iz,i,typ,ilay[i],irow[i],icol[i],BC[ilay[i],irow[i],icol[i]],wells[ilay[i],irow[i],icol[i]]
                if opt=='Pht3d' and typ[iz]==-1: continue
                if opt=='Mt3dms' and  typ[iz]==-1 and BCmt[ilay[i],irow[i],icol[i]]!=-1: continue
                s= str(ilay[i]+1).rjust(9)+' '+str(ir2[i]+1).rjust(9)+' '+str(icol[i]+1).rjust(9)
                lpts[iz].append(s)
            #print mxpts,self.nper
            mxpts += npts*self.nper
        f1.write(' %6i     \n' %mxpts)

        for ip in range(self.nper):
            #print 'mfw trans',ip
            npts, s0,s1 = 0, '',''
            if flg['iRCH']: # modfi OA 3/10
                if opt =='Mt3dms' :
                    #rch = block(self.core,opt,'btn.23',intp=False,opt=None,iper=ip)#EV 04/02/20
                    rch = self.core.getValueLong('Mt3dms','btn.23',0,iper=ip) #EV 04/02/20
                    s0 = '    1\n' + self.formatBlockMt3d(rch[0],'rech')
                if opt =='Pht3d':
                    rch ,names = self.getConcRch('main','ph.5',ip)
                    nsp, s0 = len(rch),'    1\n'
                    for i in range(nsp):
                        s0 += self.formatBlockMt3d(rch[i][0],names[i])
            if flg['iEVT']: # modfi OA 3/10
                s0 += '   1\n'
                for i in range(self.nesp):
                    s0 += '       0       0        #evt\n'
            for iz in range(nzones):
                val = clist[ip,iz]; 
                if float(val)<0: 
                    val = str(-float(val));typ[iz]=15 # for mass flux
                #print 'mtpht ssm',ip,iz,val,typ[iz],len(lpts[iz])
                if opt == 'Mt3dms': val1 = ' '+str(val).rjust(9)
                elif opt == 'Pht3d' : val1 = self.concStrings[int(float(val))] 
                #print ip,iz,val
                for i in range(len(lpts[iz])):
                    s1 +=lpts[iz][i]+' '+val.rjust(9)+' '+str(typ[iz]).rjust(9)+val1+'\n'
                    npts += 1
                    mxpts += 1
            f1.write(s0+' %6i     \n'%npts +s1)

        f1.close()
        #print 'SSM written'
        
    def createConcStrings(self):
        """creates a list of strings of concentrations from pht3d solutions for the ssm"""
        self.concStrings=[] # needed to remenber conc strings not to calculate each time
        dicE=self.core.addin.pht3d.getDictSpecies()
        Chem=self.core.addin.pht3d.Base['Chemistry']
        data, rows = Chem['Solutions']['data'],Chem['Solutions']['rows']
        ncol=len(data[0]);#print 'in concstrings',data[0],rows,ncol
        gdata, grows = Chem['Gases']['data'],Chem['Gases']['rows']
        for c in range(1,ncol):
            s=''
            for kw in ['k','i','kim']: # kinetic, instantaneous, kinet immobile
                for e in dicE[kw]:
                    inde=rows.index(e)
                    s+=' %11.5e'%float(data[inde][c]) #.rjust(11)
            for e in dicE['g']: # gases
                inde = grows.index(e)
                s+=' %11.5e'%float(gdata[inde][c*2-1]) #.rjust(11)
            for kw in ['p','e','s','kp']:
                for e in dicE[kw]:
                    s+='0.'.rjust(10)
            self.concStrings.append(s)        
        
    #********************************* Ecire fichier VDF for Seawat ********************
    def WriteVdfFile(self):
        f1 = self.CreateFile('Seawatvdf.dat','w')
##        1. MTDNCONC MFNADVFD NSWTCPL IWTABLE
##        Mtdconc=0 density specified =n dens calculated from n species
##        MFNADVFD=2 centre in sapce, ><2 upstream weigthed
##        NSWTCPL max number of non linear iterations if 0 or 1 explicit coupling
##        iwtable : 0 water table correction not applied, >0 applied
        Mtdnconc,Mfnadvdf,Nswtcpl,Iwtable=1,1,1,0
        f1.write(' %9i %9i %9i %9i   Mtdnconc Mfnadvdf Nswtcpl Iwtable \n' %(Mtdnconc,Mfnadvdf,Nswtcpl,Iwtable))
##        2. DENSEMIN DENSEMAX min an dmax density if 0 not limitation
        Densemin,Densemax=0,0
        f1.write(' %9i %9i   DENSEMIN DENSEMAX \n' %(Densemin,Densemax))
##        If NSWTCPL is greater than 1, then read item 3.
##        3. DNSCRIT convergene criterion difference in density
        Dnscrit=1e-3
        if Nswtcpl>1: f1.write(' %9.4e  \n' %Dnscrit)
##        4. DENSEREF DENSESLP
        Denseref,Denseslp=1000.,.7143 # Care in kg/m3, slp for freash and sea water
        f1.write(' %9i %9i  DENSEREF DENSESLP \n' %(Denseref,Denseslp))
##        5. FIRSTDT ength of first time step
        Firstdt=1e-3
        f1.write(' %9.4e  \n' %Firstdt)
##        FOR EACH STRESS PERIOD (read items 6 and 7 only if MTDNCONC = 0)
##        6. INDENSE if<0 val of dense reused form prev tstep
##          =0 Dense=ref >=1 dense read from item 7 =2 read from 7 but as conc
##        Read item 7 only if INDENSE is greater than zero
##        7. [DENSE(NCOL,NROW)] – U2DREL
##        Item 7 is read for each layer in the grid.
        nper=len(self.tlist);Indense=-1
        for ip in range(nper):
            if Mtdnconc==0: f1.write(' %9i %9i  INDENSE \n' %Indense)
        f1.close()

    #********************************* Ecire fichier Pht3d ********************
    def writePhFile(self,core,listE,parmk):
        f1 = open(self.fDir+os.sep+'Pht3d_ph.dat','w')
        # writes the first two lines
        s=''
        dicval = core.dicval['Pht3d']
        for v in dicval['ph.1'] : 
            s += str(v).rjust(10)
        s += '\n'+str(dicval['ph.2'][0]).rjust(10)+'\n'
        # determine nb of kinetic species and make a list
        nInorg=len(listE['i']); #+len(self.lists);
        nKmob = len(listE['k'])
        nKimob = len(listE['kim']);
        nMinx= len(listE['p'])
        nExch=len(listE['e'])
        nGas=len(listE['g'])
        nSurf=len(listE['s'])
        s += '%9i\n' %(nInorg)# PH3 nb inorg compounds
        # PH4 nb minerals and gases
        if self.uzt: nMin2 = nMinx
        else : nMin2 = nMinx+nGas
        s += '%9i %9i\n' %(nMin2, nGas)
        s += '%9i %9i\n' %(nExch, 0)# PH5 nb ech ions, 0 je sais pas pourquoi
        # PH6 surface complexation 
        s += '%9i\n' %nSurf
        #PH7 Record: NR_MOB_KIN NR_MIN_KIN NR_SURF_KIN NR_IMOB_KIN
        # nb of kinetic species mobiles, minerales, surfaces et substeps pour plus tard
        nKsurf,nKmin = 0,len(listE['kp'])
        s += '%9i %9i %9i %9i\n' %(nKmob, nKmin, nKsurf, nKimob)
        # PH8 : NR_OUTP_SPEC (complexes) PR_ALKALINITY_FLAG (futur)
        s += '%9i %9i\n' %(0,0)
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
                else : s+= '   \n' # OA 24/4/20
        if 'Kinetic_Minerals' in Chem:
            kp=Chem['Kinetic_Minerals']
            for nom in listE['kp']:
                iek = kp['rows'].index(nom);
                s += nom+'%5i \n '%(parmk[nom])  #param k
                for ip in range(parmk[nom]):
                    s += '%9.5e \n' %float(kp['data'][iek][ip+1])
        f1.write(s)
        f1.close()

    #  ''''''''''''''''''''''''''''''' fonction writeMatMt3d '''''''''''''''''''''''''''
    def formatVecMt3d(self, v, name=' '):
        s = ''
        if amin(v)==amax(v):
            if amin(v)<0:
                return '         0 %+9.2e  #'%(amin(v)) +name+'\n'                
            else:
                return '         0 %9.3e  #'%(amin(v)) +name+'\n'
        l = len(v);a=str(type(v[0]))
        if a[13:16]=='int': typ='I'
        else : typ='G'
        if typ=='I': fmt='      100    1      ('+str(l)+'I'+str(ln)+')'
        else : fmt='      100    1.0      ('+str(l)+'G13.5)'            
        s += fmt +'\n'
        
        if typ=='I': fmt='%'+str(ln)+'i'
        else : fmt='%+12.5e '
        for i in range(l):
            s +=fmt %v[i]
        return s+'\n'
        
    def formatMatMt3d(self,m,name=' '):
        s = ''
        if len(shape(m))<=1: return self.formatVecMt3d(m,name)
        [l,c] = shape(m);ln=3;a=str(type(m[0,0]))
        if 'int' in a: typ='I'  # OA 3/10/18
        else : typ='G'
        if amin(amin(m))==amax(amax(m)):
            if typ=='I':s = '         0 %8i #'%(amin(amin(m))) +name[:6]+'\n'
            else: 
                if amin(amin(m))<0 : s='         0 %+9.2e #' %(amin(amin(m))) +name[:6]+'\n'
                else :s='         0 %9.3e #' %(amin(amin(m))) +name[:6]+'\n'
            return s
        if typ=='I': fmt='      100    1      ('+str(c)+'I'+str(ln)+')'
        else : fmt='      100    1.0      ('+str(c)+'G13.5) #'+name[:6]          
        s += fmt+'\n'
        
        if typ=='I': fmt='%'+str(ln)+'i'
        else : fmt='%12.5e '        
        for i in range(l-1,-1,-1):
            for j in range(c):
                s+=fmt %(m[i][j])
            s+='\n'
        return s
        
    def formatBlockMt3d(self,m,name=' '):
        s = ''
        if len(shape(m))==3:
            nlay,a,b=shape(m)
            for l in range(nlay): s+=self.formatMatMt3d(m[l],name)
        else : s = self.formatMatMt3d(m,name)
        return s

    def writePhreeqc(self,core,listE):
        """this routine writes a phreeqc file where all solutions are written in
        phreqc format to be able to test their equilibrium before running pht3d
        1. tabke background, then cycle through pht3d zones
        2. get the solution number, phase number...it does not take rates
        3 write them in phreeqc format"""
        f1 = open(self.fDir+os.sep+'solutions.pqi','w')
        f1.write('Database '+self.fDir+'\pht3d_datab.dat \n')
        chem = core.addin.pht3d.Base['Chemistry']
        solu=chem['Solutions'];
        nbsol = len(solu['cols'])
        for isol in range(1,nbsol):
            f1.write('Solution '+str(isol)+' \n units mol/L \n')
            for esp in listE['i']: # go through phase list
                ie=solu['rows'].index(esp);#print esp,phases['rows'],ip,phases['data'][ip] # index of the phase
                conc=solu['data'][ie][isol] #backgr concentration of species
                f1.write(esp+' '+str(conc)+'\n')
#            f1.write('Equilibrium_Phases '+str(isol)+'\n')
#            phases=chem['Phases'];
#            nphas=mod(val,1000)/100
#            for esp in listE['p']: # go through phase list
#                ip=phases['rows'].index(esp);#print esp,phases['rows'],ip,phases['data'][ip] # index of the phase
#                IS,conc=phases['data'][ip][nphas*2+1:nphas*2+3] #backgr SI and concentration of phase
#                f1.write(esp+' '+str(IS)+' '+str(conc)+'\n')
            f1.write('end \n')
        f1.close()
        
class mtphtReader:

    """ this is the reder of UCN files """
    def __init__(self, fDir, fName):
        self.fDir = fDir
        self.fName = fName

    def readUCN(self,core,opt,tstep,iesp,specname=''): # struct 32*4+6*4+icbund(nrow*ncol)+3*4+Cnc(nrow,ncol)+4?
        """reads ucn"""
        grd  = core.addin.getFullGrid()
        ncol,nrow = grd['nx'],grd['ny']
        nlay = getNlayers(core)
        dim = core.addin.getDim()
        if dim in ['Xsection','Radial']:
            nlay=nrow*1;nrow=1
        suff1 = 'PHT3D'
        if opt=='Mt3dms' and 'UZT' not in  core.getUsedModulesList('Mt3dms'): suff1 = 'MT3D'
        if iesp<9: suff2='00'
        else : suff2='0'
        tlist = core.ttable['tlist']
        fname=self.fDir+os.sep+suff1+suff2+str(iesp+1)+'.UCN'
        try: f1=open(fname,'rb')
        except IOError: return None
        m0=zeros((nlay,nrow,ncol))+0.
        if opt=='Mt3dms': # pb Mt3dms sometimes write more time steps...
            lt=self.getMt3dTlist(f1,ncol,nrow,nlay);
            tstep=argmin(abs(tlist[tstep+1]-lt)); # try to find the closest one
        blok=44+ncol*nrow*4;
        p0 = blok*nlay*tstep;#print 'mp3read',tstep,nlay
        #print 'ucn',fname,tstep
        for l in range(nlay):
            pos=p0+blok*l+44;f1.seek(pos);#print 'ucn',nlay,ncol,nrow,pos,p0 # v212
            data = arr2('f');data.fromfile(f1, ncol*nrow)
            m0[l]=reshape(data,(nrow,ncol))
        f1.close()
        if dim in ['Xsection','Radial']: return m0 # [::-1,:,:] #OA 11/3/19 renversement, je sais pas bien
        else : return m0[:,::-1,:]

    def getPtObs(self,core,irow,icol,ilay,iper,opt,iesp=0,specname='',ss=''): #EV 23/03/20
        """a function to values of one variable at given point or points.
        irow, icol and ilay must be lists of the same length. iper is also
        a list containing the periods for which data are needed. opt is for 
        Mt3dms or Pht3d, iesp is a list containing the indice of the species 
        (for pht3d), specname is a list containing the name of the species 
        (for min3p), ss is for solute ('') or sorbed ('S' ) species. 
        """
        grd  = core.addin.getFullGrid()
        ncol, nrow = grd['nx'], grd['ny']
        nlay=getNlayers(core);
        suff1 = opt*1
        if opt=='Mt3dms': suff1 = 'MT3D'
        if iesp<9: suff2='00'
        else : suff2='0'
        tlist=core.getTlist2()
        fname=self.fDir+os.sep+suff1+suff2+str(iesp+1)+ss+'.UCN' #EV 23/03/20
        #print('fname',fname)
        try : f1 = open(fname,'rb')
        except IOError: return None
        if core.addin.getDim() in ['Xsection','Radial']:
            nlay=nrow*1;nrow=1;ilay=irow*1;irow=[0]*len(ilay)
        #print ilay,irow,icol
        npts=len(irow)
        pobs=zeros((len(iper),npts))+0.;
        lt=self.getMt3dTlist(f1,ncol,nrow,nlay)
        blok=44+ncol*nrow*4;
        #print 'mtpht obs',iper,lt
        for ip in range(len(iper)):
            if opt=='Mt3dms':
                ip2=argmin(abs(tlist[iper[ip]]-lt))#a=arange(len(lt));;ip2=max(a[lt==iper[ip]])
            elif opt=='Pht3d' :
                ip2=iper[ip]
            p0=blok*nlay*ip2;
            for i in range(npts):
                pos = p0+blok*ilay[i]+44+irow[i]*ncol*4+icol[i]*4
                f1.seek(pos)
                try:
                    data = arr2('f');data.fromfile(f1,1);
                    pobs[ip,i]=float(data[0])
                except EOFError: return pobs
        return pobs
   
    def getMt3dTlist(self,f1,ncol,nrow,nlay):
        tread=[];i=0;blok=(44+ncol*nrow*4)*nlay
        while 3>2:
            f1.seek(blok*i+12);
            try :
                a=arr2('f');a.fromfile(f1,1);#print 'mt per',i1,a
                tread.append(a[0]);i+=1 #this is the true time given by mt3d
            except EOFError:
                #print ncol,nrow,nlay,tread;
                return array(tread)