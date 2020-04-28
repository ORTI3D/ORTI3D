#import scipy.io.array_import as IOAscii
from array import array as arr2
import os
from .modflowKeywords import Mf
from .geometry import *
from .timeperiod import *
from .modflowUsg import *

class modflowWriter:
    
    def __init__(self, core,fDir, fName):
        self.core = core
        self.fDir,self.fName,self.Fkey = fDir,fName,Mf()
        self.fullPath = fDir+os.sep+fName;#print self.fullPath

    def writeModflowFiles(self, core,usgTrans={}):
        nbfor,recharge = 0,0
        self.ttable,self.usgTrans = self.core.makeTtable(),usgTrans;
        tlist = array(self.ttable['tlist'])
        self.per = tlist[1:]-tlist[:-1]
        self.nper = len(self.per)
        self.core.Zblock = makeZblock(self.core);
        dim = self.core.addin.getDim()
        nx,ny,a,b = getXYvects(self.core)
        if dim in ['2D','3D']: self.nlay = getNlayers(self.core)
        else : self.nlay = ny
        if len(self.usgTrans.items())>0: self.trans = True
        else : self.trans = False
        self.dicval = self.core.dicval['Modflow']
        self.core.setValueFromName('Modflow','NPER',self.nper) # nper may be wrong if transietn zones
        self.radfact=1.
        if self.core.addin.getDim() in ['Radial','Xsection']: self.setRadial()
        lexceptions=['dis.1','dis.4','dis.5','dis.8','rch.2','uzf.9','uzf.10','lpf.2']#EV 06/11 added 'dis.1'
        lexceptions.extend(['disu.2','disu.3','disu.6','disu.9']) 
        lexceptions.extend(['lpf.'+str(a) for a in range(8,14)]) # when wirting by layers
        lexceptions.extend(['upw.'+str(a) for a in range(7,13)])
        lexceptions.extend(['uzf.'+str(a) for a in range(2,8)])
        lexceptions.extend(['evt.'+str(a) for a in range(2,5)])
        lexceptions.extend(['sms.1a','hfb.3'])
        self.lexceptions = lexceptions
        lgex = ['DIS','DISU','LPF','RCH','EVT','UPW','UZF','SMS','HFB6']
        lnorm = ['BAS6','SIP','PCG','SOR','DE4','NWT','GMG'] 
        self.writeNamFile()
        #self.writeFiles() #OA 13/8/19 for loop below is new
        for grp in self.core.getUsedModulesList('Modflow'):
            if grp in ['WEL','DRN','RIV','MNWT','GHB','HFB6']: continue # in transientfile or specific (MNWT&HBF6) # EV 28/08/19
            elif grp in lnorm : self.writeOneFile(grp,{})
            elif grp in lgex : exec('self.write'+grp+'()')
        self.writeLmtFile()
        self.writeOcFile()
        #print 'mfwrite',self.ttable
        if 'bas.5' in list(self.ttable['Transient'].keys()): # var head
            if self.ttable['Transient']['bas.5']:
                self.writeTransientFile(core,'bas.5','chd')
        for n in ['wel','drn','riv','ghb']:
            if self.core.diczone['Modflow'].getNbZones(n+'.1')>0 or self.core.dictype['Modflow'][n+'.1'][0]=='array':
                self.writeTransientFile(core,n+'.1',n)
        if self.core.diczone['Modflow'].getNbZones('mnwt.2a')>0:
            self.writeMNwtFile(core)
        if self.core.diczone['Modflow'].getNbZones('hbf.3')>0: # EV 28/08/19
            self.writeHFB6()  # EV 28/08/19
            
    def setRadial(self):
        g = self.core.addin.getFullGrid()
        self.core.setValueFromName('Modflow','NLAY',g['ny'])
        self.core.setValueFromName('Modflow','NROW',1) 

    """ for the y dimension ipht3d matrices are the inverse of modflow ones
    for vectors, it is taken into account any time a y vector is written
    for matrices it is in the writeblock function"""
    #************************* fichier NAM ****************
    def writeNamFile(self):
        lmod=self.core.getUsedModulesList('Modflow')
        #set solver
        for n in ['PCG','DE4','SOR','SIP','SMS']:
            if n in lmod: self.solv = n
        f1=open(self.fullPath +'.nam','w')
        f1.write('LIST   2  ' + self.fName + '.lst\n')
        for i,mod in enumerate(lmod):
            if mod in ['RCH','WEL','CHD','MNWT','HFB6']: continue # EV 28/08/19
            f1.write(mod+'  '+str(i+10)+'  ' + self.fName + '.'+mod.lower()+'\n')
        f1.write('OC       26    '+ self.fName + '.oc\n')
        f1.write('DATA(BINARY)     30        ' + self.fName + '.head\n')
        f1.write('DATA(BINARY)     31        ' + self.fName + '.budget\n')
        if 'USG' not in self.core.dicaddin['Model']['group']:
            f1.write('LMT6     28    '+ self.fName + '.lmt\n')
        r0=self.core.getValueFromName('Modflow','RECH')
        if 'RCH' in lmod:
            if r0>0 or self.core.diczone['Modflow'].getNbZones('rch.2')>0:
                f1.write('RCH     34     ' + self.fName + '.rch\n')
        e0=self.core.getValueFromName('Modflow','EVT')
        if 'EVTR' in lmod: #EV 28/02/20
            if e0>0 or self.core.diczone['Modflow'].getNbZones('evt.2')>0:
                f1.write('EVT     35     ' + self.fName + '.evt\n')
        if 'WEL' in lmod:# EV 28/08/19
            if self.core.diczone['Modflow'].getNbZones('wel.1')>0:
                f1.write('WEL     36     ' + self.fName + '.wel\n')
        if 'bas.5' in list(self.ttable['Transient'].keys()):
            if self.ttable['Transient']['bas.5']:
                f1.write('CHD     37     ' + self.fName + '.chd\n')
        if 'MNWT' in lmod: # EV 28/08/19
            if self.core.diczone['Modflow'].getNbZones('mnwt.2a')>0:
                f1.write('MNW2     39     ' + self.fName + '.mnwt\n')
        if 'HFB6' in lmod:# EV 28/08/19
            if self.core.diczone['Modflow'].getNbZones('hfb.3')>0: # EV 28/08/19
                f1.write('HFB6     38     ' + self.fName + '.hfb6\n')  # EV 28/08/19  
        if self.trans: # OA added for usg
            f1.write(self.usgTrans['nam'])
            
        f1.close()

    #*********************** specific file writers ****************
            
    def writeDIS(self):
        # dis.4
        dim = self.core.addin.getDim()
        delr = self.core.dicval['Modflow']['dis.4']
        s = self.writeVecModflow(delr,'vecfloat')
        exceptDict={'dis.4':s+'\n'}
        #dis.5
        delc = self.core.dicval['Modflow']['dis.5']
        if dim in ['2D','3D'] : s = self.writeVecModflow(delc[-1::-1],'vecfloat')
        elif dim=='Radial': s = 'CONSTANT     1  '
        elif dim=='Xsection': 
            front = self.core.getValueFromName('Modflow','TOP')
            end = self.core.getValueFromName('Modflow','BOTM')
            s = 'CONSTANT    '+str(front-end)
        exceptDict['dis.5'] =s+'\n'
        #dis.8'
        tlist = array(self.ttable['tlist'])
        perlen = tlist[1:]-tlist[:-1]
        lval=self.dicval['dis.8'] # contains period sze, end time, 3 things to print
        SsTr,s1 = 'SS',''
        if lval[3]==1 : SsTr = 'Tr'
        for ip in range(self.nper):
            if ip>0: SsTr='Tr'
            s1 +=' %9.3e %9i %9.3e ' %(perlen[ip],lval[1],lval[2]) #OA 13/9/18
            s1 += SsTr.rjust(10)+'\n'
        exceptDict['dis.8'] =s1+'\n'
        self.writeOneFile('DIS',exceptDict)
        
    def writeDISU(self):
        nodelay = int(self.core.getValueFromName('Modflow','NODELAY'))
        exceptDict={'disu.6': 'CONSTANT    '+str(nodelay)+'\n'}
        
        s1 = self.core.addin.mfU.writeDisu()
        tlist = array(self.ttable['tlist'])
        perlen = tlist[1:]-tlist[:-1]
        lval=self.dicval['disu.9'] # contains period sze, end time, 3 things to print
        SsTr = 'SS'
        if lval[3]==1 : SsTr = 'Tr'
        for ip in range(self.nper):
            if ip>0: SsTr='Tr'
            s1 +=' %9.3e %9i %9.3e ' %(perlen[ip],lval[1],lval[2]) #OA 13/9/18
            s1 += SsTr.rjust(10)+'\n'
        exceptDict['disu.9'] = s1
        self.writeOneFile('DISU',exceptDict)
        
    def writeLPF(self):
        #lpf.2
        ilay=getNlayersPerMedia(self.core) # EV 23/11/2018
        val = self.core.dicval['Modflow']['lpf.2']
        lval1 = [[val[x]]*ilay[x] for x in range(len(ilay))]
        lval= [item for sublist in lval1 for item in sublist]
        s=''
        for i in range(len(lval)):
            s+=' '+str(int(lval[i])).rjust(2) 
        exceptDict={'lpf.2':s+'\n'}
        #lpf.8' writes several lines per layer
        llist = ['lpf.'+str(a) for a in range(8,14)];#take four lines
        value = []
        for l2 in llist:
            cond = self.Fkey.lines[l2]['cond'];
            if self.testCondition(cond) == False : continue
            v0 = self.core.getValueLong('Modflow',l2,0);#print 'mfw 173', l2,v0
            value.append(v0)
        val,s = self.core.dicval['Modflow']['lpf.2'],'' # EV 23/11/2018
        ilay=getNlayersPerMedia(self.core) 
        lval1 = [[val[x]]*ilay[x] for x in range(len(ilay))]
        lval= [item for sublist in lval1 for item in sublist]
        for l in range(self.nlay):
            for i in range(len(value)):
                if i==3 and lval[l]==0 : continue # specif writing for storage
                s += self.writeMatModflow(value[i][l],'arrfloat')+'\n'
        exceptDict['lpf.8'] = s
        self.writeOneFile('LPF',exceptDict)
        
    def writeRCH(self):
        s= ''
        if 'rch.2' in self.ttable:
            trch = self.ttable['rch.2']; zrch = True
        else :
            trch = ones(self.nper); zrch = False # a constant value over the domain
        for iper in range(self.nper): 
            if (iper==0) or (prod(trch[iper]==trch[iper-1])==0): #values diff than previous
                if 'rch' in self.usgTrans.keys():  s += '    0      INCONC\n'
                else :  s += '    0\n'# 0: data are written not reused from previous
                #m = block(self.core,'Modflow','rch.2',False,None,iper); #EV 04/02/20 
                m = self.core.getValueLong('Modflow','rch.2',0,iper) #EV 04/02/20
                s += self.writeMatModflow(m[0],'arrfloat')+ '\n'
                if 'rch' in self.usgTrans.keys(): s += self.usgTrans['rch'][iper]
            else:
                if 'rch' in self.usgTrans.keys():  s += '    -1   \n' # OA 21/8/19 removed Inconc
                else :  s += '    -1\n'
        exceptDict={'rch.2':s}
        if 'rch' in self.usgTrans.keys(): optionDict = {'rch.1': ' CONC  \n       1'} # one species
        else : optionDict = {}
        self.writeOneFile('RCH',exceptDict,optionDict)
        
    def writeEVT(self):
        s = ''
        for iper in range(self.nper):
            s += '    0     0    0\n' # INSURF INEVTR INEXDP INIEVT not read
            for k in range(2,5):
                #m = block(self.core,'Modflow',line[:4]+str(k),False,None,iper);#EV 04/02/20
                m = self.core.getValueLong('Modflow','evt.'+str(k),0,iper) #EV 04/02/20 #EV 28/02/20 line[:4]
                s += self.writeMatModflow(m[0],'arrfloat')+'\n'
        exceptDict={'evt.2':s}
        self.writeOneFile('EVT',exceptDict)

    def writeUPW(self):
        #upw.7' writes several lines per layer
        strt = int(line[-1])
        llist = [line[:4]+str(a) for a in range(strt,strt+6)];#take four lines
        value = []
        for l2 in llist:
            cond = self.Fkey.lines[l2]['cond'];
            if self.testCondition(cond) == False : continue
            v0 = self.core.getValueLong('Modflow',l2,0);#print 'mfw 173', l2,v0
            value.append(v0)
        val,s = self.core.dicval['Modflow']['lpf.2'],'' # EV 23/11/2018
        ilay=getNlayersPerMedia(self.core) 
        lval1 = [[val[x]]*ilay[x] for x in range(len(ilay))]
        lval= [item for sublist in lval1 for item in sublist]
        for l in range(nlay):
            for i in range(len(value)):
                if i==3 and lval[l]==0 : continue
                s += self.writeMatModflow(value[i][l],'arrfloat')+'\n'
        exceptDict['upw.7'] = s
        self.writeOneFile('UPW',exceptDict)
        
    def writeSMS(self):        
        exceptDict = {'sms.1a':''}
        a = self.core.dicval['Modflow']['sms.1a'][0];# print 'mfw 168',s, type(s)
        if a!= 0: 
            exceptDict['sms.1a'] = self.Fkey.lines['sms.1a']['kw'][i]
        self.writeOneFile('SMS',exceptDict)
        
    def writeUZF(self):
        #uzf.9' writes an array for each period with zero before
        m0, s = 0,''
        for iper in range(self.nper): 
            #m = block(self.core,'Modflow','uzf.10',False,None,iper); #EV 04/02/20
            m = self.core.getValueLong('Modflow','uzf.10',0,iper) #EV 04/02/20
            s += '    0\n' # 0: data are written not reused from previous
            s += self.writeMatModflow(m[0],'arrfloat')+'\n'
        exceptDict = {'uzf.9':s}
        # in uzf put only one value of these parameters
        s = ''
        for line in ['uzf.2','uzf.3','uzf.4','uzf.5','uzf.6','uzf.7']: 
            m = self.core.getValueLong('Modflow',line,0)
            s += self.writeMatModflow(m[0],ktyp)+'\n'  # OA 7/3/19
            exceptDict[line] = s
        self.writeOneFile(self,'UZF',exceptDict)
        
    def writeHFB6(self):
        #hfb.3': # EV 26/11/2018
        line= 'hfb.3'
        zname = self.core.diczone['Modflow'].dic[line]['name']
        val = self.core.diczone['Modflow'].dic[line]['value'] 
        nbz = len(zname)
        nbHfb=[]; s=''
        for iz in range(nbz):
            imed = self.core.diczone['Modflow'].getMediaList(line,iz)
            ilay = media2layers(self.core,imed)
            hfb = self.writeHfb(line,iz)
            nbHfb.append(len(ilay)*len(hfb))
            for n in range(len(ilay)):
                for i in range(len(hfb)):
                    s+=str(ilay[n]+1)+' '+str(hfb[i][0])+' '+str(hfb[i][1])
                    s+=' '+str(hfb[i][2])+' '+str(hfb[i][3])+' '+str(val[iz])+'\n' 
        exceptDict = {'hfb.3':' 0  0  '+str(sum(nbHfb))+'\n'+s}
        self.writeOneFile('HFB6',exceptDict)
        
    def writeOneFile(self,grp,exceptDict,optionDict={}):
        """to write one modflow file.
        reads the keyword file and prints all keywords by types : param (0D)
        vector (1D) array (2D). types are found by (dim1,dim2).."""
        lexceptions = self.lexceptions 
        ext, s = grp, ''
        if grp=='Solver': ext=self.solv
        f1=open(self.fullPath +'.'+ ext.lower(),'w')
        llist=self.Fkey.groups[grp];
        for ll in llist:
            cond=self.Fkey.lines[ll]['cond']#;print('mfw 96',ll)
            if self.testCondition(cond)==False : continue
            kwlist=self.Fkey.lines[ll]['kw']
            ktyp=self.Fkey.lines[ll]['type'];#print 'mfw',ktyp
            lval=self.dicval[ll]
            if (ll in lexceptions):
                if ll in exceptDict.keys(): s += exceptDict[ll]
                continue
            if (ktyp[0][:3]=='lay'):
                s += self.layerLines(ll)
                continue
            for ik in range(len(kwlist)):
                value=lval[ik]
                if ktyp[ik] in ['vecint','vecfloat','arrint','arrfloat']:
                    print('mfw 106',ll,shape(value))
                    value=self.core.getValueLong('Modflow',ll,ik);
                    s += self.writeBlockModflow(value,ktyp[ik]) # OA 1/8/17
                elif ktyp[ik]=='choice': # where there is a choice print the nb othe choice not value
                    s += str(value).rjust(10)
                elif ktyp[ik]=='title': # case of a title line
                    s += '#'+str(value).rjust(10)
                else : # values with strings
                    s += str(value).rjust(10);#print 'write str',ll
            if ll in optionDict.keys(): s+= ' '+optionDict[ll]+'\n'
            else : s += '\n'
        f1.write(s);f1.close()
        
    def layerLines(self,line):
        # to print laycbd and others
        ilay=getNlayersPerMedia(self.core)
        dim = self.core.addin.getDim()
        nx,ny,a,b = getXYvects(self.core)
        lval = self.core.dicval['Modflow'][line] #in 3D lval is a list of values EV 26/09/19
        s=''
        if dim == '2D' : s=str(lval[0])
        elif dim == '3D': #EV 26/09/19
            if len(lval)==len(ilay):
                lval1 = [[lval[x]]*ilay[x] for x in range(len(ilay))]
                lval= [item for sublist in lval1 for item in sublist]
                s=''
                for i in range(len(lval)):
                    s+=' '+str(int(lval[i])).rjust(2)
            else : 
                s=' '+str(lval[0]).rjust(2)
                for i in range(1,self.nlay):
                    if mod(i,40)==0: s+='\n'
                    s+=' '+str(lval[0]).rjust(2)
        else : # radial and xsection
            s=str(lval[0]).rjust(2)
            for i in range(1,ny):
                if mod(i,40)==0: s+='\n'
                s+=' '+str(lval[0]).rjust(2) 
        return s+'\n' 

    def testCondition(self,cond):
        """ test if the condition is satisfied"""
        return self.core.testCondition('Modflow',cond)
        
    #************************ file for transient data *************************************
    def writeTransientFile(self,core,line,ext):
        """this method write files that have point location (wells, variable head)
        which are transient but can be permanent for wells"""
        if core.dictype['Modflow'][line][0]=='array':
            self.writeTransientArray(core,line,ext)
        else :
            self.writeTransientZones(core,line,ext)
            
    def writeTransientArray(self,core,line,ext):
        f1=open(self.fullPath +'.'+ ext,'w')
        larr = core.dicarray['Modflow'][line];#print 'mfw215',line,len(shape(larr))
        nvar = 1
        if len(shape(larr))==4: # several variables
            nvar = len(larr)
            llay,lrow,lcol = where(larr[0]!=0);#print llay,lrow,lcol
        elif len(shape(larr))==3:
            llay,lrow,lcol = where(mat!=0);#print llay,lrow,lcol
        else :
            lrow,lcol = where(mat!=0);llay=[0]*len(lrow)
        larr = array(larr,ndmin=4)
        s = str(len(lrow)) +'\n' # we suppose here just one period
        for ip in range(self.nper):
            s += str(len(lcol)) +'\n'
            for i in range(len(lrow)):
                lay,r,c = llay[i],lrow[i],lcol[i];#print lay,r,c
                s += ' %9i %9i %9i '%(lay+1,r+1,c+1)
                for iv in range(nvar): s += ' %9.2e'%larr[iv,lay,r,c]
                s += '\n'
        f1.write(s);f1.close()

    def writeTransientZones(self,core,line,ext):
        f1=open(self.fullPath +'.'+ ext,'w')
        zlist = self.ttable[line]
        nper,nzones = shape(zlist);#print 'mfw trans nz',line,nper,nzones
        nper -=1 # there is one period less than times 
          
        lpts, k, npts,zvar,lindx = [],[],0,[],[] # OA 3/3/20 added lindx
        for iz in range(nzones): #creates a list of points for each zone
            lpts.append([])
            ilay,irow,icol,zvect = self.xyzone2Mflow(core,line,iz)#OA 25/4/19,zvect
            zvar.append(zvect)
            if ilay == None: break
            npts += len(irow)
            if core.addin.mesh != None : lindx.extend(list(irow))# OA 3/3/20                                                                           
            for i in range(len(irow)):
                if core.addin.mesh == None: ## regular grid # OA modif 4/2/19
                    lpts[iz].append(str(ilay[i]+1).rjust(9)+' '+str(irow[i]+1).rjust(9)+' '+\
                       str(icol[i]+1).rjust(9))
                else : #unstruct grid irow is the node number
                    lpts[iz].append(str(irow[i]+1).rjust(9))
            #print 'mfw transt',iz,ilay,irow,ir2,lpts
            if ext=='wel': 
                k.append(self.getPermScaled(ilay,irow,icol))
        indx = argsort(lindx)    # OA added 19/4/20        
        buff = ' %9i' %npts;#print(line,zlist)
        if ext == 'wel': buff += ' 31'#EV 25/02/20 it was 90
        #buff += ' 31'
        if ext in self.usgTrans.keys():  #OA 4/3/20
            nspec = len(self.usgTrans[ext][0,0].split())
            for i in range(nspec): buff += ' AUX C%02i' %(i+1)
            buff += '\n'
        else : buff += '\n'
        #print 'mfw transient',nper
        for ip in range(nper): # get each period
            if ext == 'wel': buff += '%9i  0  0\n' %npts
            else: buff +=  '%9i \n' %npts
            flgTr, buf1 = False,[]
            if len(unique(zlist[:,iz]))>1 : flagTr = True
            for iz in range(nzones): # and each zones
                val = zlist[ip,iz] # the value of the variable for the period
                if ext in self.usgTrans.keys(): 
                    soption = self.usgTrans[ext][ip,iz]
                else : soption =''
                a = core.diczone['Modflow'].getValue(line,'value',iz)
                if '$' in a: vparms = a.split('$')[1]
                vnext = zlist[min(ip+1,nper-1),iz] #OA 7/8/17 pb of Chd
                npz = len(lpts[iz])
                for pt in range(npz): # for each zone the list of points
                    if len(unique(zvar[iz]))>1:  # OA 25/4/19 for polyV
                        zbase = float(zvar[iz][pt])
                    else : 
                        zbase = 0 # OA 25/4/19
                    if ext=='wel': 
                        s1=lpts[iz][pt]+' %9.3e '%(float(val)*k[iz][pt]) #EV 20/02/19
                    elif ext=='chd': 
                        s1=lpts[iz][pt]+' %9.3e %9.3e ' %(float(val)+zbase,float(vnext)+zbase)
                    elif ext=='drn': # elevation adding zbase (polyV), then cond.
                        v1,v2 = val.split();#print(iz,pt,zbase,v1,v2)
                        s1=lpts[iz][pt]+' %9.3e %9.3e ' %(float(v1)+zbase,float(v2))
                    elif ext=='ghb':
                        a,cond = vparms.split(); #conductance steady (time & space)
                        if flgTr: hd = float(val)
                        else : hd = float(val.split()[0])
                        s1=lpts[iz][pt]+' %9.3e %9.3e '%(hd+zbase,float(cond))
                    elif ext=='riv':
                        a,cond,botm = vparms.split()
                        if flgTr: stage = float(val)
                        else : stage = float(val.split()[0])
                        s1=lpts[iz][pt]+' %9.3e %9.3e %9.3e ' %(stage,float(cond),float(botm)+zbase)
                    buf1.append(s1+soption)
            if len(lindx)>0: buff += '\n'.join(array(buf1)[indx]) + '\n'
            else : buff += '\n'.join(buf1) + '\n'
        f1.write(buff)
        f1.close()
        
    def xyzone2Mflow(self,core,line,iz):
        """returns a list of layers, rows and cols from zones that will fit to modflow
        standards"""
        if line == 'pcb.2':modName='MfUsgTrans' # OA 2/3/20
        else : modName = 'Modflow'                                                                                            
        dicz = core.diczone[modName].dic[line] # OA /3/20
        coo = dicz['coords'][iz]
        if coo != '': xy = coo
        if xy == '': return None,None,None
        if len(xy[0])==2: x,y = list(zip(*xy));z=x*1 # OA 25/4/19 for 3 coords
        else : x,y,z = list(zip(*xy))
        imed = core.diczone[modName].getMediaList(line,iz) # OA 4/3/20
        ilay = media2layers(core,imed)
        dm = core.addin.getDim()
        if core.addin.mesh == None: # OA 4/3/20
            icol,irow,zmat = zone2index(core,x,y,z) # OA 25/4/19
            nx,ny,xvect,yvect = getXYvects(core)
            if isclosed(core,x,y) : 
                irow,icol = where(fillZone(nx,ny,icol,irow,zmat)); # OA 27/4/19 a replaced by zmat
            n0 = len(icol)
            if dm in ['3D','2D']:
                icol, irow, zmat = list(icol)*len(ilay),list(irow)*len(ilay),list(zmat)*len(ilay) # OA 27/4/19 added zmat
                ilay = list(ilay)*n0;ilay.sort()
            if dm in ['Xsection','Radial']:
                irow1=[0]*len(irow);ilay=[ny-x-1 for x in irow]
            else : 
                irow1=[ny-x-1 for x in irow]
        else : # usg
            irow = where(zmesh(core,dicz,0,iz)==1)[0] # OA 22/2/20
            n0,irow1 = len(irow),[]
            ncell_lay = core.addin.mfU.getNumber('elements')
            for il in ilay: irow1.extend(list(irow+il*ncell_lay))
            ilay,icol = list(ilay)*n0,None
#        if core.addin.mesh !=None and core.getValueFromName('Modflow','MshType')<1:
#            irow1 = array(irow1)*nx+array(icol) # uses square grid ref to connect to unstrcut rect                                                                                  
        if len(xy[0]) ==2: zmat = 0 # OA 25/4/19 return 0 if not polyV
        return ilay,irow1,icol,zmat

    def getPermScaled(self,ilay,irow,icol):
        """return the permeability for a list of layer, col rows scaled by the
        sum of permeability for this list"""
        K = self.core.getValueLong('Modflow','lpf.8',0)
        #print('mfw l400',shape(K),ilay,irow,icol)
        grd = self.core.addin.getFullGrid()
        dx=grd['dx'];dy=grd['dy'];ny=grd['ny']
        zb = self.core.Zblock
        thick = zb[:-1]-zb[1:]
        ka=ones(len(ilay))*0.;#print 'mfi permsc',shape(ka),shape(K),ilay,irow,icol
        mh = self.core.addin.mesh
        for i in range(len(ilay)):
            if self.core.addin.getDim() in ['Xsection','Radial']: 
                vol=dx[icol[i]]*dy[ny-ilay[i]-1]
                ka[i]=K[ilay[i],irow[i],icol[i]]*vol
            else : 
                if mh == None: # OA 19/4/20
                    vol=dx[icol[i]]*dy[irow[i]]*thick[ilay[i],irow[i],icol[i]]
                    ka[i]=K[ilay[i],irow[i],icol[i]]*vol
                else : 
                    ncell = mh.getNumber('elements') # OA 3/10/18 moved from above
                    irow1 = mod(irow[i],ncell)
                    vol = mh.carea[irow1]*thick[ilay[i],irow1] # irow is the cell nb
                    ka[i] = K[ilay[i],irow1]*vol # OA modif 1/8 error in calcul
        return ka/sum(ka)
        
    #*************************** fichier MNWT multinode ***********************
    def writeMNwtFile(self,core):
        f1=open(self.fullPath +'.mnwt','w') 
        # write the general properties 
        # write the parameters for each well in mnwt.2a layer
        zones = core.diczone['Modflow'].dic['mnwt.2a']
        nzones = len(zones['name'])
        s  = '#text \n'+str(nzones)+' 0 0 \n' # iwl2cb mnwprnt set to 0
        zlist = self.ttable['mnwt.2a']
        nper,n1 = shape(zlist);#print 'mfw trans nz',line,nper,nzones
        # loss type 
        #Ltype = ['NONE','THIEM','SKIN','GENERAL','SPECIFYcwc']
        for iz in range(nzones):
            a,parms,val = zones['value'][iz].split('$')
            lparm = parms.split('\n')
            s += zones['name'][iz]+' '+lparm[0] + '\n' # well name & Nnodes
            s += lparm[1]+ ' 0 0 0 0 \n' # loss type PUMPLOC Qlimit PPFLAG PUMPCAP all set to 0
            if lparm[1] != 0:
                s += ' '.join(lparm[2:6])+'\n' # Rw, or Rw Kskin, or ...for the whole well
            s += lparm[6]+' '  # layer nb or ztop
            if int(lparm[0])<0 : s+= lparm[7]+' '  # zbott if present
            #coo = zones['coords'][iz]
            #x,y = zip(*coo); z=x*1
            #icol,irow,a = zone2index(core,x,y,z)
            ilay,irow,icol,zvect = self.xyzone2Mflow(core,'mnwt.2a',iz) # EV 22/07/2019
            ir,ic = irow[0],icol[0]
            s += str(icol[0]+1)+' '+str(irow[0]+1)+'\n' # add icol and irow : they are lists
        for ip in range(nper):
            s += str(nzones)+'\n'
            s += zones['name'][0]+' '+str(zlist[ip,0])+'\n'
        f1.write(s)
        f1.close()
        
    #*************************** fichier LMT ***********************
    def writeLmtFile(self):
        
        f1=open(self.fullPath +'.lmt','w')        
        f1.write('OUTPUT_FILE_NAME    '+self.fName+'.flo \n')
        f1.write('OUTPUT_FILE_UNIT    333 \n')
        f1.write('OUTPUT_FILE_HEADER  standard \n')
        f1.write('OUTPUT_FILE_FORMAT  unformatted \n')
        f1.close()
        
    #*************************** fichier OC ***********************
    def writeOcFile(self):

        f1=open(self.fullPath +'.oc','w')     
        s = 'HEAD SAVE UNIT 30 \n'
        if len(self.usgTrans.items())>0: s += 'CONC SAVE UNIT 101 \n'
        s += 'Compact Budget \n'
        nstp=int(self.core.getValueFromName('Modflow','NSTP'));#print 'mfwrite 334',self.nper,nstp
        if self.nper>1:
            for p in range(self.nper):
                s += 'Period %5i Step %5i \n' %(p+1,nstp)
                s += 'Save Head \n'
                s += 'Save Budget \n'  #EV 25/02/20
                if len(self.usgTrans.items())>0: 
                    s += 'Save Conc \n'  # OA 30/7/19 conc
        else : s += 'Period 1 Step 1 \nSave Head\nSave Budget \n'
        f1.write(s);f1.close()
        
    #------------------------- fonction  writevect, writemat -------------------
    def writeVecModflow(self, v,ktyp):
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

    def writeMatModflow(self, m, ktyp):
        #print 'mfw',shape(m),m
        if len(shape(m))==1: return self.writeVecModflow(m,ktyp)
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

    def writeBlockModflow(self,m,ktyp):
        #print shape(m),m
        s = ''
        if len(shape(m))==3:
            nlay,a,b=shape(m);
            for l in range(nlay):
                s += self.writeMatModflow(m[l],ktyp)
                if l<nlay-1: s += '\n'
        elif self.core.addin.mesh != None: # unstructured case write nlay vectors
            nlay,a=shape(m);
            for l in range(nlay):
                s += self.writeVecModflow(m[l],ktyp)
                if l<nlay-1: s += '\n'       
        else : 
            s = self.writeMatModflow(m,ktyp)
        return s
        
    #------------------------- fonction  write HBF -------------------
    def writeHfb(self,line,iz):       # EV 26/11/2018  
        ilay,irow,icol,zvect = self.xyzone2Mflow(self.core,line,iz) # EV 22/07/2019
        #print('ilay',ilay,'irow',irow,'icol',icol)
        hbf=[]
        for i in range(len(irow)):
            if i==0 :
                if irow[i]==irow[i+1]-1 and icol[i]==icol[i+1]: # vertical line
                    ir1=irow[i] ; ic1=icol[i]
                    ir2=ir1 ; ic2=ic1+1
                    hbf.append([ir1,ic1,ir2,ic2])
                elif irow[i]==irow[i+1] and icol[i]==icol[i+1]-1: #horizontal line
                    ir1=irow[i] ; ic1=icol[i]
                    ir2=ir1+1 ; ic2=ic1
                    hbf.append([ir1,ic1,ir2,ic2])
                else :
                    ir1=irow[i] ; ic1=icol[i]
                    ir2=ir1 ; ic2=ic1+1
                    hbf.append([ir1,ic1,ir2,ic2])
            if i!=0 :
                if irow[i]==irow[i-1]+1 and icol[i]==icol[i-1]: # vertical line
                    ir1=irow[i] ; ic1=icol[i]
                    ir2=ir1 ; ic2=ic1+1
                    hbf.append([ir1,ic1,ir2,ic2])
                elif irow[i]==irow[i-1] and icol[i]==icol[i-1]+1: #horizontal line
                    ir1=irow[i] ; ic1=icol[i]
                    ir2=ir1+1 ; ic2=ic1
                    hbf.append([ir1,ic1,ir2,ic2])
                elif irow[i]>irow[i-1] and icol[i]>icol[i-1]:   #decreasing line
                    ir1=irow[i]-1 ; ic1=icol[i]
                    ir2=irow[i] ; ic2=ic1
                    hbf.append([ir1,ic1,ir2,ic2])
                    ir1=irow[i] ; ic1=icol[i]
                    ir2=ir1 ; ic2=ic1+1
                    hbf.append([ir1,ic1,ir2,ic2])
                elif irow[i]>irow[i-1] and icol[i]<icol[i-1]: # growing line
                    ir1=irow[i]-1 ; ic1=icol[i]+1
                    ir2=irow[i] ; ic2=ic1
                    hbf.append([ir1,ic1,ir2,ic2])
                    ir1=irow[i] ; ic1=icol[i]
                    ir2=ir1 ; ic2=ic1+1
                    hbf.append([ir1,ic1,ir2,ic2])  
                elif irow[i]==irow[i-1] and icol[i]<icol[i-1]:  # horizontal line in the middle
                    ir1=irow[i] ; ic1=icol[i]+1
                    ir2=ir1+1 ; ic2=ic1
                    hbf.append([ir1,ic1,ir2,ic2])
                    ir1=irow[i] ; ic1=icol[i]
                    ir2=ir1+1 ; ic2=ic1
                    hbf.append([ir1,ic1,ir2,ic2]) 
        #print(hbf)
        return hbf     
        
""" --------------------------------------------------------------------
------------------------------------------------------------------------
---------------                 Reading Modflow data      ----------------
----------------------------------------------------------------------------
 """
class modflowReader:
    
    def __init__(self, fDir, fName):
        """ on recupere le nom du projet pour effectuer l'ouverture du projet a partir de ce nom """
        self.fDir,self.fName = fDir,fName

    def readHeadFile(self, core,iper=0):
        """ read .head file 
        in free flow Thksat from flo file must be added (not done)"""    
        nlay,ncol,nrow = self.getGeom(core)       #OA 4/3/20                                     
        if core.mfUnstruct and core.getValueFromName('Modflow','MshType')>0:#OA 4/3/20   
            nlay,ncell = getNlayers(core),core.addin.mfU.getNumber('elements')
            hd=zeros((nlay,ncell));#print('mfw 491', shape(hd))
        else :
            ncell = ncol*nrow
            hd=zeros((nlay,nrow,ncol))
        try : f1 = open(self.fDir+os.sep+self.fName+'.head','rb')
        except IOError: return None
        blok=44+ncell*4; # v210 60
        for il in range(nlay):
            f1.seek(iper*nlay*blok+blok*il+44) #vpmwin
            data = arr2('f')
            data.fromfile(f1,ncell)
            if core.mfUnstruct  and core.getValueFromName('Modflow','MshType')>0: 
                hd[il] = data
            else : 
                m = reshape(data,(nrow,ncol)) #
                hd[il] = m[::-1] #=1::=1
        f1.close()  
        #modify the head if free and in 3D
        # if core.addin.getDim() in ['3D','Xsection','Radial']:   
        #     if core.addin.getModelType()=='Unconfined':  # OA 2/10/19 fre-> unconfined
        #         hd = self.getHeadFree(core,hd)
        return hd
        
    def getGeom(self,core):
        grd = core.addin.getFullGrid()
        ncol, nrow = grd['nx'], grd['ny']
        nlay=getNlayers(core);#print iper, nlay,ncol,nrow
        if core.addin.getDim() in ['Xsection','Radial']:
            nlay=nrow;nrow=1
        return nlay,ncol,nrow
        
    def getThksat(self,core,iper=0):
        try : f1 = open(self.fDir+os.sep+self.fName+'.flo','rb')
        except IOError : return None
        ncol,nrow,nlay,blok,part=self.getPart() 
        f1.seek(11+36+iper*part+36 )
        data = arr2('f')
        data.fromfile(f1,nlay*ncol*nrow)        
        thksat = reshape(data,(nlay,nrow,ncol))  
        thksat[thksat<0]=0
        thksat[thksat>1e5]=0;#print 'mfw 409 thksat',thksat
        return thksat
        
    def getHeadFree(self,core,head):
        """the head in the layer where it is 0 must be replaced by the one
        of the layer below"""
        hd1 = head[1:]*1;#print 'mfw 414 hd',head[:-1]==0
        hdry = core.dicval['Modflow']['lpf.1'][1]
        hd2 = hd1*(abs(head[:-1]-hdry)<1e-3)
        hd3 = head*1;hd3[:-1]=head[:-1]+hd2;#print 'mfw 414 hd',head,hd1,hd2,hd3
        return hd3

    def readFloFile(self, core,iper=0):
        """ read flo file and gives back Darcy velocities"""
        grd = core.addin.getFullGrid()
        dx, dy = grd['dx'], grd['dy'];
        dxm,dym=meshgrid(dx,dy)
        thick = self.getThickness(core,iper);
        try : f1 = open(self.fDir+os.sep+self.fName+'.flo','rb')
        except IOError : return None
        ncol,nrow,nlay,blok,part=self.getPart()
        if core.addin.getDim() in ['Xsection','Radial']:
            dxm=array(dx);thick = reshape(dym,(nlay,1,ncol));dym=1.;
        l0=11;#36 if for 5 keywords size 4 plus header 16 char
        pos = l0+36+iper*part+blok+36       
        f1.seek(pos);data = arr2('f');data.fromfile(f1,nlay*ncol*nrow)        
        vx = reshape(data,(nlay,nrow,ncol))
        # trouver position des vitesses y
        pos = l0+36+iper*part+blok*2+36
        f1.seek(pos);data = arr2('f');data.fromfile(f1,nlay*ncol*nrow)       
        vy = reshape(data,(nlay,nrow,ncol));bal=0.0
        # trouver position des vitesses z (nexiste que si plus d'un layer)
        bal=0.
        if nlay>1: # add vz
            nb=2
            if nrow>1: nb=3 #presence of y velo
            pos = l0+36+iper*part+blok*nb+36
            f1.seek(pos);data = arr2('f');data.fromfile(f1,nlay*ncol*nrow)       
            m0 = reshape(data,(nlay,nrow,ncol));vz=m0*0.
            for l in range(nlay): vz[l] = m0[l]/dxm/dym
            vz=concatenate([vz[:1,:,:],vz],axis=0)
        f1.close();
        # rows are ordered from top to bottom in modflow so invert them here
        if core.addin.getDim() in ['Xsection','Radial']:
            vx=vx[::-1,:,:]*1;#vz=vz[::-1,:,:]*1
        vx=vx[:,::-1,:]/dym/thick;vy=-vy[:,::-1,:]/dxm/thick; #retourner les vecteurs
        # as vx start from right of 1st cell we need to add one in the first col
        vx=concatenate([vx[:,:,:1],vx],axis=2);vx[:,:,-1]=vx[:,:,-2]
        # ssame for vy start at lower face, which is the last one now (inversion)
        vy=concatenate([vy,vy[:,-1:,:]],axis=1)
        # seems that vy is surrounded by 0
        vy[:,:,0]=vy[:,:,1];vy[:,:,-1]=vy[:,:,-2];vy[:,0,:]=vy[:,1,:]
        #print 'mfred l 436',shape(vx),shape(vy)
        if nlay>1 : return vx,vy,-vz
        else : return vx,vy,None
        
    def getThickness(self,core,iper): 
        #if type(iper)==type([5]): iper=iper[0] # takes the thick only for the 1st tstep  #EV 23/03/20 
        zb = core.Zblock
        thk = zb[:-1,:,:]-zb[1:,:,:]
        thkMat=array([thk]*len(iper))
        dim = core.addin.getDim()
        if dim in ['Xsection','Radial']:
            grd = core.addin.getFullGrid()
            nx,ny,dy = grd['nx'],grd['ny'],grd['dy'];
            ep =float(core.getValueFromName('Modflow','TOP'))-float(core.getValueFromName('Modflow','BOTM'))
            thk=ones((ny,1,nx))*ep
            thkMat=array([thk]*len(iper))  #EV 23/03/20 
            return thkMat
        if core.addin.getModelType()=='Unconfined': 
            #print('ok')
            thkMat=[]
            for i in range(len(iper)):  #EV 23/03/20 
                hd=self.readHeadFile(core,iper[i])
                if dim =='3D': hd=hd[0]
                thk[0,:,:]=hd-zb[1,:,:]
                thkMat.append(array(thk))
        return thkMat
    
    def getThicknessZone(self,core,iper,layers,ix,iy):
        #print('lay',layers)
        dim = core.addin.getDim()
        thm = self.getThickness(core,iper);#print 'mflread 474',shape(thm),thm # only 2D up to now
        thMat=zeros((len(iper),len(layers)))
        ny=len(thm[0][0]) #EV 01/04/20
        for t in range(len(iper)):  #EV 23/03/20 
            th =[] ; 
            if dim in ['Xsection','Radial']:
                for i in range(len(ix)):
                    th.append(thm[t][iy[i],0,ix[i]]) # revert for different orientation of layers
            else :
                for i in range(len(ix)):
                    th.append(thm[t][layers[i],(ny-iy[i]-1),ix[i]])  #EV 01/04/20
            thMat[t,:]=th  # OA 11/4/20 desindented
        #print('thMat',thMat)
        return thMat

    def getLocalTransientV(self,core,infos,thick,ilay,irow,icol,iper):
        """a method to get the darcy velocity at a given location and a given period from
        the flo file. 
        returns two velocities at the cell boundaries in each direction
        only in x,y but for the correct layer"""
        grd = core.addin.getFullGrid()
        dx,dy,nx,ny = grd['dx'],grd['dy'],grd['nx'],grd['ny']
        try : f1 = open(self.fDir+os.sep+self.fName+'.flo','rb')
        except IOError : return None
        ncol,nrow,nlay,blok,part=infos;#print 'mfw 495',ncol,nrow,nlay
        l0=11;
        #print irow,icol,iper
        pos1 = l0+36+iper*part+blok+36
        pos2 = ilay*ncol*nrow*4+(ny-irow-1)*ncol*4+(icol-1)*4 # ny because modflow and ipht3d ordered differently   
        f1.seek(pos1+pos2);vx0 = arr2('f');vx0.fromfile(f1,2)        
        f1.seek(pos1+blok+pos2-(nx-1)*4);vy0 = arr2('f');vy0.fromfile(f1,nx+1);vy0=[vy0[-1],vy0[0]]      
        if nlay>1:
            nb = 1
            if nrow>1: nb=2
            f1.seek(pos1+blok*nb+pos2);vz0 = arr2('f');vz0.fromfile(f1,2)       
        f1.close();
        if core.addin.getDim() in ['Xsection','Radial']:
            return [vx0/dy[irow]/1.,0-vz0/dx[icol]/1.]
        vx1 = vx0/dy[irow]/thick[ilay,irow,icol]
        vy1 = vy0/dx[icol]/thick[ilay,irow,icol]
        #if nlay>1 : return [vx1,vy1,0-vz0/dy[irow]/dx[icol]]
        #else : 
        return [vx1,-vy1]
        
    def readWcontent(self,core,iper=0):
        try : f1 = open(self.fDir+os.sep+self.fName+'.flo','rb')
        except IOError : return None
        # this is a file produced by uzf
        ncol,nrow,nlay,blok,part = self.getPart2(f1)
        title = 11+21*4
        pos = title+iper*part+5*4+16     # OA 9/3/19
        f1.seek(pos);data = arr2('f');data.fromfile(f1,nlay*ncol*nrow)   
        uzsat = reshape(data,(nlay,nrow,ncol))
        #pos += blok*4*4 # Thksat is 4 bloc further
        #f1.seek(pos);data = arr2('f');data.fromfile(f1,nlay*ncol*nrow)   
        #thks = reshape(data,(nlay,nrow,ncol))
        return uzsat # OA 9/3/19 removed ::-1
        
    def getPart2(self,f1):
        # for UZF variables Wcontent, UzFlux, UzSto, GwOut,ThKsat, Qxx, Qzz, Sto, Cnh
        nvar = 8 # CNH not included, it has different shape
        # simplified version : no wells, no rech... should be included into getpart
        f1.seek(11);data = arr2('i');data.fromfile(f1,26);
        ncol,nrow,nlay = data[-3:]
        ncnh = data[6];l2=5*4+16 # OA 9/3/19 5 values iper isubper..+title
        blok=l2+ncol*nrow*nlay*4;
        part = nvar*blok+l2+4+ncnh*4*4
        return ncol,nrow,nlay,blok,part
        
    def getPtObs(self,core,irow,icol,ilay,iper,typ):
        """for typ=flux return fluxes for a series of obs cell
        for typ=head return the heads
        iper is a list of periods indices"""
        try : f1 = open(self.fDir+os.sep+self.fName+'.flo','rb')
        except IOError : return None
        nper=len(iper);
        if core.addin.getDim() in ['Xsection','Radial']:
            ilay=irow*1;irow=[0]*len(ilay) #[-1::-1]
        #print 'mfw getpt',typ,irow,icol
        if typ=='Head': 
            return self.getHeadPtObs(core,irow,icol,ilay,iper)
        elif typ=='Wcontent':
            return self.getWcontentPtObs(f1,core,irow,icol,ilay,iper)
        ncol,nrow,nlay,blok,part = self.getPart()
        blok2 = ncol*nrow*4
        l0=11
        qx= zeros((nper,len(irow)));qy=qx*0.
        for ip in range(nper):
            posx = l0+36+iper[ip]*part+blok+5*4+16
            posy = posx+blok
            for i in range(len(irow)):
                pos2 = blok2*ilay[i]+irow[i]*ncol*4+icol[i]*4
                f1.seek(posx+pos2);data = arr2('f');data.fromfile(f1,1)
                qx[ip,i] = float(data[0])
                f1.seek(posy+pos2);data = arr2('f');data.fromfile(f1,1)        
                qy[ip,i] = float(data[0])
        return qx,qy

    def getHeadPtObs(self,core,irow,icol,ilay,iper):
        try : f1 = open(self.fDir+os.sep+self.fName+'.head','rb')
        except IOError: return None
        grd = core.addin.getFullGrid()
        ncol, nrow = grd['nx'], grd['ny']
        nlay = getNlayers(core);
        if core.addin.getDim() in ['Xsection','Radial']:
            nlay=nrow*1;nrow=1;# already done above ilay=irow*1;irow=[0]*len(ilay)
        nper=len(iper)
        hd=zeros((nper,len(irow)))
        f1.seek(32);data=arr2('i');
        blok=44+nrow*ncol*4;#print 'mfw gethedpt ilay,icol,row',ilay,icol,irow
        for ip in range(nper):
            for i in range(len(irow)):
                pos=44+iper[ip]*nlay*blok+ilay[i]*blok+irow[i]*ncol*4+icol[i]*4;
                f1.seek(pos)
                data = arr2('f');data.fromfile(f1,1);#print iper[ip],irow[i],icol[i],data
                hd[ip,i]=float(data[0])
        f1.close()
        return hd
        
    def getWcontentPtObs(self,f1,core,irow,icol,ilay,iper):
        grd = core.addin.getFullGrid()
        ncol, nrow = grd['nx'], grd['ny']
        nlay = getNlayers(core);
        #print nlay,ncol,nrow,ilay,icol,irow
        if core.addin.getDim() in ['Xsection','Radial']:
            nlay=nrow*1;nrow=1;ilay=irow*1;irow=[0]*len(ilay)
        nper=len(iper)
        ncol,nrow,nlay,part = self.getPart2(f1)
        wc = zeros((nper,len(irow)))
        l0 = 11+26*4
        for ip in range(nper):
            pos = l0+iper[ip]*part+5*4+16
            for i in range(len(irow)):
                pos2=ilay[i]*nrow*ncol*4+irow[i]*ncol*4+icol[i]*4;
                f1.seek(pos+pos2)
                data = arr2('f');data.fromfile(f1,1);#print iper[ip],irow[i],icol[i],data
                wc[ip,i]=float(data[0])
        return wc
        
    def getPart(self):
        f1 = open(self.fDir+os.sep+self.fName+'.flo','rb')
        l0=11;f1.seek(l0);data = arr2('i');data.fromfile(f1,14);
        mwel,mdrn,mrch,mevt,mriv,mghb,mchd,mss,nper,kp,kst,ncol,nrow,nlay = data
        blok=5*4+16+ncol*nrow*nlay*4;l2=5*4+16+4;l0+=36
        nwel=0;mss=1; # even in SS it seems to have sto?
        ichd,iwel,irch,idrn=sign(mchd),sign(mwel),sign(mrch),sign(mdrn)
        ievt,iriv,ighb=sign(mevt),sign(mriv),sign(mghb)
        icol,irow,ilay=sign(ncol-1),sign(nrow-1),sign(nlay-1);
        part=blok+icol*blok+irow*blok+ilay*blok+mss*blok
        if nper>1:
            f1.seek(l0+part+l2-4);data = arr2('i');data.fromfile(f1,1)
            nchd = data[0]
            part += l2+nchd*16
        if iwel>0:
            f1.seek(l0+part+l2-4);data = arr2('i');data.fromfile(f1,1)
            nwel = data[0];part += l2+nwel*16
        if idrn>0:
            f1.seek(l0+part+l2-4);data = arr2('i');data.fromfile(f1,1)
            ndrn =  data[0]      
            part += l2+ndrn*16
        if irch>0: part += l2-4+ncol*nrow*4*2
        if ievt>0: part += l2-4+ncol*nrow*4*2
        if iriv>0:
            f1.seek(l0+part+l2-4);data = arr2('i');data.fromfile(f1,1)
            nriv = data[0];part += l2+nriv*16
        if ighb>0:
            f1.seek(l0+part+l2-4);data = arr2('i');data.fromfile(f1,1)
            nghb = data[0];part += l2+nghb*16
        return ncol,nrow,nlay,blok,part
        
    def readFlowMesh(self,core,mesh):
        '''a simplified reader of the budget file
        if we assume flux enters on face 1 and flows out on face 2 and 3, the velocity vector
        starts from the point on face 1 corresponding to the proportion of flux
        leaving face 2 and 3
        and stops at the instersection point between face 2 and 3
        also works if flux is entering on two cells and leaving on one'''
        nlay = getNmedia(core) ### !!! nlayer = nmedia
        try : f1 = open(self.fDir+os.sep+self.fName+'.budget','rb')
        except IOError: return None
        f1.seek(0);data =arr2('i');data.fromfile(f1,2);kstep,kper = data
        f1.seek(24);data =arr2('i');data.fromfile(f1,4);nval,step,icode,imeth = data
        if icode>=0: 
            data =arr2('f');data.fromfile(f1,nval);pos = 24+4*4+nval*4
        else :
        # method =2
            data =arr2('f');data.fromfile(f1,3); delt,pertim,totim = data
            data = arr2('i'); data.fromfile(f1,1);nlist=data[0]
            data = arr2('i'); data.fromfile(f1,nlist*2) # should have the 2nd value as 'f'
            pos = 24+4*4+4*4+nlist*8+4
        # get flow faces valid only when nothing else is printed
        f1.seek(pos+4);
        data =arr2('c');data.fromfile(f1,16) # this is the 'flow ja face  ' name
        data =arr2('i');data.fromfile(f1,3); njag,step,icode = data # not sure
        if icode>=0:
            data =arr2('f');data.fromfile(f1,njag);             
        else :
            data =arr2('f');data.fromfile(f1,4); delt,pertim,totim,a = data
            data =arr2('f');data.fromfile(f1,njag); 
        ncell,ncelltot = mesh.ncell_lay,mesh.ncell
        ln = [len(mesh.cneighb[ic])+1 for ic in range(ncelltot)] # +1 because flux contain 1 more
        pos = [0]; pos.extend(cumsum(ln))
        vx, vy, vz = zeros((nlay,ncell)),zeros((nlay,ncell)),zeros((nlay,ncell))
        l1=[0,1,0,2] # a list to transform sum of faces to common point
        i3d = (ncelltot>ncell)*1
        for ic in range(ncelltot):
            il = int(floor(ic/ncell))
            ic1 = int(ic-il*ncell)
            q = array(data[pos[ic]:pos[ic+1]])
            if q[0]!=0: 
                #print(ic)
                continue # not fair but there seems to be some values where ln does not correpsond to data!!!
            q1 = q[1:ln[ic1]-i3d] # use ic1 to keep it 2d
            s0 = sign(q1)
            idx = list(where(s0==sum(s0))[0])
            if (len(q1)==2) or (len(idx)==0) or (sum(abs(s0))<=1): continue
            inod = mesh.elements[ic1,-3:]
            inod1 = r_[inod, inod[0]]
            xno,yno = mesh.nodes[inod1,1],mesh.nodes[inod1,2]
            prop = q1[idx]/sum(q1[idx])  # proportion of flux that leaves on each face
            id1 = where(s0!=sum(s0))[0][0]    #  face with only one flux
            ifa = id1-1; 
            if ifa<0: ifa=2
            iprop = idx.index(ifa) # face to search for the proportion
            xp1 = xno[id1]+prop[iprop]*(xno[id1+1]-xno[id1]) # x point 
            yp1 = yno[id1]+prop[iprop]*(yno[id1+1]-yno[id1]) # y point
            ip2 = l1[sum(idx)] # pt commom to faces 
            xp2,yp2 = xno[ip2],yno[ip2]
            v0, v1 = (xp2-xp1)*s0[id1], (yp2-yp1)*s0[id1] # vx, vy but no module
            d = sqrt((xp2-xp1)**2+(yp2-yp1)**2)
            vangle = arctan(v1/v0)
            v = abs(q1[id1]*cos(vangle+pi/2-mesh.angles[ic1][id1])/mesh.fahl[ic][id1])
            vx[il,ic1],vy[il,ic1] = v0/d*v, v1/d*v
            if nlay>1 :
                vz[il,ic1] = q[-1]/mesh.carea[ic1]
        return vx,vy,vz

        
