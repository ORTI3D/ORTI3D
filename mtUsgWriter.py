# -*- coding: cp1252 -*-
from array import array as arr2
import os,time
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
        usgTrans = {'active':True}
        self.ttable = self.core.makeTtable();#print 'mtpht ttable',self.ttable
        self.dim = self.core.addin.getDim()
        self.core.updateDicts()
        usgTrans['nam'] = self.writeNamFile(opt)
        self.mfloW.writeModflowFiles(self.core,usgTrans=usgTrans)
        tlist = array(self.ttable['tlist'])
        self.per = tlist[1:]-tlist[:-1]
        self.nper = len(self.per)#;print('writempht l.26',self.nper)
        mcomp,ncomp,gcomp = listEsp['mcomp'],listEsp['ncomp'],listEsp['gcomp']
        self.nesp = ncomp
        nkim = len(listEsp['kim'])
        self.writeFiles(opt)
        if opt=='Pht3d':
            self.writePhFile(self.core,listEsp,parmk)
            self.writePhreeqc(self.core,listEsp);
        if 'pcb.2' in list(self.core.diczone['MfUsgTrans'].dic.keys()):
            self.writePcbFile(self.core)
        return 
        
    def writeNamFile(self,opt):
        s='BCT  41 '+self.fName+'.bct\n'
        if 'pcb.2' in list(self.core.diczone['MfUsgTrans'].dic.keys()):
            s += 'PCB  42 '+self.fName+'.pcb\n'
        if opt=='Pht3d':
            s += ' PHC  64    Pht3d_ph.dat\n'
        s += 'DATA(BINARY) 101  '+self.fName+'.conc\n'
        return s

    #*********************** generic file writer ****************
    def writeFiles(self,opt):
        """to write all modflow usg transport files.
        reads the keyword file and prints all keywords by types : param (0D)
        vector (1D) array (2D). types are found by (dim1,dim2).."""
        lexceptions, s = [],''
        llist=self.Mkey.groups['BCT'];#print n1,name
        for ll in llist:
            cond=self.Mkey.lines[ll]['cond'];#print('mtw 50',ll)
            if self.testCondition(cond)==False : continue
            kwlist=self.Mkey.lines[ll]['kw']
            ktyp=self.Mkey.lines[ll]['type']
            lval=self.core.dicval['MfUsgTrans'][ll];
            if ktyp[0] in ['vecint','vecfloat','arrint','arrfloat']:
                s += self.writeArray(opt,ll,0)
            elif ktyp[0]=='title': # case of a title line
                s += '#'+str(ktyp[0]).rjust(10)+'\n'
            else : # classical keywords
                for ik in range(len(kwlist)):
                    if ik<len(lval): s += str(lval[ik]).rjust(10)
                    else : s += '0'.rjust(10)
                s += '\n'
        f1=open(self.fDir+os.sep+self.fName +'.bct','w')
        f1.write(s);f1.close()
            #print grp+' written'

    def testCondition(self,cond):
        """ test if the condition is satisfied"""
        return self.core.testCondition('MfUsgTrans',cond)
        
    def writeArray(self,opt,line,ik):
        """writes arrays, need specific treatment for btn concentrations if pht3d
        and also for react modules of MfUsgTrans"""
        grd = self.core.addin.getFullGrid()
        dx,ny = array(grd['dx']),int(grd['ny'])
        if (opt=='Pht3d') and (line in ['btn.13','rct.2c']):
            s = ''
            self.Conc, self.Names = self.getConcInit('main','ph.3',iper=0);
            nspec = len(self.Names)
            for i in range(nspec):
                s += self.formatBlockMusg(self.Conc[i],self.Names[i])
        else: # normal print of array
            arr = self.core.getValueLong('MfUsgTrans',line,ik);
            s = self.formatBlockMusg(arr,line)
        return s 
            
    #************************ file for transient data ******************
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
        

    #  '''''''''''''''' fonction writeBlockMusg '''''''''''''''''''''''''''
    def writeVecModflow(self, v,line):
        #print shape(v),amin(v),amax(v)
        l=len(v);ln=3;s=''
        a=str(type(v[0]))
        if 'int' in a: typ='I'
        else : typ='G'
        if amin(v)==amax(v):
            if typ=='I': s += 'CONSTANT     %9i  \n' %amin(v)
            else : s += 'CONSTANT     %9.5e  \n' %amin(v)
            return s
        if typ=='I': fmt='1    ('+str(l)+'I'+str(ln)
        else : fmt='0    ('+str(l)+'G12.4'           
        s += 'INTERNAL     '+fmt+')     3 \n'
        
        if typ=='I': fmt='%'+str(ln)+'i'
        else : fmt='%+11.4e '            

        for i in range(l):
            s += fmt %v[i]
        return s+'\n'


    def formatBlockMusg(self,m,line):
        s=''
        if len(shape(m))==2:
            nlay,a = shape(m)
            for l in range(nlay):
                s += self.writeVecModflow(m[l],line)
                if l<nlay-1: s += '\n'
        else :
            s += self.writeVecModflow(m,line)            
        return s
              
class mtUsgReader:

    """ this is the reder of UCN files """
    def __init__(self, fDir, fName):
        self.fDir = fDir
        self.fName = fName

    def readUCN(self,core,opt,iper,iesp,specname=''): 
        """ read .conc file, here opt, iesp, specname are not used
        in free flow Thksat from flo file must be added (not done)"""    
        if core.mfUnstruct:
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
            f1.seek(iper*nlay*blok+blok*il+44) #vpmwin
            data = arr2('f')
            data.fromfile(f1,ncell)
            if core.mfUnstruct: 
                cnc[il] = data
            else : 
                m = reshape(data,(nrow,ncol)) #
                cnc[il] = m[::-1] #=1::=1
        f1.close()  
        return cnc
