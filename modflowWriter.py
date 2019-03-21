#import scipy.io.array_import as IOAscii
from array import array as arr2
import os
from .modflowKeywords import Mf
from .geometry import *
from .timeperiod import *
from .modflowUNS import *

class modflowWriter:
    
    def __init__(self, core,fDir, fName):
        self.core = core
        self.fDir,self.fName,self.Fkey = fDir,fName,Mf()
        self.fullPath = fDir+os.sep+fName;#print self.fullPath

    def writeModflowFiles(self, core):
        nbfor,recharge = 0,0
        self.ttable = core.makeTtable();#print 'mfw',self.ttable
        tlist = array(self.ttable['tlist'])
        self.per = tlist[1:]-tlist[:-1]
        self.nper = len(self.per)
        self.core.Zblock = makeZblock(self.core);#print 'mfw,zb',shape(self.core.Zblock)
        self.dicval = core.dicval['Modflow']
        self.core.setValueFromName('Modflow','NPER',self.nper) # nper may be wrong if transietn zones
        self.radfact=1.
        if self.core.addin.getDim() in ['Radial','Xsection']: self.setRadial()
        self.writeNamFile()
        self.writeFiles()
        self.writeLmtFile()
        self.writeOcFile()
        #print 'mfwrite',self.ttable
        if 'bas.5' in list(self.ttable['Transient'].keys()): #in bas.5 initial head there are variable zones
            if self.ttable['Transient']['bas.5']:
                self.writeTransientFile(core,'bas.5','chd')
        for n in ['wel','drn','riv','ghb']:
            if self.core.diczone['Modflow'].getNbZones(n+'.1')>0 or self.core.dictype['Modflow'][n+'.1'][0]=='array':
                self.writeTransientFile(core,n+'.1',n)
        if self.core.diczone['Modflow'].getNbZones('mnwt.2a')>0:
            self.writeMNwtFile(core)
            
    def setRadial(self):
        g = self.core.addin.getFullGrid()
        self.core.setValueFromName('Modflow','NLAY',g['ny'])
        self.core.setValueFromName('Modflow','NROW',1) 
#        for n in ['LAYCBD','LAYTYP','LAYAVG','LAYVKA','LAYWET']:
#            lc = self.core.getValueFromName('Modflow',n);print 'mfw setr',lc,type(lc)
#            self.core.setValueFromName('Modflow',n,[lc]*int(g['ny']))

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
            if mod in ['RCH','WEL','CHD','MNWT']: continue
            f1.write(mod+'  '+str(i+10)+'  ' + self.fName + '.'+mod.lower()+'\n')
        f1.write('OC       26    '+ self.fName + '.oc\n')
        f1.write('DATA(BINARY)     30        ' + self.fName + '.head\n')
        f1.write('DATA(BINARY)     31        ' + self.fName + '.budget\n')

        if self.core.mfUnstruct==False:
            f1.write('LMT6     28    '+ self.fName + '.lmt\n')
        r0=self.core.getValueFromName('Modflow','RECH')
        if 'RCH' in lmod:
            if r0>0 or self.core.diczone['Modflow'].getNbZones('rch.2')>0:
                f1.write('RCH     34     ' + self.fName + '.rch\n')
        e0=self.core.getValueFromName('Modflow','EVT')
        if 'EVT' in lmod:
            if e0>0 or self.core.diczone['Modflow'].getNbZones('evt.2')>0:
                f1.write('EVT     35     ' + self.fName + '.evt\n')
        if self.core.diczone['Modflow'].getNbZones('wel.1')>0:
            f1.write('WEL     36     ' + self.fName + '.wel\n')
        if 'bas.5' in list(self.ttable['Transient'].keys()):
            if self.ttable['Transient']['bas.5']:
                f1.write('CHD     37     ' + self.fName + '.chd\n')
        if self.core.diczone['Modflow'].getNbZones('mnwt.2a')>0:
            f1.write('MNW2     39     ' + self.fName + '.mnwt\n')
        f1.close()

    #*********************** generic file writer ****************
    def writeFiles(self):
        """to write all modflow file.
        reads the keyword file and prints all keywords by types : param (0D)
        vector (1D) array (2D). types are found by (dim1,dim2).."""
        lexceptions=['dis.1','dis.4','dis.5','dis.8','rch.2','uzf.9','uzf.10','lpf.2']#EV 06/11 added 'dis.1'
        lexceptions.extend(['disu.2','disu.3','disu.6','disu.9']) 
        lexceptions.extend(['lpf.'+str(a) for a in range(8,14)]) # when wirting by layers
        lexceptions.extend(['upw.'+str(a) for a in range(7,13)])
        lexceptions.extend(['uzf.'+str(a) for a in range(2,8)])
        lexceptions.extend(['evt.'+str(a) for a in range(2,5)])
        lexceptions.extend(['sms.1a','hfb.3'])
        for grp in self.core.getUsedModulesList('Modflow'):
            if grp in ['WEL','DRN','RIV','MNWT','GHB']: continue # WEL is written by transientfile
            ext = grp
            if grp=='Solver': ext=self.solv
            f1=open(self.fullPath +'.'+ ext.lower(),'w')
            llist=self.Fkey.groups[grp];#print('llist',llist)
            for ll in llist:
                cond=self.Fkey.lines[ll]['cond'];#print 'mfw 96',ll
                if self.testCondition(cond)==False : continue
                kwlist=self.Fkey.lines[ll]['kw']
                ktyp=self.Fkey.lines[ll]['type'];#print 'mfw',ktyp
                lval=self.dicval[ll]
                if (ll in lexceptions) or (ktyp[0][:3]=='lay'):
                    self.writeExceptions(ll,kwlist,ktyp[0],f1)
                    continue
                for ik in range(len(kwlist)):
                    value=lval[ik]
                    if ktyp[ik] in ['vecint','vecfloat','arrint','arrfloat']:
                        value=self.core.getValueLong('Modflow',ll,ik);#print 'mfw 106',ll,shape(value)
                        self.writeBlockModflow(value,f1,ktyp[ik]) # OA 1/8/17
                    elif ktyp[ik]=='choice': # where there is a choice print the nb othe choice not value
                        f1.write(str(value).rjust(10))
                    elif ktyp[ik]=='title': # case of a title line
                        f1.write('#'+str(value).rjust(10))
                    else : # values with strings
                        f1.write(str(value).rjust(10));#print 'write str',ll
                f1.write('\n')
            f1.close()
            #print grp+' written'
                    
    def writeExceptions(self,line,kwlist,ktyp,f1):
        """to write some things modflow wants in a ...specific format"""
        dim = self.core.addin.getDim()
        nx,ny,a,b = getXYvects(self.core)
        if dim in ['2D','3D']: nlay = getNlayers(self.core)
        else : nlay = ny
        tlist = array(self.ttable['tlist'])
        nper = len(tlist)-1
        
        if line == 'dis.4': # pb for radial and xsection with layers
            delr = self.core.dicval['Modflow'][line]
            self.writeVecModflow(delr,f1,'vecfloat')
            f1.write('\n')
        
        if line == 'dis.5': # row height needs to be inverted
            delc = self.core.dicval['Modflow'][line]
            if dim in ['2D','3D'] : self.writeVecModflow(delc[-1::-1],f1,'vecfloat')
            elif dim=='Radial': f1.write('CONSTANT     1  ')
            elif dim=='Xsection': 
                front = self.core.getValueFromName('Modflow','TOP')
                end = self.core.getValueFromName('Modflow','BOTM')
                f1.write('CONSTANT    '+str(front-end))
            f1.write('\n')
                    
        if line == 'disu.6':
            nodelay = int(self.core.getValueFromName('Modflow','NODELAY'))
            f1.write('CONSTANT    '+str(nodelay)+'\n')
        
        if line in ['dis.8','disu.9']: # periods characteristics 4 values per period
            s1 = ''
            if line == 'disu.9':
                s1 = self.core.addin.mfU.writeDisu()
            perlen = tlist[1:]-tlist[:-1]
            lval=self.dicval[line] # contains period sze, end time, 3 things to print
            SsTr = 'SS'
            if lval[3]==1 : SsTr = 'Tr'
            for ip in range(nper):
                if ip>0: SsTr='Tr'
                s1 +=' %9.3e %9i %9.3e ' %(perlen[ip],lval[1],lval[2]) #OA 13/9/18
                s1 += SsTr.rjust(10)+'\n'
            f1.write(s1)
            
        if line=='sms.1a':
            s = self.core.dicval['Modflow'][line][0];# print 'mfw 168',s, type(s)
            if s!= 0: 
                f1.write(self.Fkey.lines[line]['kw'][i])
            
        if line in ['lpf.8','upw.7']: # writes several lines per layer
            strt = int(line[-1])
            llist = [line[:4]+str(a) for a in range(strt,strt+6)];#take four lines
            value = []
            for l2 in llist:
                cond = self.Fkey.lines[l2]['cond'];#print 'mfw cond',l2,cond,self.testCondition(cond)
                if self.testCondition(cond) == False : continue
                v0 = self.core.getValueLong('Modflow',l2,0);#print 'mfw 173', l2,v0
                value.append(v0)
            val = self.core.dicval['Modflow']['lpf.2'] # EV 23/11/2018
            ilay=getNlayersPerMedia(self.core) 
            lval1 = [[val[x]]*ilay[x] for x in range(len(ilay))]
            lval= [item for sublist in lval1 for item in sublist]
            for l in range(nlay):
                for i in range(len(value)):
                    #print 'mfw 183',l,i,shape(value[i][l])
                    if i==3 and lval[l]==0 : continue
                    self.writeMatModflow(value[i][l],f1,'arrfloat');f1.write('\n')
                    
        if line in ['rch.2','uzf.9']: # writes an array for each period with zero before
            if line=='uzf.9': line2='uzf.10'
            elif line == 'rch.2': line2=line*1
            #print 'mfw rch',self.core.ttable[line2]
            m0=0
            for iper in range(nper): # OA 14/2/2019 modif below if rech remain similar
                #if iper>0: m0 = block(self.core,'Modflow',line2,False,None,iper); #EV 21/03/2019 remove OA modif 
                m = block(self.core,'Modflow',line2,False,None,iper);
                #if sum(m0-m)==0:
                    #f1.write('    -1\n')
               # else :
                f1.write('    0\n') # 0: data are written not reused from previous
                self.writeMatModflow(m[0],f1,'arrfloat');f1.write('\n')

        if line == 'evt.2': # writes an array for each period with zero before
            for iper in range(nper):
                f1.write('    0     0    0\n') # INSURF INEVTR INEXDP INIEVT not read
                for k in range(2,5):
                    m = block(self.core,'Modflow',line[:4]+str(k),False,None,iper);#print iper,shape(m),amax(m),where(m==1)
                    self.writeMatModflow(m[0],f1,'arrfloat');f1.write('\n')
                
        if line in ['uzf.2','uzf.3','uzf.4','uzf.5','uzf.6','uzf.7']: # in uzf put only one value of these parameters
            m = self.core.getValueLong('Modflow',line,0)
            self.writeMatModflow(m[0],f1,ktyp);f1.write('\n')  # OA 7/3/19
            
        if line == 'lpf.2':
            ilay=getNlayersPerMedia(self.core) # EV 23/11/2018
            val = self.core.dicval['Modflow'][line]
            lval1 = [[val[x]]*ilay[x] for x in range(len(ilay))]
            lval= [item for sublist in lval1 for item in sublist]
            s=''
            for i in range(len(lval)):
                s+=' '+str(int(lval[i])).rjust(2) 
            f1.write(s+'\n') 
            
        #if line in ['lpf.11']

        if ktyp[:3] == 'lay': # to print laycbd and others
            lval = self.core.dicval['Modflow'][line][0] # in 3D, should be a list of values
            s=''
            if dim == '2D' : s=str(lval)
            elif dim == '3D': #presently there is only one value for all layers
                s=' '+str(lval).rjust(2)
                for i in range(1,nlay):
                    if mod(i,40)==0: s+='\n'
                    s+=' '+str(lval).rjust(2)
            else : # radial and xsection
                s=str(lval).rjust(2)
                for i in range(1,ny):
                    if mod(i,40)==0: s+='\n'
                    s+=' '+str(lval).rjust(2) 
            f1.write(s+'\n')   
        
        if line=='hfb.3': # EV 26/11/2018
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
            f1.write(' 0  0  '+str(sum(nbHfb))+'\n')
            f1.write(s)               
            
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
          
        lpts, k, npts = [],[],0 # lists of points for each zone
        for iz in range(nzones): #creates a list of points for each zone
            lpts.append([])
            ilay,irow,icol = self.xyzone2Mflow(core,line,iz)
            if ilay == None: break
            npts += len(irow)
            for i in range(len(irow)):
                if core.addin.mesh == None: # regular grid # OA modif 4/2/19
                    lpts[iz].append(str(ilay[i]+1).rjust(9)+' '+str(irow[i]+1).rjust(9)+' '+\
                       str(icol[i]+1).rjust(9))
                else : #unstruct grid irow is the node number
                    lpts[iz].append(str(irow[i]+1).rjust(9))
            #print 'mfw transt',iz,ilay,irow,ir2,lpts
            if ext=='wel': k.append(self.getPermScaled(ilay,irow,icol))
            
        buff = ' %9i   \n' %npts
        #print 'mfw transient',nper
        for ip in range(nper): # gets each period
            buff += ' %9i   \n' %npts
            for iz in range(nzones): # and each zones
                val = zlist[ip,iz] # the value of the curent variable for the period
                vnext = zlist[min(ip+1,nper-1),iz] #OA 7/8/17 pb of Chd
                npz = len(lpts[iz])
                for pt in range(npz): # for each zone the list of points
                    if ext=='wel': 
                        buff+= lpts[iz][pt]+' %+9.2e'%(float(val)*k[iz][pt])+'\n' #EV 20/02/19
                    elif ext=='chd': 
                        buff+=lpts[iz][pt]+' %+9.5g%+9.5g  \n'%(float(val),float(vnext))
                    elif ext in ['drn','ghb']:
                        v1,v2 = val.split()
                        buff+=lpts[iz][pt]+' %+9.6g %+9.6g  \n'%(float(v1),float(v2))
                    elif ext=='riv':
                        v1,v2,v3 = val.split()
                        buff+=lpts[iz][pt]+' %+9.6g %+9.6g %+9.6g  \n'%(float(v1),float(v2),float(v3))
        f1.write(buff)
        f1.close()
        
    def xyzone2Mflow(self,core,line,iz):
        """returns a list of layers, rows and cols from zones that will fit to modflow
        standards"""
        dicz = core.diczone['Modflow'].dic[line] #print zone['name'][iz]
        coo = dicz['coords'][iz]
        if coo != '': xy = coo
        if xy == '': return None,None,None
        x,y = list(zip(*xy))
        z=x*1
        imed = core.diczone['Modflow'].getMediaList(line,iz)
        ilay = media2layers(core,imed)
        dm = core.addin.getDim()
        if core.addin.mesh == None: # regular grid
            icol,irow,a = zone2index(core,x,y,z)
            nx,ny,xvect,yvect = getXYvects(core)
            if isclosed(core,x,y) : 
                irow,icol = where(fillZone(nx,ny,icol,irow,a)); # 22/3/17
            n0 = len(icol)
            if dm in ['3D','2D']:
                icol, irow = list(icol)*len(ilay),list(irow)*len(ilay)
                ilay = list(ilay)*n0;ilay.sort()
            if dm in ['Xsection','Radial']:
                irow1=[0]*len(irow);ilay=[ny-x-1 for x in irow]
            else : 
                irow1=[ny-x-1 for x in irow]
        else : # usg
            irow = zmesh(core,dicz,0,iz) # media=0 will calc it below
            n0,irow1 = len(irow),[]
            ncell_lay = core.addin.mfU.getNumber('elements')
            for il in ilay: irow1.extend(list(irow+il*ncell_lay))
            ilay,icol = list(ilay)*n0,None
        return ilay,irow1,icol

    def getPermScaled(self,ilay,irow,icol):
        """return the permeability for a list of layer, col rows scaled by the
        sum of permeability for this list"""
        K = self.core.getValueLong('Modflow','lpf.8',0)
        #print 'mfw l264',shape(K),ilay,irow,icol
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
                if mh == None: # struct grid
                    vol=dx[icol[i]]*dy[irow[i]]*thick[ilay[i],irow[i],icol[i]]
                    ka[i]=K[ilay[i],irow[i],icol[i]]*vol
                else : 
                    ncell = mh.getNumber('elements') # OA 3/10/18 moved from above
                    irow1 = mod(irow[i],ncell)
                    vol = mh.carea1[irow[i]]*thick[ilay[i],irow1] # irow is the cell nb
                    ka[i] = K[ilay[i],irow1*vol]
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
            ilay,irow,icol = self.xyzone2Mflow(core,'mnwt.2a',iz)
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
        f1.write('HEAD SAVE UNIT 30 \nCompact Budget \n')
        nstp=int(self.core.getValueFromName('Modflow','NSTP'));#print 'mfwrite 334',self.nper,nstp
        if self.nper>1:
            for p in range(self.nper):
                f1.write('Period %5i Step %5i \n' %(p+1,nstp))
                f1.write('Save Head \n')
        else : f1.write('Period 1 Step 1 \nSave Head\nSave Budget \n')
        f1.close()
        
    #------------------------- fonction  writevect, writemat -------------------
    def writeVecModflow(self, v, f1,ktyp):
        #print shape(v),amin(v),amax(v)
        l=len(v);ln=3
        if ktyp[3:]=='int': typ='I' #OA 1/8/17
        else : typ='G'
        if amin(v)==amax(v):
            f1.write('CONSTANT     %9.5e  ' %amin(v))
            return
        if typ=='I': fmt='1    ('+str(l)+'I'+str(ln)
        else : fmt='0    ('+str(l)+'G12.4'           
        f1.write('INTERNAL     '+fmt+')     3 \n' )
        
        if typ=='I': fmt='%'+str(ln)+'i'
        else : fmt='%+11.4e '            

        for i in range(l):
            f1.write(fmt %v[i])
        #f1.write('\n')

    def writeMatModflow(self, m, f1,ktyp):
        #print 'mfw',shape(m),m
        if len(shape(m))==1: return self.writeVecModflow(m,f1,ktyp)
        [l,c] = shape(m);ln=3
        if ktyp[3:]=='int': typ='I' #OA 1/8/17
        else : typ='G'
        if amin(amin(m))==amax(amax(m)):
            if typ=='I': f1.write('CONSTANT     %9i  ' %(amin(amin(m))))
            else : f1.write('CONSTANT     %9.5e  ' %(amin(amin(m))))
            return
        if typ=='I':
            fmt='1    ('+str(c)+'I'+str(ln)
        else :
            fmt='0    ('+str(c)+'G12.4' #+str(ln)            
        f1.write('INTERNAL     '+fmt+')     3  \n' )      
        if typ=='I':
            fmt='%'+str(ln)+'i'
        else :
            fmt='%+11.4e ' #'+str(ln)+'e '            
        s=''
        for i in range(l-1,-1,-1): # to write the rows from top to bottom
            for j in range(c):
                s+=fmt %(m[i][j])
            s+='\n'
        f1.write(s[:-1]) 

    def writeBlockModflow(self,m,f1,ktyp):
        #print shape(m),m
        if len(shape(m))==3:
            nlay,a,b=shape(m);
            for l in range(nlay):
                self.writeMatModflow(m[l],f1,ktyp)
                if l<nlay-1: f1.write('\n')
        elif self.core.addin.mesh != None: # unstructured case write nlay vectors
            nlay,a=shape(m);
            for l in range(nlay):
                self.writeVecModflow(m[l],f1,ktyp)
                if l<nlay-1: f1.write('\n')        
        else : self.writeMatModflow(m,f1,ktyp)
        
    #------------------------- fonction  write HBF -------------------
    def writeHfb(self,line,iz):       # EV 26/11/2018  
        ilay,irow,icol = self.xyzone2Mflow(self.core,line,iz)
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
        if core.mfUnstruct:
            nlay,ncell = getNlayers(core),core.addin.mfU.getNumber('elements') # only 1 layer up to now
            hd=zeros((nlay,ncell));#print('mfw 491', shape(hd))
        else :
            nlay,ncol,nrow = self.getGeom(core)
            ncell = ncol*nrow
            hd=zeros((nlay,nrow,ncol))
        try : f1 = open(self.fDir+os.sep+self.fName+'.head','rb')
        except IOError: return None
        blok=44+ncell*4; # v210 60
        for il in range(nlay):
            f1.seek(iper*nlay*blok+blok*il+44) #vpmwin
            data = arr2('f')
            data.fromfile(f1,ncell)
            if core.mfUnstruct: 
                hd[il] = data
            else : 
                m = reshape(data,(nrow,ncol)) #
                hd[il] = m[::-1] #=1::=1
        f1.close()  
        #modify the head if free and in 3D
#        if core.addin.getDim() in ['3D','Xsection','Radial']:   
#            if core.addin.getModelType()=='free':
#                hd = self.getHeadFree(hd)
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
        
    def getHeadFree(self,head):
        """the head in the layer where it is 0 must be replaced by the one
        of the layer below"""
        hd1 = head[1:]*1;#print 'mfw 414 hd',head[:-1]==0
        hd2 = hd1*(head[:-1]==0)
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
        if type(iper)==type([5]): iper=iper[0] # takes the thick only for the 1st tstep
        zb = core.Zblock
        thk = zb[:-1,:,:]-zb[1:,:,:]
        dim = core.addin.getDim()
        if dim in ['Xsection','Radial']:
            grd = core.addin.getFullGrid()
            nx,ny,dy = grd['nx'],grd['ny'],grd['dy'];
            ep =float(core.getValueFromName('Modflow','TOP'))-float(core.getValueFromName('Modflow','BOTM'))
            thk=ones((ny,1,nx))*ep
            return thk
        if core.addin.getModelType()=='free': 
            hd=self.readHeadFile(core,iper);#print hd
            if dim =='3D': hd=hd[0]
            thk[0,:,:]=hd-zb[1,:,:]
        return thk
    
    def getThicknessZone(self,core,iper,layers,ix,iy):
        dim = core.addin.getDim()
        thm = self.getThickness(core,iper);#print 'mflread 474',shape(thm),thm # only 2D up to now
        th =[]
        if dim in ['Xsection','Radial']:
            for i in range(len(ix)):
                th.append(thm[iy[i],0,ix[i]]) # revert for different orientation of layers
        else :
            for i in range(len(ix)):
                th.append(thm[layers[i],iy[i],ix[i]])
        return th

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
                print(ic)
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

        
