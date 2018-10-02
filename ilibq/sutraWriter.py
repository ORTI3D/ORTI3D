# -*- coding: utf-8 -*-
"""
Created on Sat Jan 03 21:53:13 2015

@author: olive
"""
from array import array as arr2
import os,time
import sutraKeywords as Skey
from geometry import *
import time

class sutraWriter:

    def __init__(self,core,fDir, fName):
        self.core = core
        self.fDir,self.fName = fDir,fName

    def writeSutraFiles(self):
        fn = self.fDir+os.sep+self.fName
        # test transient state
        self.transient,self.t0 = False,time.time()
        linelist = ['ssm.17','ssm.18','ssm.19','ssm.20']
        self.ttable = self.core.makeTtable()
        for l in linelist :
            if self.ttable['Transient'].has_key(l):
                if self.ttable['Transient'][l]: self.transient = True
        # write the fil file
        f1=open(self.fDir+os.sep+'sutra.fil','w')
        for i,n in enumerate(['INP','ICS','LST']):
            f1.write('\''+n+'\' '+str(50+i)+' \''+fn+'.'+n.lower()+'\' \n')
        f1.write('\'NOD\' 54 \''+fn+'.nod\' \n')
        if self.transient:
            f1.write('\'BCS\' 55 \''+fn+'.bcs\' \n')
        f1.close()
        grpInput=['GLOB','SOLV','OUT','MED']
        #exceptions
        lexcept = ['glob.6b','glob.7','med.14b']
        lexcept.extend(['med.15'+x for x in ['b','c','d','e','f','g','h','i']])
        # need to write the nb of BC in glob 3
        self.setDims()
        self.setNbFixedValues()
        f1=open(fn+'.inp','w')
        s = ''
        for grp in grpInput:
            llist=Skey.groups[grp];#print n1,name
            s += self.writeLines(llist,lexcept)
        grp = 'SSM'
        llist=Skey.groups[grp]
        for line in llist:
            cond=Skey.lines[line]['cond']
            if self.testCondition(cond)==False : continue
            s += '#'+ line +'\n'
            if line != 'ssm.22':
                s += self.writeBCs(line)[1] # gets only the text values, not the nb
            else : 
                s += self.writeIncidence()
        f1.write(s)
        f1.close()
        ################" writes ics
        f1=open(fn+'.ics','w')
        s = self.writeIcs(); print('Ics written')
        f1.write(s)
        f1.close()
        ################" writes bcs
        if self.transient:
            f1=open(fn+'.bcs','w')
            s = self.writeBcs(); print('Bcs written')
            f1.write(s)
            f1.close()
        
    def setDims(self):
        grd = self.core.addin.getFullGrid()
        self.nx,self.ny = int(grd['nx']),int(grd['ny'])
        self.TriDim = False
        if self.core.addin.getDim()=='3D': self.TriDim = True
        self.nlay = getNlayers(self.core)
        self.Pinput = self.core.getValueFromName('Sutra','PINPUT')

    def setNbFixedValues(self):
        # just get the number of fixed conditions (does not write)
        nsop = self.writeBCs('ssm.17')[0] 
        self.core.setValueFromName('Sutra','NSOP',nsop)
        nsou = self.writeBCs('ssm.18')[0] 
        self.core.setValueFromName('Sutra','NSOU',nsou)
        npbc = self.writeBCs('ssm.19')[0] 
        self.core.setValueFromName('Sutra','NPBC',npbc)
        nubc = self.writeBCs('ssm.20')[0] 
        self.core.setValueFromName('Sutra','NUBC',nubc)
        
    def testCondition(self,cond):
        """ test if the condition is satisfied"""
        return self.core.testCondition('Sutra',cond)        
        
    def writeLines(self,llist,lexcept):
        s = ''
        for line in llist:
            cond=Skey.lines[line]['cond'];#print 'sutw 35',line
            if self.testCondition(cond)==False : continue
            s += '#'+ line +'\n'
            if line in lexcept : 
                s += self.writeExceptions(line);continue
            kwlist = Skey.lines[line]['kw']
            ktyp = Skey.lines[line]['type']
            kdetail = Skey.lines[line]['detail']
            lval=self.core.dicval['Sutra'][line];#print 'mtw 77',self.core.dicval['Mt3dms'],lval,kwlist,ktyp
            if ktyp[0] in ['vecint','vecfloat','arrint','arrfloat']:
                s += str(lval[0])+'\n' # !!! just one value
            elif ktyp[0]=='title': # case of a title line
                s += str(ktyp[0])+'\n'
            else : # classical keywords
                s1 = ''
                for ik in range(len(kwlist)):
                    if ktyp[ik] == 'choice':
                        s1 += '\''+kdetail[ik][lval[ik]+1].upper()+'\' '
                    else : 
                        s1 += str(lval[ik])+' '
                if line in ['glob.2a','glob.2b']: #stupid sutra syntax
                    s1 = s1.replace('\' \'',' ')
                s += s1+'\n'
            if line == 'glob.6b': s += '-\n' # to end the time list
        return s
        
    def writeExceptions(self,line):
        if line == 'med.14b': 
            s = self.writeNodes();print ('nodes written')
        elif line == 'med.15b': 
            s = self.writeElements();print ('elts written')
        elif line == 'glob.6b':
            s = self.writeTime()
        else :
             s= ''
        return s
        
    def writeTime(self):
        s = ''
        kdetail = Skey.lines['glob.6b']['detail']
        for ik in range(3): # frist three keywords
            s += '\''+kdetail[ik][1].upper()+'\' '
        s += str(self.core.getValueFromName('Sutra','SCALT'))+' '
        tl0 = self.ttable['tlist']
        stepdiv = int(self.core.getValueFromName('Sutra','STEPDIV'))
        tlist,t0,s1 = [],float(tl0[1])/stepdiv,'';#print t0,tl0
        for i,tf in enumerate(tl0[1:]):
            r1 = arange(t0,tf,(tf-t0)/stepdiv)
            s1 += ' '.join([str(a) for a in r1])+'\n'
            t0 = tf
        s1 += str(tl0[-1])
        ntime = (len(tl0)-1)*stepdiv+1
        s += str(ntime)+' '+s1+'\n - \n'
        return s
        
    def writeNodes(self):
        # up to now only in 2D
        # nreg,x,y,thickness,poro
        poro = self.core.dicval['Sutra']['med.14b'][0]
        nx,ny,xvect,yvect = getXYvects(self.core)
        znode = self.getZnode()
        i,s = 1,''
        nl = 1
        if self.TriDim : nl = self.nlay+1
        #print 'sutw l 164',nx,ny,nl,len(xvect),len(yvect)
        for il in range(nl):
            for iy,y in enumerate(yvect):
                for ix,x in enumerate(xvect):
                    if self.TriDim :
                        s += str(i)+' 0 '+str(x)[:7]+' '+str(y)[:7]+' '+str(znode[il,iy,ix])[:7]+\
                        ' '+str(poro)+'\n'
                    else :
                        s += str(i)+' 0 '+str(x)[:7]+' '+str(y)[:7]+' 1. '+str(poro)+'\n'
                    i += 1
        #print 'nodes',time.time()-self.t0
        return s
        
    def getZnode(self): 
        zblock = makeZblock(self.core);#print 'suw l 178 zblock',zblock
        if self.core.addin.getDim() in ['Radial','Xsection']: 
            zn = zblock
            #zn = (zblock[1:,:,:]+zblock[:-1,:,:])/2
            zn = (zn[:,:,1:]+zn[:,:,:-1])/2
            #zn = concatenate([zn[:1,:,:],zn,zn[-1:,:,:]],axis=0)
            znode = concatenate([zn[:,:,:1],zn,zn[:,:,-1:]],axis=2)
            return znode
        else : 
            zn = (zblock[:,1:,:]+zblock[:,:-1,:])/2
            zn = (zn[:,:,1:]+zn[:,:,:-1])/2
            zn = concatenate([zn[:,:1,:],zn,zn[:,-1:,:]],axis=1)
            znode = concatenate([zn[:,:,:1],zn,zn[:,:,-1:]],axis=2)
            if self.TriDim : return znode
            else : return znode[0]*0+1
        
    def writeElements(self):
        # up to now only one value, no anisotropy
        # lreg, pmax,pmin,angl1,,almax,almin,atmax,atmin
        dval = self.core.dicval['Sutra']
        perm = dval['med.15c'][0]
        aL,aT = dval['med.15g'][0],dval['med.15h'][0]
        nl = 1
        if self.TriDim : nl = self.nlay
        nel = self.nx*self.ny*nl
        s = ''
        for i in range(nel):
            if self.TriDim :
                s += str(i+1)+' 0 '+str(perm)+' '+str(perm)+' '+str(perm)+\
                    ' 0. 0. 0. '+str(aL)+' '+str(aL)+' '+str(aL)+' '+\
                    str(aT)+' '+str(aT)+' '+str(aT)+'\n'
            else :
                s += str(i+1)+' 0 '+str(perm)+' '+str(perm)+' 0. '+\
                    str(aL)+' '+str(aL)+' '+str(aT)+' '+str(aT)+'\n'
        #print 'elts',time.time()-self.t0
        return s
        
    def writeBCs(self,line):
        if line not in self.core.diczone['Sutra'].dic.keys(): 
            return 0,''
        diczone = self.core.diczone['Sutra'].dic[line]
        nb,s = 0,''
        znode = ravel(self.getZnode());
        for iz in range(len(diczone['coords'])): 
            nodes = self.zone2nodes(diczone['coords'][iz],diczone['media'][iz])
            nb += len(nodes);#print len(nodes)
            #val = diczone['value'][iz].split() # here it is two values
            val = self.ttable[line][0][iz].split();#print 'suwri l231',val # 0 for the first time step
            for n in nodes:
                if self.Pinput and line=='ssm.19': 
                    press = (float(val[0])+znode[n])*9.81*1e3
                    s += str(n)+' '+str(press)[:8]+' '+val[1]+'\n'
                else :
                    s += str(n)+' '+val[0]+' '+val[1]+'\n'
        s += '0\n'
        return nb,s
        
    def zone2nodes(self,coords,medialist):
        #up to now just for one line with 2 points in x or y dir
        # finds the x,y,z of a zone
        if type(medialist)==type(1): medialist= [medialist] # sometime there is just one layer
        x,y = zip(*coords)
        z = x*1
        # returns indices of cells below a zone as arrays
        ix,iy,iz = zone2index(self.core,x,y,z)
        # finds nodes from cell indices
        nn0 = iy*(self.nx+1)+ix+1;#print 'suw l245',shape(nn0),nn0
        #nn0l = list(nn0)
        if self.TriDim :
            laylist = []
            for m in medialist: laylist.extend(media2layers(self.core,m))
            nn2d = (self.nx+1)*(self.ny+1)
            #creates an array contiang nodes for all layers
            for i,lay in enumerate(laylist):
                if i==0: nn = list(nn0+nn2d*lay)
                else : nn.extend(list(nn0+nn2d*lay))
        else :
            nn = list(nn0)
        return nn

    def writeIncidence(self):
        # for a simple 2D square grid
        iel,s = 0,'\'INCIDENCE\' \n'
        nx, ny = self.nx,self.ny
        nl = 1
        if self.TriDim : nl = self.nlay+1
        node_lay = (nx+1)*(ny+1)
        nel = nx*ny*nl
        for il in range(nl):
            for iy in range(ny):
                for ix in range(nx):
                    iel += 1
                    n0 = il*node_lay+iy*(nx+1)+ix+1
                    nn = [n0,n0+1,n0+nx+2,n0+nx+1]
                    s += str(iel)+' '
                    if self.TriDim :
                        s += ' '.join([str(x+node_lay) for x in nn])+' '
                    s += ' '.join([str(x) for x in nn])
                    s+= '\n'
        #print 'incid',time.time()-self.t0
        return s

    def writeIcs(self):
        nodes,values = [],[]
        press = self.mextend(self.core.getValueLong('Sutra','ics.2',0));
        #print shape(press),press
        U = self.mextend(self.core.getValueLong('Sutra','ics.3',0));
        znode = ravel(self.getZnode());#print 'sutwics',len(press),len(U),len(znode)
        if self.core.diczone['Sutra'].dic.has_key('ssm.19'):
            diczone = self.core.diczone['Sutra'].dic['ssm.19']
            for iz in range(len(diczone['coords'])): 
                nn = self.zone2nodes(diczone['coords'][iz],diczone['media'][iz])
                p,u = diczone['value'][iz].split()
                press[nn] = float(p)
                U[nn]= float(u)
        if self.Pinput:
            press = (press+znode)*9.81*1e3
        #print 'suw l232',nodes,values
        sP,sU,i = '','',0
        nl = 1
        if self.TriDim : nl =self.nlay+1
        for iz in range(nl): 
            for iy in range(self.ny+1):
                for ix in range(self.nx+1):
                    sP += str(press[i])+' '
                    sU += str(U[i])+' '
                    i += 1
                    if mod(i,20)==0: 
                        sP+= '\n';sU += '\n'
        s = '0\n\'NONUNIFORM\'\n'+sP+'\n\'NONUNIFORM\'\n'+sU+'\n'
        #print 'Ics',time.time()-self.t0
        return s
        
    def mextend(self,m):
        m = concatenate((m,m[:,:,-1:]),axis=2)
        dm = self.core.addin.getDim();#print shape(m),dm
        if dm not in ['Radial','Xsection']:
            m = concatenate((m,m[:,-1:,:]),axis=1)
        if self.TriDim or (dm in ['Radial','Xsection']): 
            return ravel(concatenate((m,m[-1:,:,:]),axis=0))
        else :
            return ravel(m)
        
    def writeBcs(self):
        """writes the Bcs file if there are variation of BC along time"""
        # get data from ttable
        cdict,zones,nbnodes = {},{},[]
        linelist = ['ssm.17','ssm.18','ssm.19','ssm.20']
        # get the list of zones and number of variable BC (they will remain the same for all time steps)
        for line in linelist:
            if self.ttable.has_key(line):
                cdict[line] = self.ttable[line];#print "mfw transient",line,clist
                zones[line] = self.core.diczone['Sutra'].dic[line]
                nz = len(zones[line]['coords'])
                nn = 0
                for iz in range(nz):
                    nn += len(self.zone2nodes(zones[line]['coords'][iz],zones[line]['media'][iz]))
                nbnodes.append(nn)
            else :
                nbnodes.append(0)
        sBCnbs = ' '.join([str(a) for a in nbnodes])+'\n'
        tlist = array(self.ttable['tlist'])
        znode = ravel(self.getZnode());#print 'suw l 232',nodes,len(znode)
        stepdiv = int(self.core.getValueFromName('Sutra','STEPDIV')) # stepdivide
        per = range(len(tlist)*stepdiv)
        s = '\'TIME_STEPS\' \n'
        # loop on time steps
        for ip in per:
        # writes for a given time step its name and the number of nodes for each BC type
            s += '\'ts'+str(ip)+'\'  '+sBCnbs
        # loop on BCs (ssm.17 to 20)
            for line in linelist:
                if line not in zones.keys(): continue
        # loop on zones
                for iz in range(len(zones[line]['coords'])): 
                    nodes = self.zone2nodes(zones[line]['coords'][iz],zones[line]['media'][iz])
                    val = cdict[line][ip/stepdiv,iz].split();#print 'suwrite l375',val
                    for n in nodes:
                        if self.Pinput and line=='ssm.19': 
                            press = (float(val[0])+znode[n])*9.81*1e3
                            s += str(n)+' '+str(press)[:8]+' '+val[1]+'\n'
                        else :
                            s += str(n)+' '+val[0]+' '+val[1]+'\n'
                s += '0 \n'
        return s
        
class sutraReader:
    
    def __init__(self,fDir, fName):
        self.fDir,self.fName = fDir,fName
        self.fullPath = fDir+os.sep+fName;#print self.fullPath
        self.read = False
        self.nodeData = None # to store data 5 dims : col,time,z,y,x
        self.dim = 2
        
    def readNodeFile(self,core,icol,tstep):
        #read the node file and transform it to an array (for square grids)
        dm = core.addin.getDim()
        if dm=='3D': self.dim = 3
        if self.read == False : # file not read before
            nx = core.getValueFromName('Sutra','NN1')
            ny = core.getValueFromName('Sutra','NN2')
            nlay = core.getValueFromName('Sutra','NN3')
            nn = nx*ny*nlay
            ntimes = int(core.getValueFromName('Sutra','NTLIST')); 
#            m = loadtxt(self.fullPath+'.nod',comments='#')
#            self.nodeData = zeros((nc-self.dim,ntimes,nlay,ny,nx))
#            for ic in range(nc-self.dim):
#                self.nodeData[ic] = reshape(m[:,ic+self.dim],(ntimes,nlay,ny,nx))
            lout = self.loadTxtFile(self.fullPath+'.nod')
            nt = len(lout)
            if dm in ['Radial','Xsection']:
                nlay =ny; ny = 1
            self.nodeData = zeros((2,nt,nlay,ny,nx))
            for it in range (nt):
                self.nodeData[0,it] = reshape(lout[it][0],(nlay,ny,nx))
                self.nodeData[1,it] = reshape(lout[it][1],(nlay,ny,nx))
            self.read = True
        return self.nodeData[icol,tstep]
        
    def loadTxtFile(self,fname):
        '''more complex but faster than loadtxt'''
        f1=open(fname,'r')
        a=f1.read()
        f1.close()
        b=a.split('## TIME STEP    ')[1:]
        lout = []
        a1,a2,a3,a4 = 33,47,48,62
        if self.dim==3: a1,a2,a3,a4 = 48,62,63,78
        for s in b[:-1]:
            s1=s.split('\n')[3:-3]
            l1=[float(a[a1:a2]) for a in s1]
            l2=[float(a[a3:a4]) for a in s1]
            lout.append((l1,l2))
        # for last time step
        s1=b[-1].split('\n')[3:]
        l1=[float(a[a1:a2]) for a in s1[:-1]]
        l2=[float(a[a3:a4]) for a in s1[:-1]]
        lout.append((l1,l2))
        return lout
    
    def readHeadFile(self,core,tstep):
        return self.readNodeFile(core,0,tstep) # pressure is the first column
    def readUCN(self,core,opt,tstep,iesp,spname):
        return self.readNodeFile(core,1,tstep) # conc is the 2nd column
        
    def getPtObs(self,core,irow,icol,ilay,iper,opt,iesp=0):
        """a function to values of one variable at given point or points.
        irow, icol and ilay must be lists of the same length. iper is also
        a list containing the periods for which data are needed."""
        # only 2D now
        #print irow,icol,ilay,iper,opt
        npts=len(irow);#print self.read,iper, shape(self.nodeData)
        pobs=zeros((len(iper),npts))+0.;
        if opt in ['Head','flux'] : nc = 0 # flux is worong just for test
        elif opt == 'Tracer': nc = 1
        for tstep in range(len(iper)):
            m = self.readNodeFile(core,nc,tstep)
            for i in range(npts):
                pobs[tstep,i]=m[0,irow[i],icol[i]]
        return pobs

    
