# -*- coding: utf-8 -*-
"""
Created on Wed Apr 15 18:21:30 2015

@author: olive
"""
import os
from .geometry import *

class opgeoWriter:

    def __init__(self,core,fDir,fName):
        self.core = core
        self.fDir,self.fName = fDir,fName
        self.fullPath = fDir+os.sep+fName;#print self.fullPath

    def writeFiles(self,core,opt):
        '''opt : flow or trans'''
        self.core = core
        self.opt = opt.lower()
        self.opgeo = core.addin.opgeo
        self.meshtype = core.getValueFromName('OpgeoFlow','O_GRID')
        self.nbCurves = 1
        self.curvefile = open(self.fullPath+'.rfd','w')
        #0 rectangular, 1 triangular
        lsuff=['gli','bc','ic','mfp','msh','mmp','msp','st','num','out','pcs','tim']
        if self.opt == 'flow':
            self.pcs = ['GROUNDWATER_FLOW']
            self.var = ['HEAD']
        elif self.opt == 'trans':
            self.pcs = ['GROUNDWATER_FLOW','MASS_TRANSPORT']
            self.var = ['HEAD','ConsTracer']
            lsuff.append('mcp')
        for suff in lsuff:
            fname = self.fullPath+'.'+suff
            f1 = open(fname,'w')
            exec('s = self.write'+suff.title()+'()')
            f1.write(s)
            f1.close()
            if s=='': os.system('rm '+fname)
        self.curvefile.write('#STOP')
        self.curvefile.close()
        
    def writeGli(self):
        #p_list = []
        i,s1,s2 = 0,'#POINTS \n',''
        llist = ['flow.'+str(a) for a in range(1,11)]
        if self.opt == 'trans': 
            llist.extend(['trans.1','trans.2'])
        ptnames = []
        for line in llist:
            if line[:4]=='flow': mod = 'OpgeoFlow'
            elif line[:5]=='trans' : mod = 'OpgeoTrans'
            nbz = 0
            if line in self.core.diczone[mod].dic:
                dicz = self.core.diczone[mod].dic[line]
                nbz = len(dicz['name'])
            for iz in range(nbz):
                a = ''
                if len(dicz['coords'][iz])==1:
                    coo = dicz['coords'][iz][0]
                    name = dicz['name'][iz]
                    if name in ptnames: continue
                    s1 += str(i)+' '+str(coo[0])[:7]+'  '+str(coo[1])[:7]+' 0.  0. 0. 0. 0. '
                    s1 +='$NAME '+name+'\n' # points are named
                    i+=1
                    ptnames.append(name)
                else :
                    for coo in dicz['coords'][iz]:
                        s1 += str(i)+' '+str(coo[0])[:7]+'  '+str(coo[1])[:7]+' 0.  0. 0. 0. 0. \n'
                        a += ' '+str(i) # makes the list of points for the polyline
                        i += 1
                    s2 += '#POLYLINE \n $NAME \n   '+dicz['name'][iz]+' \n $TYPE \n -1 \n'
                    s2 += ' $EPSILON \n 0.5 \n $POINTS \n'
                    s2 += a + '\n'
        return s1 + s2 + '#STOP \n'
      
    def writeBc(self):
        '''write the bc values'''
        s = ''
        for nb in ['flow.2']: #,'flow.3','flow.4'
            if nb not in list(self.core.diczone['OpgeoFlow'].dic.keys()): continue
            dicz = self.core.diczone['OpgeoFlow'].dic[nb]
            nbz,s = len(dicz['name']),''
            for iz in range(nbz):
                s += '#BOUNDARY_CONDITION \n $PCS_TYPE \n   GROUNDWATER_FLOW \n'
                s += ' $PRIMARY_VARIABLE \n   HEAD \n $GEO_TYPE \n  POLYLINE '
                s += dicz['name'][iz]+'\n'
                s += ' $DIS_TYPE \n'
                val = dicz['value'][iz]
                s += self.transient(val)
            if self.opt == 'trans':    
                dicz = self.core.diczone['OpgeoTrans'].dic['trans.2']
                nbz = len(dicz['name'])
                for iz in range(nbz):
                    s += '#BOUNDARY_CONDITION \n $PCS_TYPE \n   MASS_TRANSPORT \n'
                    s += ' $PRIMARY_VARIABLE \n ConsTracer \n  $GEO_TYPE \n  '
                    if len(dicz['coords'][iz])==1: ptyp = 'POINT '
                    else : ptyp = 'POLYLINE '
                    s += ptyp + dicz['name'][iz] +'\n'
                    s += ' $DIS_TYPE \n'
                    val = dicz['value'][iz]
                    s += self.transient(val)
        return s+'#STOP \n'

    def transient(self,val):
        s=''        
        if len(val.split())==1: 
            s += 'CONSTANT '+str(val)+'\n'
        else :
            s += 'CONSTANT 1\n'
            s += ' $TIM_TYPE\n  CURVE  '+str(self.nbCurves)+'\n'
            self.curvefile.write('#CURVES \n '+val+' \n')
            self.nbCurves += 1
        return s
        
    def writeIc(self):
        '''write the initial conditions'''
        s = ''
        kwl = [('GROUNDWATER_FLOW','HEAD',self.core.dicval['OpgeoFlow']['flow.1'])]
        if self.opt == 'trans':    
            kwl.append(('MASS_TRANSPORT','ConsTracer',self.core.dicval['OpgeoTrans']['trans.1']))
        for k in kwl:
            s += '#INITIAL_CONDITION \n $PCS_TYPE \n '+ k[0] +'\n'
            s += '$PRIMARY_VARIABLE \n '+k[1] + '\n'
            s += '$GEO_TYPE \n  DOMAIN \n $DIS_TYPE \n '
            s += 'CONSTANT '+str(k[2][0])+' \n'
        return s+'#STOP \n'
        
    def writeMcp(self): 
        s='#COMPONENT_PROPERTIES \n $NAME \n  ConsTracer \n $MOBILE \n   1\n'
        s+=' $DIFFUSION \n 1 '+ str(self.core.dicval['OpgeoTrans']['trans.8'][0])+'\n'
        return s+'#STOP \n'
        
    def writeMfp(self):
        '''fluid properties'''
        s = '#FLUID_PROPERTIES \n $FLUID_TYPE \n  LIQUID \n'
        s += ' $PCS_TYPE \n  HEAD \n'
        if self.opt=='trans': s += 'consTracer \n'
        s += ' $DENSITY \n  1 1.000000e+003 \n'
        s += ' $VISCOSITY \n  1 1.000000e-003 \n#STOP'
        return s
        
    def writeMmp(self):
        s = '#MEDIUM_PROPERTIES  \n $GEOMETRY_DIMENSION \n  2 \n'
        s += ' $GEOMETRY_AREA \n  1.000000e+000 \n $GEO_TYPE \n  DOMAIN \n'
        s += ' $POROSITY \n  1 '+ str(self.core.dicval['OpgeoFlow']['flow.8'][0]) +'\n'
        s += ' $TORTUOSITY \n  1   1.000000e+000 \n'
        s += ' $STORAGE  \n  1 '+ str(self.core.dicval['OpgeoFlow']['flow.7'][0]) +'\n'
        if self.core.dictype['OpgeoFlow']['flow.5'][0] in ['formula','interpolate']:
            s += ' $PERMEABILITY_DISTRIBUTION \n     permeabilities.txt \n'  
            self.writeExtFile('OpgeoFlow','flow.5','permeabilities.txt')
        else :
            s += ' $PERMEABILITY_TENSOR \n'
            s += '  ISOTROPIC '+ str(self.core.dicval['OpgeoFlow']['flow.5'][0]) +'\n'
        s += ' $MASS_DISPERSION \n  1 '
        s += str(self.core.dicval['OpgeoTrans']['trans.5'][0])+' '
        s += str(self.core.dicval['OpgeoTrans']['trans.6'][0])+' \n'
        return s+'#STOP \n'
        
    def writeExtFile(self,modName,variable, fname):
        """writes an external file for the given variable"""
        s = '#MEDIUM_PROPERTIES_DISTRIBUTED\n $MSH_TYPE\n  GROUNDWATER_FLOW \n'
        s += '$MMP_TYPE \n  PERMEABILITY \n'
        s += '$DIS_TYPE\n  NEAREST_VALUE ; or GEOMETRIC_MEAN\n'
        s += '$CONVERSION_FACTOR \n  1.0 \n $DATA \n'   
        # creates the x,y,z,k  matrix
        V = self.core.getValueLong(modName,variable,0)
        b = zeros((self.opgeo.nel,4))
        xc, yc = self.opgeo.elcenters
        b[:,0],b[:,1],b[:,3] = xc,yc,V
        s += self.opgeo.arr2string(b)+'\n#STOP\n'
        f1 = open(self.fDir+os.sep+fname,'w');f1.write(s);f1.close()
        
    def writeSt(self):
        """write teh source term, up to now only flow wells"""
        s = ''
        if 'flow.9' in self.core.diczone['OpgeoFlow'].dic:
            dicz = self.core.diczone['OpgeoFlow'].dic['flow.9']
            nbwells = len(dicz['value'])
            for iw in range(nbwells):
                s += '#SOURCE_TERM\n $PCS_TYPE\n  GROUNDWATER_FLOW\n'
                s += '$PRIMARY_VARIABLE\n  HEAD\n $GEO_TYPE\n  POINT '
                s += dicz['name'][iw] +'\n'
                s += ' $DIS_TYPE\n   '
                val = dicz['value'][iw]
                s += self.transient(val)
            s +'#STOP'
#                
#        if self.opt == 'trans' and self.core.diczone['OpgeoTrans'].dic.has_key('trans.2'):    
#            dicz = self.core.diczone['OpgeoTrans'].dic['trans.2']
#            nbwells = len(dicz['value'])
#            for iw in range(nbwells):
#                s += '#SOURCE_TERM\n $PCS_TYPE\n  MASS_TRANSPORT\n'
#                s += '$PRIMARY_VARIABLE\n  ConsTracer\n $GEO_TYPE\n  '
#                if len(dicz['coords'][iw])==1: ptyp = 'POINT '
#                else : ptyp = 'POLYLINE '
#                s += ptyp+dicz['name'][iw] +'\n'
#                s += ' $DIS_TYPE\n  CONSTANT  '
#                s += str(dicz['value'][iw])+'\n'            
            
        return s
        
    def writeMsh(self):
        coords = self.core.addin.getFullGrid()
        s = '#FEM_MSH \n $PCS_TYPE \n  LIQUID_FLOW \n $NODES \n'
        if self.meshtype==0: # rectangular mesh
            s += self.writeNodesCoords(coords)
            s += '$ELEMENTS \n'
            s += self.writeIncidence(coords)
        elif self.meshtype==1: # triangular mesh
            s += str(self.opgeo.nnod)+'\n'+self.opgeo.nodestring+'\n'
            s += '$ELEMENTS \n'
            s += str(self.opgeo.nel)+'\n'+self.opgeo.elementstring+'\n'            
        s += '$LAYER \n  1 \n #STOP'
        return s
        
    def writeMsp(self):
        s = '#SOLID_PROPERTIES \n   $DENSITY \n   1 2500 \n'
        s += ' $THERMAL \n  EXPANSION \n   1.0e-5 \n'
        s += '  CAPACITY \n  1 6000.\n  CONDUCTIVITY\n   1 5 \n'
        return s+'#STOP \n'
        
    def writeNum(self):
        s = ''
        for p in self.pcs:
            s += '#NUMERICS \n  $PCS_TYPE \n  '+p+'\n $LINEAR_SOLVER \n'
            s += '; method error_tolerance max_iterations theta precond storage \n'
            s += '2      1 1.e-10       1000           1.0   100     4 \n'
            s += '$ELE_GAUSS_POINTS \n   2 \n'
            if self.opt=='trans':
                s += ' $NON_LINEAR_SOLVER \n'
                s += '; method error_tolerance max_iterations relaxation \n'
                s += 'PICARD 1e-3            25             0.0 \n'
        return s+'#STOP \n'  

    def writeOut(self):
        s = ''
        s += '#OUTPUT \n $NOD_VALUES \n'
        for v in self.var: s+= v+'\n'
        s += ' $GEO_TYPE\n  DOMAIN \n $DAT_TYPE \n  TECPLOT \n $TIM_TYPE \n'
        t = self.core.dicaddin['Time'];#print t
        tend,tstep = float(t['final'][0]),float(t['steps'][0])
        tlist = arange(0,tend+tstep/10,tstep)
        for t in tlist: s+= str(t) + '\n'
        return s+'#STOP \n'  
        
    def writePcs(self):
        s = ''
        for p in self.pcs:
            s+='#PROCESS \n $PCS_TYPE \n '+p+ '\n'
            s += '$NUM_TYP \n  NEW \n'
        return s+'#STOP \n'
        
    def writeTim(self):
        itunit = self.core.dicval['OpgeoFlow']['domn.7'][0]
        tname = self.core.dickword['OpgeoFlow'].lines['domn.7']['detail'][0][itunit]
        s = '#TIME_STEPPING \n $PCS_TYPE \n'
        for p in self.pcs: s+= '  '+p + '\n'
        s += ' $TIME_UNIT \n  '+tname + '\n'
        t = self.core.dicaddin['Time']
        s += ' $TIME_START \n  0.0 \n $TIME_END \n  '+t['final'][0]+ '\n'
        #s += ' $TIME_STEPS \n  500  8640 \n
        s += ' $TIME_CONTROL \n  SELF_ADAPTIVE \n 1000 1.1 \n 1000 0.9 \n'
        s += '  MAX_TIME_STEP \n  '+t['steps'][0]+'\n'
        minstep = float(t['steps'][0])/1000
        s += '  MIN_TIME_STEP \n  '+str(minstep)+'\n'
        s += '  INITIAL_STEP_SIZE \n  '+str(minstep)+'\n'
        return s+'#STOP \n'
            
    def writeNodesCoords(self,coords):
        '''writes nodes coords from a square grid, x0,y0,z0 coords of
        the starting point, dx,dy,dz grid cells dimensions
        up to now in 2D'''
        s = ''
        x0,y0,z0 = coords['x0'],coords['y0'],0
        dx,dy,dz = coords['dx'],coords['dy'],[0]
        xvect = concatenate(([x0],x0+cumsum(dx)))
        if len(dy)==1: yvect = [y0]
        else : yvect = concatenate(([y0],y0+cumsum(dy)))
        if len(dz)==1: zvect = [z0]
        else : zvect = concatenate(([z0],z0+cumsum(dz)))
        i = 0
        s += str(len(xvect)*len(yvect)*len(zvect))+'\n'
        for z in zvect:
            for y in yvect:
                for x in xvect:
                    s += str(i)+' '+str(x)+' '+str(y)+' '+str(z)+'\n'
                    i += 1
        return s
    
    def writeIncidence(self,gcoords):
        # for a simple 2D square grid
        nx, nz, nl = len(gcoords['dx']),len(gcoords['dy']),1
        nel = nx*nz
        iel,s = 0,str(nel)+'\n'
        media = zone2mesh(self.core,'OpgeoFlow','flow.5',opt='zon') # only uses permeability zones
        #print 'ogw l220',shape(media),nz,nx
        for iz in range(nz):
            for ix in range(nx):
                n0 = iz*(nx+1)+ix
                nn = [n0,n0+1,n0+nx+2,n0+nx+1]
                s += str(iel)+' '+str(media[iel])+' quad '
                s += ' '.join([str(x) for x in nn])
                s+= '\n'
                iel += 1
        #print 'incid',time.time()-self.t0
        return s
    
'''////////////////////////////////////////////////////////////////
/////////////////////////////// READER     ///////////////////////
////////////////////////////////////////////////////////////////
'''    

class opgeoReader:
    def __init__(self,fDir,fName,option):
        #self.core = core
        self.fDir,self.fName = fDir,fName
        self.fullPath = self.fDir + os.sep+ self.fName
        self.read = False
        self.option = option
        
    def readNodeFile(self,core,ivar,istep):
        """reads a file with two variable with size of grid known"""
        meshtype = core.getValueFromName('OpgeoFlow','O_GRID')
        nnod = core.addin.opgeo.nnod
        if self.option in ['flow','heat']: nvar = 1
        elif self.option =='trans': nvar =2
        if meshtype==0: mtype='quad' #rect
        else : mtype ='tri'
        if self.read == False:
            f=open(self.fullPath+'_domain_'+mtype+'.tec','r')
            a=f.read()
            f.close()
            b=a.split('VARIABLE')[1:] #separates each time step
            ldata = []
            for s in b:
                s1 = s.split('\n1  ')[0] #removes the cell indices
                s2 = s1.split('\n')[3:nnod+3] # removes the first three lines and separates lines
                for iv in range(nvar):
                    l1 = [float(x.split()[3+iv]) for x in s2]
                    data = array(l1,ndmin=2)
                    ldata.append(data)
            nlay = 1
            nt=int(len(ldata)/nvar)
            self.data = zeros((nvar,nt,nlay,nnod));#print nl,nc,nt,shape(data)
            for it in range(nt):
                for iv in range(nvar):
                    self.data[iv,it,0] = ldata[it*nvar+iv]
            self.read = True
        return self.data[ivar,istep]

    def readHeadFile(self,core,tstep):
        hd = self.readNodeFile(core,0,tstep) # pressure is the first column
        #hd = (hd[:,:,1:]+hd[:,:,:-1])/2
        #hd = (hd[:,1:,:]+hd[:,:-1,:])/2
        return hd
        
    def readUCN(self,core,opt,tstep,nb,name):
        cnc = self.readNodeFile(core,1,tstep) # cnc or temp is the 1st var
        #cnc = (cnc[:,:,1:]+cnc[:,:,:-1])/2
        #cnc = (cnc[:,1:,:]+cnc[:,:-1,:])/2
        return cnc
