#from config import *
import os
import xml.dom.minidom as xdom
from geometry import *
#from wxDialogs import *
#from qtDialogs import *
from functools import partial

class impFile:
    """a generic class to import and export different types of files
    each method can be called with a file"""
    def __init__(self,gui,core): #,modName):
        self.core,self.gui = core,gui
        #self.dickword = core.dickword[modName];#print self.dickword
        
    def impGrid(self,fileDir=None,fileName=None):
        """import a grid 2 or 3D to inut a variable"""
        if fileDir==None : fileDir,fileName=self.fileDialog()
        pass
        
    def impZones(self,fileDir,fileName,modName,line):
        """ import a tabulated file with cols : 'name',name of zone,'value',
        value attached to the zone,'coord', and as many colums as coordinates
        organissed x1 y1 x2 y2...
        then add the zoens to the aquifer and visu"""
        fullName = fileDir+os.sep+fileName+'.txt'
        f1=open(fullName,'r');
        name,val,med,coords=[],[],[],[]
        for ll in f1:
            n=ll[ll.find('Name')+len('Name'):ll.rfind('Value')]
            v=ll[ll.find('Value')+len('Value'):ll.rfind('Media')]
            m=ll[ll.find('Media')+len('Media'):ll.rfind('Coord')]
            c=ll[ll.find('Coord')+len('Coord'):len(ll)]
            name.append(n.strip()) 
            val.append(v.strip().replace('\\n','\n')) # OA 23/9/20 added \\n
            if not ',' in m :
                med.append(int(m.strip())) #EV 19/02/20
            else : 
                m=m.strip()
                m=m.split(',')
                m = [int(i) for i in m]
                med.append(m)
            c0=c.strip();a=[] 
            c0=c0.split()
            for j in range(0,len(c0),2):
                a.append((float(c0[j]),float(c0[j+1])))
            coords.append(a)
            #l1=ll.split('\t');
            #name.append(l1[1]) 
            #val.append(l1[3].replace('\\n','\n'))
            #if len(l1[5])==1 : #EV 08/20
            #if not ',' in l1[5] :
                #med.append(int(l1[5])) #EV 19/02/20
            #else : 
                #m=l1[5].split(',')
                #m = [int(i) for i in m]
                #med.append(m)
            #c0=l1[7:];a=[] 
            #for j in range(0,len(c0),2):
                #a.append((float(c0[j]),float(c0[j+1])))
            #coords.append(a)
        for i in range(len(name)):
            dicz = self.core.diczone[modName]
            dicz.addZone(line)
            iz = dicz.getNbZones(line)-1
            dicz.dic[line]['name'][iz] = name[i]
            dicz.dic[line]['value'][iz] = val[i]
            dicz.dic[line]['media'][iz] = med[i]
            dicz.dic[line]['coords'][iz] = coords[i]
            self.gui.visu.addZone(med[i],name[i], val[i], coords[i],False)
        self.gui.visu.redraw()

    def impZones_old(self,fileDir,fileName,modName,line):
        """ import a tabulated file with cols : 'name',name of zone,'value',
        value attached to the zone,'coord', and as many colums as coordinates
        organissed x1 y1 x2 y2...
        then add the zoens to the aquifer and visu"""
        fullName = fileDir+os.sep+fileName+'.txt'
        f1=open(fullName,'r');
        name,val,med,coords=[],[],[],[]
        for ll in f1:
            l1=ll.split('\t');
            name.append(l1[1])
            val.append(l1[3])
            med.append(l1[5])
            c0=l1[7:];a=[]
            for j in range(0,len(c0),2):
                a.append((float(c0[j]),float(c0[j+1])))
            coords.append(a)
        for i in range(len(name)):
            dicz = self.core.diczone[modName]
            dicz.addZone(line)
            iz = dicz.getNbZones(line)-1
            dicz.dic[line]['name'][iz] = name[i]
            dicz.dic[line]['value'][iz] = val[i]
            dicz.dic[line]['media'][iz] = med[i]
            dicz.dic[line]['coords'][iz] = coords[i]
            self.gui.visu.addZone(med[i],name[i], val[i], coords[i])
        self.gui.visu.redraw()

    def impTabFile(self,fullName,titleNbCol=1,titleNbRow=1): # EV 03/20
        """ opens a file with 1st line titles and each line with a name in 1st
        column and send it backs as a dict with rows, cols, data
        titleNbCol indicates the number of columns that don't need to be read
        in each line, normally 1"""
        f1=open(fullName,'r'); #utils//
        if titleNbRow >1:
            tiCol=[]
            for i in range(titleNbRow):
                tiCol.append(f1.readline().split())
        else :tiCol=f1.readline().split()
        if titleNbCol>0: tiCol=tiCol[titleNbCol:]#;print ('ticol',tiCol)
        nC=len(tiCol)
        dat0=[];tiL=[]
        for ll in f1:
            l1=ll.split();#l1[-1]=l1[-1][:-1]
            #print('l1',l1)
            titl = ''
            if titleNbCol>0: titl = l1[0]
            tiL.append(titl);dat0.append(l1[titleNbCol:])
        nl=len(tiL) #; print('nl',nl)
        data=zeros((nl,nC))*0.;#print 'in model',dat0
        for l in range(nl):
            for c in range(nC):
                data[l,c]=float(dat0[l][c])
        return {'cols':tiCol,'rows':tiL,'data':data}
        
    def impVersion1(self,fileDir=None,fileName=None):
        zonkwd = {'Permeabilite':('Modflow','lpf.8'),'Potentiel':('Modflow','bas.5'),'Recharge':('Modflow','rch.2'),
            'Mur':('Modflow','dis.7'),'Toit':('Modflow','dis.6'),'Forages':('Modflow','wel.1'),
            'Porosite':('Mt3dms','btn.11'),'Transport':('Mt3dms','btn.13'),'PHT3D':('Pht3d','ph.4'),
            'PH_RECH':('Pht3d','ph.5')}
        filename = fileDir+os.sep+fileName + '.ipht'
        f1 = file(filename, 'r');doc = f1.read();f1.close()
        dom=xdom.parseString(doc)        
        dicts=dom.getElementsByTagName("dict")
        #basezone = {'number':[],'name':[],'coords':[],'media':[],'value':[],'type':[]}
        for d in dicts:
            dname = d.getElementsByTagName("name")[0].childNodes[0].data
            if dname != 'Zones': continue
            keys = d.getElementsByTagName("key")
            for k in keys: # groups of zones
                kname = k.getElementsByTagName("name")[0].childNodes[0].data
                if kname not in list(zonkwd.keys()): continue
                modName,line = zonkwd[kname]
                k0 = k.getElementsByTagName("content")[0].childNodes[0].data
                zgroup = eval(k0)
                if len(zgroup)==0 : continue
                # kdata is a dict of zones 
                diczone = self.core.diczone[modName]
                self.core.dictype[modName][line]=['zone']
                for i,z in enumerate(zgroup):
                    diczone.addZone(line)
                    diczone.dic[line]['number'][i]=i
                    diczone.dic[line]['media'][i]=0
                    diczone.dic[line]['name'][i]=z['nom']
                    diczone.dic[line]['coords'][i]=z['xy']
                    val = z['val']
                    if type(val) in [type([5]),type([5.])]: val = '\n'.join(str(val))
                    diczone.dic[line]['value'][i]=val
    
    def impAsciiGrid(self,fileDir,fileName): #EV 04/02/19
        f1 = open(fileDir+os.sep+fileName,'r')
        s = f1.read();f1.close()
        s1 = s.split('\n')
        # first lines
        i,dct=0,{}
        while len(s1[i].split())==2:
            a1,a2 = s1[i].split()
            dct[a1.strip().lower()] = a2.strip()
            i+=1
        l0=[]
        for j in range(int(dct['nrows'])): 
            l0.append(s1[i+j].split())
        return array(l0)
    
    def impGridVar(self,fileDir,fileName): #EV 02/04/20 #EV 28/04/20 #OA 17/10/20
        '''this imports a var file, 1st line 1 for a matrix oriented along increasing y
        and -1 for a modflow oriented grid (inverse to y), then spacing in x direction
        spacing in y direction and then the matrix of values
        if there are several variables, the number of variables is wrtten after ysign'''
        f1 = open(fileDir+os.sep+fileName,'r')
        s = f1.readlines()
        cols = [float(val) for val in s[1].split()]
        rows = [float(val) for val in s[2].split()]
        nc,nr = len(cols),len(rows)
        if len(s[0].split())==1: # OA modified 17/10/20
            ysign = int(s[0]);nvar = 1
        else :
            ysign, nvar = [int(a) for a in s[0].split()]
        arr = zeros((nvar,nr,nc))
        for iv in range(nvar): # OA modified 17/10/20
            for ir in range(nr):
                arr[iv,ir,:] = [float(x) for x in s[iv*nr+ir+3].split()]
        if nvar==1: arr=arr[0]
        return ysign,cols,rows,arr # OA added 13/6/20

    def imp3DgeomDis(self,core,fileName):
        '''
        this imports geometry as gvar from a disrw file produced by gWvistas 
        (will be modified to be more general later)
        this works only for a regular grid
        '''  
        f1=open(fileName+'.disrw','r')
        for i in range(4): f1.readline()
        a=f1.readline();b=a.split()
        nlay,nrow,ncol=int(b[0]),int(b[1]),int(b[2]);print (nlay,nrow,ncol)
        f1.readline();
        #reading col width
        a=f1.readline();b=a.split()
        nc=0;wx=[]
        if int(b[0])==0: 
            wx=ones(ncol)*float(b[1].split('(')[0])
        else:
            while nc<ncol:
                a=f1.readline().split();wx.extend(a);nc+=len(a)
            wx=array(wx).astype('float')
        #reading row height
        nr=0;hy=[]
        f1.readline();b=a.split()
        if int(b[0])==0: 
            hy=ones(nrow)*float(b[1].split('(')[0])
        else:
            while nr<nrow:
                a=f1.readline().split();hy.extend(a);nr+=len(a)
            hy=array(hy).astype('float')
        #reading layer top
        ztop=zeros((nlay,nrow,ncol));flgTop=zeros(nlay)
        for ilay in range(nlay):
            a=f1.readline()
            if int(a.split()[0])==11:
                flgTop[ilay]=1
                for irow in range(nrow):
                    nc=0;z=[]
                    while nc<ncol:
                        a=f1.readline().split();z.extend(a);nc+=len(a)
                    ztop[ilay,irow]=array(z).astype('float')
        zbot=zeros((nrow,ncol));flgBot=0
        a=f1.readline()
        if int(a.split()[0])==11:
            flgBot=1
            for irow in range(nrow):
                nc=0;z=[]
                while nc<ncol:
                    a=f1.readline().split();z.extend(a);nc+=len(a)
                zbot[irow]=array(z).astype('float')
        f1.close()
        # writing 
        for ilay in range(nlay):
            if flgTop[ilay]==1:
                f2=open(core.fileDir+os.sep+'top'+str(ilay)+'.gvar','w')
                f2.write('-1\n')
                f2.write(' '.join(wx.astype('str'))+'\n')
                f2.write(' '.join(hy.astype('str'))+'\n')
                for irow in range(nrow):
                    f2.write(' '.join(ztop[ilay,irow].astype('str'))+'\n')
                f2.close()
        if flgBot==1:
            f2=open(core.fileDir+os.sep+'botm.gvar','w')
            f2.write('-1\n')
            f2.write(' '.join(wx.astype('str'))+'\n')
            f2.write(' '.join(hy.astype('str'))+'\n')
            for irow in range(nrow):
                f2.write(' '.join(zbot[irow].astype('str'))+'\n')
            f2.close()
        return True
                  
''' below not used
'''
class impAsciiModflow:
    """this class imports a whole modflow model, and sets it as an
    iqpht3d model. It reads each file according to the keyword dictionnary
    and gets the value according to conditions and type of values
    dic1 is the dict that stores the vaues during reading
    After that the grid cannot be modified"""
    def __init__(self,core,fileDir=None,fileName=None):
        #if fileDir==None : fileDir,fileName=self.fileDialog()
        self.fileDir,self.fileBase = fileDir,fileName.split('.')[0]
        self.dickword = core.dickword['Modflow'];#print self.dickword
        self.lexceptions=['dis.8','1pf.8','1pf.9','1pf.10','1pf.11','rch.2']
        self.core = core
        
    def readVec(self,line,size,typ):
        """reads a vector, which may be on several lines
        """
        l0=self.s_in[self.indx].split();
        self.indx += 1
        if l0[0]=='CONSTANT': vec= float(l0[1])
        else :
            nread=0;l0=[]
            nval=size[0];#print size,nval
            while nread<nval:
                l1=self.s_in[self.indx].split();self.indx += 1
                #print 'readline',nread
                l0.extend(l1)
                nread+=len(l1)
            if typ[3:] in ['int','float']: vec = list(array(l0).astype(typ[3:]))
            else : vec= list(l0)
        self.core.dicval['Modflow'][line] = vec

    def readArray(self,line,size,typ):
        """reads an array with the modflow formulation and transforms it to
         an array in ipht3d
        """
        typ1=typ[3:]
        #nlay=self.core.getValueFromName('Modflow','NLAY')
        vals,arr = [],zeros((size))
        if len(size)==2: arr = self.getMat(size,typ1)
        else : 
            for il in range(size[0]): arr[il] = self.getMat(size,typ1);#print il
        #print 'impexp 144',line,size, typ,typ1,shape(arr),type(arr[0,0])
        self.core.dicarray['Modflow'][line] = arr;#print sum(arr)
        self.core.dictype['Modflow'][line] = ['array']
        
    def getMat(self,size,typ1):
        """reads a matrix in the ascii file, size is the global size of the arrays (with layers)
        typ1 is float or int
        """
        l0=self.s_in[self.indx];self.indx += 1;#print 'impexp 151',l0
        if l0[:10]=='CONSTANT  ' or ' 0' in l0[:10]: 
            v1=float(l0[11:20])
            mat = v1*ones(size[-2:])
        else : 
            nread=0;l0=[] # start the reading of several lines
            nval=size[-2]*size[-1];#print size,nval
            while nread<nval: # read several lines
                l1=self.s_in[self.indx].split();self.indx += 1;
                l0.extend(l1)
                nread+=len(l1)
            m = array(l0).astype(typ1)
            mat = reshape(m,size[-2:]);#print 'impexp 155',shape(mat) # this is the read matrix
        return mat[-1::-1,:] # python has different orientaiton than modflow
                    
    def readStrLay(self,size,typL):
        nval=size[0]
        nlines=nval/40+1
        l0=[]
        for n in range(nlines):
            lval = array(self.s_in[self.indx].split()).astype(typL)
            l0.extend(list(lval))
            self.indx += 1
        return l0
        
    def readTransient(self,nvar):
   #def readTransient(f1,nvar)
        """this allows to read a transient file: wel, drn, riv...
        returns a matrix (or list of matrices)
        """
        def lst2mat(shap,indx,vec):
            """fills an array with values read from indx: 1st layer, 2nd row, 3rd col
            and the vector of values for these positions
            """
            nlay,nrow,ncol = shap
            il,ir,ic = indx[:,0]-1,indx[:,1]-1,indx[:,2]-1 # indices in mdflow start at 1
            mat = zeros((nlay,nrow,ncol))
            ind2 = il*nrow*ncol+ir*ncol+ic
            put(mat,ind2,vec); #print shap,ind2,vec,mat
            return mat[:,-1::-1,:] #differnt orienation modflow python
        shap=self.core.dicval['Modflow']['dis.2'][:3];#print 'impexp l 191',shap
        b = self.s_in[self.indx:]
        b.remove('')
        nbmax = int(b[0][:10]); b= b[1:] # get max nb of points and remove line 
        nbper = int(b[0][:10]); b= b[1:] # get nb of points/period and remove line 
        b=b[:nbper] # takes only the points of the 1st perdio (up to now)
        if len(b)==0: return None # for some cases there are no data (stupid writer)
        b1=[[s[i*10:i*10+10] for i in range(3+nvar)] for s in b] # cut each line by 10 characters
        arr=array(b1,ndmin=2)
        lstarr = []
        indx = arr[:,:3].astype('int')
        for iv in range(nvar):
            vec = arr[:,3+iv].astype('float')
            lstarr.append(lst2mat(shap,indx,vec))
        #print 'impexp 200',nvar,len(lstarr)
        return lstarr
        
    def readExceptions(self,line,size,typ):
        #nper=self.core.dicval['Modflow']['dis.2'][3]
        if line =='dis.8': #read the time from periods
            vals = self.s_in[self.indx].split()
            self.core.dicval['Modflow'][line] = vals
            ntimes = len(self.s_in)-self.indx-1
            # we suppose all periods have the same length!!!!
            step = float(vals[0])
            d={'final': str(ntimes*step),'steps': str(step),'mode':'linear'}
            self.core.dicaddin['Time'] = d
            
        if line=='1pf.8':
            llines = ['lpf.8','lpf.9','lpf10','lpf11']
            for ll2 in llines:
                if self.testCondition(self.dickword.lines[ll2]['cond'])==False:
                    llines.remove(ll2)
            nbv = len(llines)
            siz2 = size.insert(0,nbv)
            m1 = zeros(siz2)
            for il in range(nlayers):
                for iv in range(nbv):
                    m1[iv][il] = self.getMat(size,typ1)
            for i,ll in enumerate(llines) :
                self.core.dicarray['Modflow'][ll] = m1[i]
        if line=='rch.2':
            # !! ONE PERIOD now, rech is a 2d array for several periods
            self.indx += 1 # first line : data read or same as previous
            self.readArray(line,size[1:],typ[0]) # here we juste read a 2d variable, one layer
            
    def testCondition(self,cond):
        bool = self.core.testCondition('Modflow',cond);#print cond, bool
        return bool
    
    def setTypeList(self,typlist):
        """reads a simple list of values on one line, not so simple as formats can be different"""
        def test(typlist,nbv,lval):
            lout = []
            for i in range(nbv):
                if typlist[i] =='int': lout.append(int(lval[i])) # pb for choice
                elif typlist[i]=='float': lout.append(float(lval[i]))
                else : lout.append(lval[i])
            return lout
        nbv=len(typlist)
        val=self.s_in[self.indx];self.indx += 1
        try : 
            lval = val.split()
            lout = test(typlist,nbv,lval)
        except ValueError : 
            lval = [val[i*10:i*10+10] for i in range(nbv)]
            lout = test(typlist,nbv,lval)
        return lout 
            
    def readFile(self,lines,filename):
        '''reads any type of modflow file and according to the type of value
        defined in keywords, chooses the type of reader and send the values
        to the core'''
        f1=open(filename,'r')
        self.s_in = f1.read().split('\n')
        f1.close()
        self.indx = 0
        while self.s_in[self.indx][0] in ["#",'P']: self.indx += 1 # reading the firs tlines
        core = self.core
        md = 'Modflow'
        tr_nvar = {'drn':2,'ghb':2,'riv':3,'wel':1,'chd':2} #nb of variable for transient files
        for ll in lines:
            typ=self.dickword.lines[ll]['type']
            kw = self.dickword.lines[ll]['kw'][0]
            cond = self.dickword.lines[ll]['cond']
            #varname = self.dickword.lines[ll]['comm'][:4]; #print ll,typ[0][:3],typ[0][:3]=='arr'
            size = self.core.getSizeKw(md,kw);#print ll,size,typ
            if ll in self.lexceptions :
                self.readExceptions(ll,size,typ)
                continue
            if self.testCondition(cond)==False :
                continue
            if ll in ['drn.1','riv.1','wel.1','ghb.1','chd.1']:
                lstarr=self.readTransient(tr_nvar[ll[:3]])
                if lstarr!= None: 
                    core.dicarray[md][ll] = lstarr
                    core.dictype[md][ll] = ['array']
                continue
            if typ[0][:3]=='vec':
                self.readVec(ll,size,typ[0])
            elif typ[0][:3]=='arr':
                self.readArray(ll,size,typ[0])
            elif typ[0][:3]=='lay':
                typL = typ[0][3:]
                core.dicval[md][ll]=self.readStrLay(size,typL)
            elif typ[0]=='title':
                pass #already read above
            else :
                a= self.setTypeList(typ);#a list of values
                core.dicval[md][ll] = a
            #self.changeStoredValues(ll)
        f1.close() 
        
    def readLine(self):
        s = self.s_in[self.indx];self.indx+=1
        return s
        
    # now read the whole stuff
    def readAll(self):
        md = 'Modflow'
        dic1={}
        grpList=self.dickword.grpList
        filePath=self.fileDir+os.sep+self.fileBase
        f1=open(filePath+'.nam')
        txtnam=''.join(f1.readlines())
        f1.close()
        for g in grpList: # sweep through the existing modeuls in iPht3d
            g1=g.swapcase()
            if txtnam.find(g)==-1 and txtnam.find(g1)==-1: continue # the module is not implemented in iPht3d
            #print g
            self.core.dicaddin['usedM_Modflow'][1][grpList.index(g)]=True # add true to the usedmodule of modflow
            filename = filePath+'.'+g[:3]
            lines = self.dickword.groups[g]
            dic1=self.readFile(lines,filename)
        dx,dy = self.core.dicval[md]['dis.4'],self.core.dicval[md]['dis.5']
        nlay,nrow,ncol=self.core.dicval['Modflow']['dis.2'][:3]
        try : float(dx);dx=[dx]*ncol;self.core.dicval[md]['dis.4']=dx
        except TypeError: pass
        try: float(dy);dy=[dy]*nrow;self.core.dicval[md]['dis.5']=dy
        except TypeError: pass
        grd = {'x0':0,'y0':0,'x1':sum(dx),'y1':sum(dy),'nx':len(dx),'ny':len(dy),'dx':'fixed','dy':'fixed'}
        self.core.dicaddin['Grid'] = grd
        if nlay>1:
            self.core.dicaddin['Model']['dimension']='3D'
            self.core.dicaddin['3D']['topMedia']=list(range(nlay,0,-1))
        print ('done')
        return dic1
