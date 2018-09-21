"""PEST utility for iPht3D :
steps :
1. ipht3d seek for the pest_parm.txt file
the parameters have to start with MF, MT or PH for modflow, mt3d or pht3d
This file must be in the same folder than the ipht3d file
for permeability : use MFKz1 where z1 is the name of the zone of interest
for dispersivity : MTPaL or MTPaT MTPaZ
for source zone coordinates in transport : write MTx +name of the zone and MTy...
the format of each line is the same as in the pst file
observation model input are not read from this file
they will be produced

2. from this file ipht3d produces the template files (.tpl) for pest
and a batch file (runmod.bat)

3. ipht3d reads files starting as pest_obs that contain the observed values
one line per observation zone, can contain layer ($l7)  name :pest_obs.txt
start with line headers that are the zone namen, time and observed variables names


4. ipht3d writes in the same folder a script to read data and write them in a correct
format in a txt file (out_..)

5. then pest is run by the user in a cmd window (using runmod.bat)
for each pest run, the python script is run"""

# all imports
import os,sys
from pylab import loadtxt,savetxt
from config import *

class Pest:
    def __init__(self,core):
    # get all important parameters from the current model (names and timelist)
        self.core = core
        cfg = Config(self.core)
        if core.gui != None:
            self.dialogs = cfg.dialogs
              
    def getDicBack1(self):
        '''creates a dic to choose the lines that will be used to vary parameters'''
        if self.core.dicaddin['Pback1'] != {}:
            return self.core.dicaddin['Pback1']
        else :
            dic = {'Modflow':{},'Mt3dms':{},'Pht3d':{}}
            for md in dic.keys() :
                l0 = self.core.dickword[md].lines.keys()
                dic[md] = zip(l0,[False]*len(l0))
            return dic

    def getDicBack2(self):
        if self.core.dicaddin['Pback2'] != {}:
            return self.core.dicaddin['Pback2']
        else :
            dic = {'Modflow':{},'Mt3dms':{},'Pht3d':{}}
            cols = ['Use','media','value','min','max','Group']
            for md in dic.keys() :
                #print 'pest 58 dict',self.core.dicaddin['Pback1'][md]
                a = self.core.dicaddin['Pback1'][md]
                lines = [x[0] for x in a if x[1]!=0]
                data = []
                for i,line in enumerate(lines):
                    val = float(self.core.dicval[md][line][0])
                    data.append([False,0,val,val/2,val*2,'G1'])
                dic[md] = {'cols':cols,'rows':lines,'data':data}
            return dic
            
    def getDicZones1(self):
        '''creates a dic to choose the lines that will be used to vary parameters'''
        if self.core.dicaddin['Pzones1'] != {}:
            return self.core.dicaddin['Pzones1']
        else :
            dic = {'Modflow':{},'Mt3dms':{},'Pht3d':{}}
            for md in dic.keys() :
                l0 = self.core.diczone[md].dic.keys()
                dic[md] = zip(l0,[False]*len(l0))
            return dic

    def getDicZones2(self):
        '''create a table from the selected zones (added 7/4/18)
        dicline is a dict for the lines the be used, may com from a dialog
        dicPzones contain the zones info to start Pest later'''
        dicout,dicline = {},self.core.dicaddin['Pzones1']
        if self.core.dicaddin['Pzones2'] != {}:
            dicPzones = self.core.dicaddin['Pzones2']
        else :
            dicPzones = {}
        #print dicPzones,dicline
        cols = ['Use','media','value','min','max','Transf','Group','UseX','value','min','max',\
            'UseY','value','min','max']
        for md in dicline.keys():
            a = dicline[md];
            llist = [x[0] for x in a if x[1]!=0] 
            for line in llist : 
                if line not in dicPzones.keys(): dicPzones[line] = {'rows':[]}
                dicz = self.core.diczone[md].dic[line]
                nzone = len(dicz['name'])
                dicout[line] = {'cols': cols,'rows': dicz['name'],'data':['a']*nzone}
                xmn,xmx,ymn,ymx = 1e8, -1e8, 1e8, -1e8
                for coo in dicz['coords']:
                    x,y = zip(*coo)
                    xmn,xmx = min(xmn,amin(x)),max(xmx,amax(x))
                    ymn,ymx = min(ymn,amin(y)),max(ymx,amax(y))
                dx, dy = (xmx-xmn)/100,(ymx-ymn)/100
                for i,n in enumerate(dicz['name']):
                    val = float(dicz['value'][i]);print line,n
                    x,y = zip(*dicz['coords'][i])
                    xm,ym = mean(x),mean(y)
                    if n not in dicPzones[line]['rows']: # a zone has been added
                        vmin,vmax,xmin,xmax,ymin,ymax = str(val/5),str(val*5),str(xm-dx)[:6],str(xm+dx)[:6],str(ym-dy)[:6],str(ym+dy)[:6]
                        lst = [False,dicz['media'][i],str(val),vmin,vmax,'none','G2']
                        lst.extend([False,str(xm)[:6],xmin,xmax]) # x coords
                        lst.extend([False,str(ym)[:6],ymin,ymax]) # y coords
                    else : # get the corresponding limits from dicin
                        i_old = dicPzones[line]['rows'].index(n)
                        lst = dicPzones[line]['data'][i_old]
                        lst[2],lst[8],lst[12] = val,str(xm)[:6],str(ym)[:6] # take current values of the zone
                    dicout[line]['data'][i] = lst
        #print dicout
        return dicout
                               
    def dic2parms(self,dicPback,dicPzones):
        '''transform the data in the dictionnary from the dialog to a format
        similar to the pest_parm file
        format for background ['MF,MT,PH]_line_bk_m_[media nb] (i media =-1 all)
        format for zones ['MF,MT,PH]_line_zo_[v,x,y]_zname
        if PH [k,s,p,m,u,i]'''
        self.pnames=[]; # list of parameter names
        self.ptrans=[]; # list of parameter transformation
        self.pchglim=[] #list of parameter change limit
        self.pvalbnd=[] #list of parameter values - lower bound - upper bound
        self.pgrp = [] #list of parameter group
        self.pparm=[] #list of parameter scale - offset - dercom
        self.ptied=[] #tied parameters
        prunstring ='' # the base to write pest_run.txt
        for md in dicPback.keys(): # here the keys are the models
            #if dicPback[md] != {}: pref = md[:2].upper()
            for i,line in enumerate(dicPback[md]['rows']):
                dp = dicPback[md]['data'][i]
                if dp[0]:
                    self.pnames.append(line+'_bm'+str(dp[1]))
                    self.pvalbnd.append(dp[2:5])
                    self.ptrans.append('none')
                    self.pgrp.append(dp[5])
        for line in dicPzones.keys(): # here the keys are the lines
            cols = dicPzones[line]['cols']
            i1,ix,iy = cols.index('Use'),cols.index('UseX'),cols.index('UseY')
            for i,zname in enumerate(dicPzones[line]['rows']):
                dp = dicPzones[line]['data'][i]
                if dp[i1]: #zone selected
                    self.pnames.append(line+'_zv'+zname)
                    self.pvalbnd.append(dp[i1+2:i1+5])
                    self.ptrans.append(dp[5]);self.pgrp.append(dp[6])
                if dp[ix]: #x zone selected
                    self.pnames.append(line+'_zx'+zname)
                    self.pvalbnd.append(dp[ix+1:ix+4])
                    self.ptrans.append(dp[5]);self.pgrp.append(dp[6])
                if dp[iy]: #y zone selected
                    self.pnames.append(line+'_zy'+zname)
                    self.pvalbnd.append(dp[iy+1:iy+4])
                    self.ptrans.append(dp[5]);self.pgrp.append(dp[6])                    
        for i in range(len(self.pnames)):
            self.pchglim.append('factor')
            self.pparm.append(['1','0','1']);self.ptied.append('')
            prunstring += self.pnames[i]+' '+str(self.pvalbnd[i][0])+'\n'
        self.nparm = len(self.pnames);#print self.pnames
        self.nparmgrp = len(unique(self.pgrp))
        f1=open(self.fdir+os.sep+'pest_run.txt','w');f1.write(prunstring);f1.close()

        
    def writeTpl(self):
        # write the template file
        self.prtMF,self.prtMT,self.prtPH=False,False,False
        s='ptf @ \n'
        for name in self.pnames:
            line = name.split('_')[0]
            s+=name+' @'+name+'@ \n'
            if line in self.core.dickword['Modflow'].lines: self.prtMF=True
            if line in self.core.dickword['Mt3dms'].lines: self.prtMT=True
            if line in self.core.dickword['Pht3d'].lines: self.prtPH=True
        f1=open(self.fdir+os.sep+'pest_tpl.txt','w')
        f1.write(s);f1.close()
        
    def writeBat(self):
    # produce the runmod.bat file (from ev)
        if self.core.dicval['Pest']['sys.1']==[0]:
            s='python scriptPest1.py \n'
            if self.prtMF: s+=self.mdir+'\\bin\\mf2k_Pmwin.exe '+self.fname+'\n'
            if self.prtMT:
                if 'VDF' in self.core.getUsedModulesList('Mt3dms'):
                    s+=self.mdir+'\\bin\\swt_v4.exe Mt3dms \n'
                else : s+=self.mdir+'\\bin\\mt3dms5b.exe Mt3dms \n'
            if self.prtPH: s+=self.mdir+'\\bin\\pht3dv217.exe Pht3d \n'
            s+='python scriptPest2.py'
            f1=open(self.fdir+os.sep+'runmod.bat','w')
            f1.write(s);f1.close()
            print 'runmod.bat written'
        if self.core.dicval['Pest']['sys.1']==[1]:
            s='#!/bin/sh \n'
            s+='## \n'
            s+='cd ../model \n'
            s+='## \n'
            s+='python scriptPest1.py \n'
            s+='## \n'
            if self.prtMF: s+='mf2k '+self.fname+'\n'
            if self.prtMT:
                if 'VDF' in self.core.getUsedModulesList('Mt3dms'):
                    s+='swtv4 Mt3dms \n'
                else : s+='mt3dms Mt3dms \n'
            if self.prtPH: s+='pht3dv217 Pht3d \n'
            s+='## \n'
            s+='python scriptPest2.py'
            f1=open(self.fdir+os.sep+'runmodel','w')
            f1.write(s);f1.close()
            print 'runmodel written'
        #else : ###### implement linux
        
                    
    def getObsPt(self):
    # get the observation files and gathers data in one dict from ev
        os.chdir(self.fdir)
        self.ospec=[]; # list of observed species = head tracer or chemical species
        self.onames=[]; # list of observation points or zones
        self.obs=[] #list that will contain the time and  value of observed data
        self.oweight=[] #dict that will contain the weight
        self.ogrp=[] #dict that will contain the group 
        f1=open('pest_obs.txt','r') #a file with all obs, 1st col obs point, 2nd time, then variables
        titl=f1.readline()
        for l in f1:
            l1 = l.split()
            self.ospec.append(l1[0])
            self.onames.append(l1[1])
            self.obs.append(l1[2:4])
            self.oweight.append(l1[4])
            self.ogrp.append(l1[5])
        self.nobs = len(self.onames)  # number of observation
        self.nobsgrp = len(unique(self.ogrp)) # number off observation group
        f1.close()
        print 'observation read'

    def writeInst(self):
    # write the instruction files (from ev)
        os.chdir(self.fdir)
        s='pif @ \n'
        n=1
        f1=open('pest.ins','w')       
        for i in range(len(self.onames)):
            self.ncol=len(self.obs[i])
            s+='l1 '
            s+='['+str(self.ospec[i][:2])+str(self.onames[i])+'_'+str(i)+']'+str(11)+':'+str(22)+' '
            s+=' \n';n+=1
        f1.write(s);f1.close()
        print 'instruction file written'

    def writePyscript(self):
    # writes the pyton script to write data before model run
        os.chdir(self.mdir+os.sep+'ilibq')
        f1=open('tplPestScript1.txt','r')
        s=f1.read();f1.close()
        os.chdir(self.fdir)
        if self.core.dicval['Pest']['sys.1']==[0]: #windows
            s=s.replace('ppmdir',self.mdir)
        else : 
            s=s.replace('ppmdir',str(self.core.dicval['Pest']['sys.2'][0]))
        s=s.replace('ppfdir',self.fdir)
        s=s.replace('ppfname',self.fname)
        f1=open(self.fdir+os.sep+'scriptPest1.py','w')
        f1.write(s);f1.close()
    # writes the pyton script to retrieve data after model run
        os.chdir(self.mdir+os.sep+'ilibq')
        f1=open('tplPestScript2.txt','r')
        s=f1.read();f1.close()
        os.chdir(self.fdir)
        if self.core.dicval['Pest']['sys.1']==[0]:
            s=s.replace('ppmdir',self.mdir)
        else : 
            s=s.replace('ppmdir',str(self.core.dicval['Pest']['sys.2'][0]))
        s=s.replace('ppfdir',self.fdir)
        s=s.replace('ppfname',self.fname)
        tlist=[]
        for i in range(len(self.obs)):
            tlist.append([self.obs[i][0]])
        s=s.replace('pptime',str(tlist))
        s=s.replace('pponames',str([str(a) for a in self.onames]))
        s=s.replace('ppospec',str(self.ospec))
        if self.core.dicval['Pest']['obs.1']==[0]:
            s=s.replace('ppobstrans',str("f1.write('%.11f '%(d1[it,0]))"))
        if self.core.dicval['Pest']['obs.1']==[1]:
            s=s.replace('ppobstrans',str("f1.write('%.11f '%(d1[it,0])) if d1[it,0] == 0.0 else f1.write('%.11f '%(log10(d1[it,0])))"))
        if self.core.dicval['Pest']['obs.1']==[2]:
            s=s.replace('ppobstrans',str("f1.write('%.11f '%(sqrt(abs(d1[it,0]))))"))
        f1=open(self.fdir+os.sep+'scriptPest2.py','w')
        f1.write(s);f1.close()
        print 'pyscripts written'
        
    def writePst(self):
    # Writes the pst file
        self.ntfiles,self.nifiles=1,1
    # Control Data Section
        s='pcf \n* control data \n'
        if self.core.dicval['Pest']['ctd.1'][0]==0:  ## l1
            s+= 'restart estimation\n'
        else : 
            s+= 'norestart estimation\n'
        s+=str(self.nparm)+' '+str(self.nobs)+' '+str(self.nparmgrp)+' 0 '+str(self.nobsgrp)+'\n' ## l2
        s+=str(self.ntfiles)+' '+str(self.nifiles)+' single point 1 0 0 \n' ## l3
        s+='10.0 2.0 0.3 0.03 10\n' ## l4 RLAMBDAl RLAMFAC !PHIRATSUF !PHIREDLAM !NUMLAM
        s+=' '.join([str(x) for x in self.core.dicval['Pest']['ctd.4']])+'\n' ## l5 RELPARMAX FACPARMAX FACORIG 
        s+=str(self.core.dicval['Pest']['ctd.5'][0])+'\n'  ## l6 PHIREDSWH
        s+=' '.join([str(x) for x in self.core.dicval['Pest']['ctd.6']])+'\n'  ## l7 NOPTMAX PHIREDSTP NPHISTP NPHINORED RELPARSTP NRELPAR
        s+='1 1 1 \n' # l8 ICOV ICOR IEIG
    # SVD section
        if self.core.dicval['Pest']['svd.1'][0]>0:
            s+='* singular value decomposition \n'
            s+=str(self.core.dicval['Pest']['svd.1'][0])+'\n' ## l1
            if self.core.dicval['Pest']['svd.2'][0]=='' :  ## l2
                SingVal=[x for x in range(len(self.ptrans)) if self.ptrans[x]!='tied']
                s+=str(len(SingVal))+' '+str(self.core.dicval['Pest']['svd.2'][1])+'\n'
            else : 
                s+=str(int(self.core.dicval['Pest']['svd.2'][0]))+' '+str(self.core.dicval['Pest']['svd.2'][1])+'\n'
            s+=str(self.core.dicval['Pest']['svd.3'][0])+'\n' ## l3
    # Parameter groups section
        s+='* parameter groups \n'
        parm_inc=['relative','absolute','rel_to_max']
        deriv=['switch','always_2','always_3','switch_5','always_5']
        central=['parabolic','outside_pts','best_fit','minvar','maxprec']
        lgrp = unique(self.pgrp)
        for lg in lgrp:  # nb of groups of parameters
            s+=str(lg)+' '
            s+=(str(parm_inc[self.core.dicval['Pest']['pgr.2'][1]])+' '+
                str(self.core.dicval['Pest']['pgr.2'][2])+' '+
                str(self.core.dicval['Pest']['pgr.2'][3])+' '+
                str(deriv[self.core.dicval['Pest']['pgr.2'][4]])+' '+
                str(self.core.dicval['Pest']['pgr.2'][5])+' '+
                str(central[self.core.dicval['Pest']['pgr.2'][6]])+'\n') ## INCTYP DERINC DERINCLB FORCEN DERINCMUL DERMTHD
    # Parameter data
        s+='* parameter data \n'
        for i in range(len(self.pnames)):
            s+=str(self.pnames[i])+'\t'+self.ptrans[i]+'\t'+self.pchglim[i]+'\t'
            s+=str(self.pvalbnd[i][0])+'\t'+str(self.pvalbnd[i][1])+'\t'+str(self.pvalbnd[i][2])+'\t'+self.pgrp[i]+'\t'
            s+=str(self.pparm[i][0])+'\t'+str(self.pparm[i][1])+'\t'+str(self.pparm[i][2])
            if self.ptied[i] !='' : 
                s+=str(self.ptied[i][0])+' '+str(self.ptied[i][1])+'\n'
            else : s+='\n'
    # Observation group
        s+='* observation groups \n'
        s+= '\n'.join(unique(self.ogrp))+'\n'
    # Observation data
        s+='* observation data \n'
        for i in range(len(self.onames)):
            s+=self.ospec[i][:2]+self.onames[i]+'_'+str(i)+' '
            if self.ospec[i] not in ['Head','Part']:
                if self.core.dicval['Pest']['obs.1']==[0]:
                    s+=str(float(self.obs[i][1]))+' '
                elif self.core.dicval['Pest']['obs.1']==[1]:
                    s+=str(log10(float(self.obs[i][1])))+' '
                elif self.core.dicval['Pest']['obs.1']==[2]:
                    s+=str(sqrt(float(self.obs[i][1])))+' ' 
            else :
                s+=str(float(self.obs[i][1]))+' '
            s+=self.oweight[i]+' '
            s+=self.ogrp[i]+'\n'
    # Model command line & model input/output
        if self.core.dicval['Pest']['sys.1']==[0]:
            s+='* model command line \n runmod.bat \n'
            s+='* model input/output \n'
            s+='pest_tpl.txt  pest_run.txt \n'
            s+= 'pest.ins pest_out.txt \n'
        else :
            s+='* model command line \n./runmodel \n'
            s+='* model input/output \n'
            s+='tpl/pest_tpl.txt  ../model/pest_run.txt \n'
            s+= 'ins/pest.ins ../model/pest_out.txt \n'
        f1=open(self.fdir+os.sep+self.fname+'.pst','w')
        f1.write(s);f1.close()
        print 'Pst file written'
        
    def writeRegPst(self):
        self.Reg = False
        if self.core.dicval['Pest']['ctd.1'][1]==2:
            s = self.mdir+os.sep+'bin'+os.sep+'addreg1 '+self.fname+'.pst '+self.fname+'r.pst'
            os.system(s)#+ " & pause")
            self.Reg = True
            f1=open(self.fname+'r.pst','r')
            s=f1.read();f1.close()
            reg1 = (str(self.core.dicval['Pest']['reg.1'][0])+' '+
                    str(self.core.dicval['Pest']['reg.1'][1])+' '+
                    str(self.core.dicval['Pest']['reg.1'][2]))
            if self.core.dicval['Pest']['reg.1'][3] == 1 : reg1+=' memsave'
            if self.core.dicval['Pest']['reg.1'][3] == 2 : reg1+=' nomemsave'
            s=s.replace('1.0000000E-10  1.0500000E-10  0.1000000',reg1)
            reg2 = (str(self.core.dicval['Pest']['reg.2'][0])+' '+
                    str(self.core.dicval['Pest']['reg.2'][1])+' '+
                    str(self.core.dicval['Pest']['reg.2'][2]))
            if self.core.dicval['Pest']['reg.2'][3] == 1 : reg2+=' linreg'
            if self.core.dicval['Pest']['reg.2'][3] == 2 : reg2+=' nonlinreg'
            if self.core.dicval['Pest']['reg.2'][4] == 1 : reg2+=' continue'
            if self.core.dicval['Pest']['reg.2'][4] == 2 : reg2+=' nocontinue' 
            s=s.replace('1.0   1.0e-10    1.0e10',reg2)
            reg3 = (str(self.core.dicval['Pest']['reg.3'][0])+' '+
                    str(self.core.dicval['Pest']['reg.3'][1])+' '+
                    str(self.core.dicval['Pest']['reg.3'][2])+' ')
            if self.core.dicval['Pest']['reg.3'][3] != '' :
                reg3+=str(int(self.core.dicval['Pest']['reg.3'][3]))+' '
            if self.core.dicval['Pest']['reg.3'][3] != '' :
                reg3+=str(float(self.core.dicval['Pest']['reg.3'][4]))+' '
            if self.core.dicval['Pest']['reg.3'][4] != '':
                reg3+=str(float(self.core.dicval['Pest']['reg.3'][5]))
            s=s.replace('1.3   1.0e-2     1',reg3)
            f1=open(self.fname+'r.pst','w')
            f1.write(s);f1.close()
            print 'Reg pst file written'

    def writeFiles(self):
        d0,self.fname = self.core.fileDir,self.core.fileName
        d1=d0.split(os.sep)
        self.fdir = '//'.join(d1)
        d0=self.core.baseDir # main directory
        d1=d0.split(os.sep)
        self.mdir = '//'.join(d1)
        sys.path.append(self.mdir)

        self.dic2parms(self.core.dicaddin['Pback2'],self.core.dicaddin['Pzones2'])
        self.writeTpl()
        self.writeBat()
        self.getObsPt()
        self.writeInst()
        self.writePyscript()
        self.writePst()
        self.writeRegPst()

