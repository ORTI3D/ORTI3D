# -*- coding: utf-8 -*-
"""
Created on Sun Oct 19 12:54:49 2014

@author: olive
"""
from array import array as arr2
import os
from pylab import savetxt,loadtxt
from .geometry import *
from .timeperiod import *

class min3pWriter:
    
    def __init__(self, core,fDir, fName):
        self.core = core
        self.fDir,self.fName = fDir,fName
        self.fullPath = fDir+os.sep+fName;#print self.fullPath
        self.fileKeys = {'poro.1':['porosity field','por'], # to write in distrib files
                'flow.1':['hydraulic conductivity field','hyc'],
                'flow.2':['specific storage coefficient','spstor'],
                'flow.4':['hydraulic conductivity field','hyc'], # although it is permeability here
                'engp.3':['energy balance parameters','energybal'],
                }
        self.addKey = {}
        
    def writeMin3pFiles(self, core, option):
        #option Flow, Trans or Chem
        self.core,self.min3p = core,core.addin.min3p
        self.chem = self.core.addin.chem
        if option == 'Trans' and core.dicval['Min3pFlow']['glo.1'][6]==1: 
            option = 'Heat' #in transport is engery balcne is active, this is Heat
        if core.addin.getDim()=='Radial':
            self.addKey['spat']=['radial coordinates']
        f1=open(self.fullPath +'.dat','w')
        self.ttable = core.makeTtable();#print 'mfw',self.ttable
        tlist = array(self.ttable['tlist'])
        self.per = tlist[1:]-tlist[:-1]
        self.nper = len(self.per)
        if core.getValueFromName('Min3pFlow','Fsteady')==0: self.nper=1
        self.min3p.nper = self.nper
        core.Zblock = makeZblock(core);#print 'mfw,zb',shape(self.core.Zblock)
        self.xsect = False
        if core.addin.getDim() in ['Xsection','Radial']: self.xsect = True
        self.domn=[]
        for n in ['Xmin','Xmax','Ymin','Ymax','Zmin','Zmax']:
            self.domn.append(core.getValueFromName('Min3pFlow',n))
        self.grid = core.addin.getFullGrid()
        s=''
        lgrp=['glo','spat','time','out','conf','poro','flow','inif','bcf']
        if core.getValueFromName('Min3pFlow','use_evapo')==1: lgrp.append('conv')
        if option == 'Trans': lgrp.extend(['cont','trans','trac','init','bct'])
        elif option == 'Heat': lgrp.extend(['cont','cone','trans','trac','init','bct','engp','inie','bce'])
        elif option == 'Chem': lgrp.extend(['cont','conc','trans','inic','bcc'])
        Fkey,Tkey,Ckey = core.dickword['Min3pFlow'],core.dickword['Min3pTrans'],core.dickword['Min3pChem']
        for grp in lgrp:
            #print 'group',grp
            if grp in Fkey.grpList:
                if grp not in self.core.addin.getUsedModulesList('Min3pFlow'):continue
                self.dicval = core.dicval['Min3pFlow']
                s+='!-------------------------\n\''+Fkey.longNames[grp]+'\'\n'
                s+=self.writeBlock('Min3pFlow',Fkey,grp,Fkey.groups[grp])
            elif grp in Tkey.grpList:
                self.dicval = core.dicval['Min3pTrans']
                s+='!-------------------------\n\''+Tkey.longNames[grp]+'\'\n'
                s+=self.writeBlock('Min3pTrans',Tkey,grp,Tkey.groups[grp])
            elif grp in Ckey.grpList:
                self.dicval = core.dicval['Min3pChem']
                s+='!-------------------------\n\''+Ckey.longNames[grp]+'\'\n'
                s+=self.writeBlock('Min3pChem',Ckey,grp,Ckey.groups[grp])
        if option == 'Chem': 
            s+= self.writeGeochem()
            self.writeDatabases()
        f1.write(s)
            
    def writeBlock(self,mod,Dict,grp,llist):
        if grp == 'spat':
            s = self.writeSpat()
        elif grp in ['glo','time','out','trac','conf','cont','cone','conc','conv']: 
            s = self.writeGeneral(mod,Dict,grp,llist)
        elif grp[:4] == 'poro':
            s = self.writeMedia(mod,Dict,llist)
        elif grp[:3] == 'ini':
            s = self.writeIni(mod,Dict,llist)
        elif grp[:2] == 'bc':
            s = self.writeBcs(mod,Dict,llist)
        else: # concerns only flow1..4 and engp3
            s = self.writeByZones(mod,Dict,llist)
        return s
            
    def writeGeneral(self,mod,Dict,grp,llist):
        s = '';#print llist
        for l in llist:
            line = l*1;#print line
#            if self.xsect and l == 'spat.2': line = 'spat.3' # for coords inversion
#            if self.xsect and l == 'spat.3': line = 'spat.2'
            cond=Dict.lines[line]['cond'];#print line
            if self.core.testCondition(mod,cond)==False : continue
            kwlist=Dict.lines[line]['kw']
            ktyp=Dict.lines[line]['type'];#print 'mfw99',line,Dict.lines[line]
            lval=self.dicval[line]
            name = Dict.lines[line]['comm']
            name=name.split('(')[0]
            # try to find the keywords that are just one line
            lprint,oneline = True,False; #print line,lval,Dict.lines[line]['detail'][0]
            if Dict.lines[line]['type'][0]=='choice':
                if Dict.lines[line]['detail'][0][0] == 'use': oneline = True
                if Dict.lines[line]['detail'][0][lval[0]+1]=='no': 
                    lprint = False # oa 26/5 do not print if there is only one line and no
            #print line,name,Dict.lines[line],lprint,oneline
            s += '\n'
            if (line[:4] not in ['spat','time']) and lprint: 
                s += '\''+name+'\'\n'
            for ik in range(len(kwlist)):
                #print line,lval
                value=lval[ik];#
                if ktyp[ik]=='choice': # where there is a choice print the nb of he choice not value
                    detail = Dict.lines[line]['detail'][ik]
                    choi = detail[value+1]# +1 because the first line is the title
                    if choi in ['true','false']: sep='.' # for the global keywords at the start
                    elif choi == 'yes': 
                        if oneline : choi,sep = '', ''
                        else : choi, sep = Dict.lines[line]['detail'][ik][0],'\'' 
                    elif choi == 'no': 
                        choi,sep = '',''
                    elif 'print' in Dict.lines[line]:
                        choi,sep = '\''+detail[0]+'\'\n\''+detail[value+1]+'\'\n',''
                    else :
                        sep =  '\''
                    if choi+sep != '': s += sep+choi+sep+'         ;'+detail[0]+'\n' # oa 26/5
                elif ktyp[ik]=='title': # case of a title line
                    s += '#'+str(value)+'\n'
                else : # values with strings
                    if line in ['trac.2']: 
                        value = '\''+self.core.gui.mainDir+os.sep+'utils\''
                    if line == 'out.1':
                        value = self.getTimeList()
                    s += str(value)+'\n'
        if grp in self.addKey:
            for n in self.addKey[grp]: s+='\''+n+'\'\n'
        return s+'\n\'done\'\n\n'
        
    def writeSpat(self):  # OA added 24/5
        '''writes spatial data'''
        def str1dim(x0,x1,nx,dx):
            if len(unique(dx))==1: # regular grid
                s1 = '1\n'+str(nx)+'\n'+str(x0)+' '+str(x1)+'\n'
            else :
                s1 = str(nx)+'\n' # nx control intervals and nx volumes
                cumul = x0
                for i in range(len(dx)):
                    s1 += '1\n  %4g  %4g' %(cumul,cumul+dx[i]) +'\n'
                    cumul += dx[i]
            return s1
        g = self.grid
        s = str1dim(g['x0'],g['x1'],g['nx'],g['dx'])
        s2 = str1dim(g['y0'],g['y1'],g['ny'],g['dy'])
        if self.xsect==False: s += s2 + '\n1\n1\n0.0 1.0\n'
        else : s += '\n1\n1\n0.0 1.0\n' + s2
        return s+'\n\'done\'\n\n'
        
    def writeMedia(self,mod,Dict,llist):
        llist1 = ['poro.1','flow.1','flow.2','engp.3']
        nbztot = self.getNbZonesMedia(mod,llist1) # total nb of zones for list of lines - domain, starts from 0
        s = str(nbztot+1)+'\n'
        s += '\'number and name of zone\'\n1 \n'
        s += '\'domain\'\n'
        poro_base = self.getValue(mod,llist1[0],None,0,'whole')
        s += poro_base +'\n'
        s += '\'extent of zone\'\n'
        s += self.getDomainCoords()
        s += '\n\'end of zone\'\n\n'
        lznames,lzcoords = [],[]  # OA 24/5
        for line in llist1: #then write zones
            if line=='engp.3': mod='Min3pTrans'
            if line not in list(self.core.diczone[mod].dic.keys()): continue # no zones
            dicz = self.core.diczone[mod].dic[line]
            nbzline = len(dicz['value']);#print line, nbzline,self.allRectZones(dicz)# nb of zones for this line
            nbz=1 # OA 31/05
            if self.allRectZones(dicz) or self.min3p.nodes!=None: # write only for rect zone or unstructured
                for iz in range(nbzline): # OA 23/5
                    zname = dicz['name'][iz]
                    if (zname == 'domain') or (zname in lznames): continue # OA 23/5 added lznames on 24/5
                    nbz += 1 # OA 23/5
                    lznames.append(zname) # OA 24/5 remember names
                    s += '\'number and name of zone\'\n'+ str(nbz) +'\n'
                    s += '\''+zname+'\''  + '\n'
                    if line[:4]=='poro': s += self.getValue(mod,line,dicz,iz,'zone')
                    else : s += poro_base
                    s += '\'extent of zone'
                    lzcoords.append(self.getCoords(dicz['coords'][iz],line))
                    s +=  lzcoords[-1]+'\n'
                    s += '\'end of zone\'\n\n'
        self.lznames,self.lzcoords = lznames,lzcoords # OA 25/5 remember the list of zones
        return s+'\'done\'\n'
        
    def writeIni(self,mod,Dict,llist):
        nbztot = self.getNbZonesMedia(mod,llist) # nb of zones for this ini
        line=llist[0]
        s = str(nbztot+1)+'\n'
        s += '\'number and name of zone\'\n1 \n'
        s += '\'domain\'\n \'initial condition\'\n'
        s += self.getValue(mod,line,None,0,'whole')
        s += '\'extent of zone\'\n'
        s += self.getDomainCoords()
        s += '\n\'end of zone\'\n\n'
        if nbztot ==0: return s+'\'done\'\n'
        dicz = self.core.diczone[mod].dic[line]
        nbzline = len(dicz['value']);#print line, nbzline,self.allRectZones(dicz)# nb of zones for this line
        if self.allRectZones(dicz) or self.min3p.nodes!=None: # write only for rect zone or unstructured
            for iz in range(nbzline): # OA 23/5
                zname = dicz['name'][iz]
                s += '\'number and name of zone\'\n'+ str(iz+2) +'\n'
                s += '\''+zname+'\''  + '\n \'initial condition\'\n'
                s += self.getValue(mod,line,dicz,iz,'zone')
                s += '\'extent of zone'
                s += self.getCoords(dicz['coords'][iz],line)+'\n'
                s += '\'end of zone\'\n\n'
        return s+'\'done\'\n'
               
    def writeByZones(self,mod,Dict,llist):
        '''writes spatial data several situations appear :
            this is written as a line list because MIn3p requires all paramters in the zones
        - some lines here are not as an array format, so just one number shall be written
        - for physical conditions if zones are present, it is written directly in distributed files, 
           except some lines that are not arrays
        '''
        nbztot = len(self.lznames);# total nb of zones in porous medium
        s= str(nbztot+1)+'\n' # starts at 0 
        for line in llist: # don't write as spatial data if not array, for engp
            if Dict.lines[line]['type'][0][:3] != 'arr': 
                s += '\''+Dict.lines[line]['comm'].split('(')[0]+'\'\n'
                s += self.getValue(mod,line,None,0,'whole')
                
        s += '\'number and name of zone\'\n1 \n\'domain\'\n'
        for line in llist:
            if Dict.lines[line]['type'][0][:3] != 'arr':  continue
            cond=Dict.lines[line]['cond'];#the conditio to test if the line shall be printed
            if self.core.testCondition(mod,cond)==False : continue
            info = Dict.lines[line]['comm'].split('(')[0]
            s += '\''+info+'\'' + '\n'
            s += self.getValue(mod,line,None,0,'whole')
            if line == 'trans.1': s+= self.getDiffusion()
        s += '\'extent of zone\'\n'+self.getDomainCoords()+'\n'
        s += '\'end of zone\'\n\n'

        flagDis = False #write read from file if necessary
        if self.min3p.nodes==None: # read from file only for structured grid
            for line in llist: 
                if (line not in list(self.core.diczone[mod].dic.keys())): continue
                dicz = self.core.diczone[mod].dic[line]
                if self.allRectZones(dicz)==False: # or self.core.dictype[mod][0]=='formula': # spatial data can be stored in a file if they contain zones
                    self.writeDistributedFile(line,self.getFileKwd(line,1))
                    s += '\'read '+self.getFileKwd(line,0)+' from file\'\n'
                    flagDis = True
            if flagDis: return s+'\'done\'\n\n'
        #print 'nbztot',nbztot
        for i,znam in enumerate(self.lznames):#then write zones, with all properties zones being the same as in poro
            s += '\'number and name of zone\'\n'+ str(i+2) +'\n';#OA 26/5 the first zone is domain
            s += '\''+znam+'\''  + '\n'
            dictop = self.core.diczone[mod].dic
            for line in llist: 
                if Dict.lines[line]['type'][0][:3] != 'arr':  continue
                cond=Dict.lines[line]['cond'];#the conditio to test if the line shall be printed
                if self.core.testCondition(mod,cond)==False : continue
                info = Dict.lines[line]['comm'].split('(')[0]
                if line not in list(dictop.keys()):  # no zones in that line
                    s += '\''+info+'\'\n'+ self.getValue(mod,line,None,0,'whole')                    
                    continue
                elif znam not in dictop[line]['name']: #there are zones but not the right one
                    s += '\''+info+'\'\n'+ self.getValue(mod,line,None,0,'whole')                    
                    continue
                else : # the zone is there                    
                    iz = dictop[line]['name'].index(znam)
                    s += '\''+info+'\'\n'+ self.getValue(mod,line,dictop[line],iz,'zone')
            s += '\'extent of zone'
            s += self.lzcoords[i] +'\n'
            s += '\'end of zone\'\n\n'

        return s+'\'done\'\n\n'
        
    def writeBcs(self,mod,Dict,llist):
        nbz,nbztot = 0,self.getNbZonesBCs(mod,llist)   #normally BC zoens are "rect" (in fact line)
        s = str(nbztot)+'\n'
        for line in llist: 
            if line not in list(self.core.diczone[mod].dic.keys()): continue # no zones
            dicz = self.core.diczone[mod].dic[line]
            nbzline = len(dicz['value'])# nb of zones for this line
            for iz in range(nbzline):
                nbz += 1
                s += '\'number and name of zone\'\n'+ str(nbz) +'\n'
                s += '\''+dicz['name'][iz]+'\''  + '\n'
                # for this zone write type of information and value
                info = Dict.lines[line]['comm'].split('(')[0]
                s += '\'boundary type\'\n'
                s += '\''+info+'\''  + '\n'
                s += self.getValue(mod,line,dicz,iz,'zone')
                # for this zone write extent
                s += '\'extent of zone'
                s += self.getCoords(dicz['coords'][iz],line,'boundary') +'\n'
                s += '\'end of zone\'\n\n'
        s += '\'done\'\n\n'
        return s
        
    def getValue(self,mod,line,dicz,iz,opt):
        """get the values for zones or for the wole domain
        allow to get the chemistry from solutions too
        """
        if mod =='Min3pFlow': Dkeys = self.core.dickword['Min3pFlow'] # added 26/5 because we need to find lines in any dict
        elif mod =='Min3pTrans': Dkeys = self.core.dickword['Min3pTrans']
        elif mod =='Min3pChem': Dkeys = self.core.dickword['Min3pChem']
        prefx=''
        if 'prefx' in Dkeys.lines[line]: prefx = Dkeys.lines[line]['prefx']
        suffx=''
        if 'suffx' in Dkeys.lines[line]: suffx = Dkeys.lines[line]['suffx']
        s = ''
        dicz = self.core.diczone[mod].dic
        # several variables, work only with zones
        #print line, Dkeys.lines[line]
        if 'names' in Dkeys.lines[line] and line in dicz: # several variables in one line, with same zone
            vals = dicz[line]['value'][iz].split('$')[1].split('\n')
            longnames = Dkeys.lines[line]['longNames'];#print line,vals,longnames
            for i in range(len(longnames)):
                lname, a = longnames[i],''
                if lname != '': a= '\''+ lname +'\'\n'
                s+=a+' %+11.4e' %float(vals[i])+'\n'
        else :
            if opt=='whole':
                vlist = self.core.dicval[mod][line]
            elif opt=='zone':
                vlist = [dicz[line]['value'][iz]]
            for v in vlist: # to write all keywords
                if line in ['bcc.1','bcc.4','inic.1']: # OA 23/5 pb for inic not read
                    if opt == 'whole': nbsolu = int(v)  # OA 25/5 0 -> v because nbsolu can be in background
                    else: nbsolu = int(dicz[line]['value'][iz])
                    v = self.getStringFromSolution(line,nbsolu) # here, this is the background
                else :
                    v = ' %+11.4e' %float(v)
                    v = v.replace('e','d')
                s += prefx + v +' '+ suffx +'\n'
        return s
        
    def getStringFromSolution(self,line,value):
        """provide the string for the solution (or mineral) composition
        given by value as 1000 (like for pht3d) or 1210 for initial conditions
        an dsolution number in boundary cond
        """
        if line[:4] == 'inic':
            nsol,nmin,nexch,nsurf = value/1000,mod(value,1000)/100,mod(value,100)/10,mod(value,10)
        else :
            nsol = value
        #lgrp = ['comp','mineral','sorption','sorption']
        #lnames = ['concentration input','mineral input','sorption parameter input','sorption parameter input']
        # solutions
        dChem = self.core.addin.chem.Base['MChemistry']['comp']
        s = ''
        #if line[:3] == 'bcc': #OA removed 23/5 needed ofr inic too
        s = '\'concentration input\' \n'
        for ir,row in enumerate(dChem['rows']):
            if dChem['data'][ir][0]: #the species is ticked
                s += str(dChem['data'][ir][nsol+2]).replace('e','d')+'  \''+\
                    dChem['data'][ir][-1]+'\'  ;'+row+'\n' # last columns contains 'free', 'ph'...
        if line[:3] == 'bcc': return s
        # linear sorption
        dChem = self.core.addin.chem.Base['MChemistry']['linear sorption']
        sl = ''
        for ir,row in enumerate(dChem['rows']):
            if dChem['data'][ir][0]: #the species is ticked
                sl += '\''+row+'\'  '+str(dChem['data'][ir][nsol+1])+' \n' 
        if sl != '': s += '\'linear sorption input\' \n' + sl +'\n'
        
        # ion exchange # added OA 23/5 modif 21/6
        dChem = self.core.addin.chem.Base['MChemistry']['exchange']
        se = '';#print 'nexch',nexch
        if dChem['data'][0][0]:
            se += str(dChem['data'][0][nexch+1])+' \n' # cec of the right exchanger
            se += str(dChem['data'][0][1])+' \n' # rho
        if se != '': 
            s += '\n\'sorption parameter input\' \n' + se +'\n'      
            s += '\'equilibrate with fixed solution composition\'\n'

        # surface parameters OA added 14/06
        dChem = self.core.addin.chem.Base['MChemistry']['sorption'];#print dChem
        ss = '';
        if len(dChem['rows'])>0: # 
            for ir,row in enumerate(dChem['rows']):
                if dChem['data'][ir][0]: #the species is ticked
                    data = dChem['data'][ir]
                    ss += '\''+data[1]+'\'  '+ str(data[nsurf+2]).replace('e','d')+' ' # name, mass
                    ss += ' '.join([str(a) for a in data[-2:]])+' ; '+dChem['rows'][ir]+' name, mass, area,density\n'
        if ss != '': 
            s += '\n\'sorption parameter input\' \n' + ss +'\n'      
            s += '\'equilibrate with fixed solution composition\'\n'

        #minerals
        dChem = self.core.addin.chem.Base['MChemistry']['mineral']
        sm = '' # 25/5 OA not to print minerals if nothing
        for ir,row in enumerate(dChem['rows']):
            if dChem['data'][ir][0]: #the species is ticked
                data = dChem['data'][ir]
                sm += str(data[nmin+1]).replace('e','d')
                sm += '  '+data[5]+' \''+data[6]+'\' ; '+dChem['rows'][ir]+'\n'
                s1 = ' '.join([str(a) for a in data[7:]])
                sm += s1.replace('e','d')+'\n'
        if sm !='':
            s += '\n\'mineral input\' \n'+sm
        return s
        
    def getNbZonesMedia(self,mod,llist):
        '''returns the total number of zones for a list of lines, only counting
        the zones that are rectangles and that are different from the domain
        and the names that are different'''
        nbz,names = 0,[] # OA 24/5 added names
        dicz = self.core.diczone[mod]
        for line in llist:
            if line in dicz.dic:
                dz1 = dicz.dic[line]
                if self.allRectZones(dz1) != True:continue # OA 22/6
                for iz in range(len (dz1['coords'])):
                    n = dz1['name'][iz]
                    if (n=='domain') or (n in names): continue # OA 24/5
                    if (self.isRectZone(dz1['coords'][iz])):
                        nbz += 1 # nb of zones
                        names.append(n) # OA 24/5
        return nbz
        
    def getNbZonesBCs(self,mod,llist):
        '''returns the total number of zones for a list of lines, only counting
        the zones that are not rectangles and that are different from the domain'''
        nbz=0
        dicz = self.core.diczone[mod]
        for line in llist:
            if line in dicz.dic:
                nbz += dicz.getNbZones(line)
        return nbz

    def getDomainCoords(self):
        s = ''
        lst=list(range(6))
        for i in lst: s+= str(self.domn[i])+'  '
        return s

    def getCoords(self,coolist,line,opt='all'):
        '''returns a formatted string of coords from the orti list, these orti coords are in 
        x and y , they are transformed if Xsection'''
        nodes = self.min3p.nodes
        if self.isRectZone(coolist): #squared zones
            s = '\'\n'
            coolist = self.coord2BC(coolist) # OA 25/5 removed the if, may apply to all
            xl,yl = list(zip(*coolist));#print 'bc l 273',xl,yl
            s0 = str(min(xl))+' '+str(max(xl))+' '
            s1 = str(min(yl))+' '+str(max(yl))+' '
            if self.xsect : s += s0+str(self.domn[2])+'  '+str(self.domn[3])+'  '+s1
            else : s += s0+s1+str(self.domn[4])+'  '+str(self.domn[5])+'  '
        else : # usntructured non squared
            s = ': '+opt+' nodes in polygon\'\n\'read data\'\n'
            xl,yl = list(zip(*coolist));#print 'bc l 273',xl,yl
            npts = len(xl)
            s += str(npts)+'\n'
            zl = [0]*npts
            if self.xsect : zl = yl*1; yl = [0]*npts
            s += '\n'.join([str(xl[i])+' '+str(yl[i])+' '+str(zl[i]) for i in range(npts)])
        return s
        
    def allRectZones(self,dicz):
        '''test if all zoness in a list are rectangles'''
        coords = dicz['coords']
        rect = True
        for iz in range(len(coords)):
            rect = self.isRectZone(coords[iz])
            if rect != True: break
        return rect
        
    def isRectZone(self,coolist):
        rect = True
        x,y = list(zip(*coolist))
        if len(unique(x))>2 or len(unique(y))>2: rect = False
        return rect
           
    def coord2BC(self,coolist):
        dx,dy = float(self.grid['dx'][0]),float(self.grid['dy'][0])
        coo2 = [0]*len(coolist)
        cooy = self.domn[2:4]
        if self.xsect: cooy = self.domn[4:6]
        for i,c in enumerate(coolist):
            x,y = c; 
            if abs(x-self.domn[0])<dx/2 : x = self.domn[0]
            if abs(x-self.domn[1])<dx/2 : x = self.domn[1]
            if abs(y-cooy[0])<dy/2 : y = cooy[0]
            if abs(y-cooy[1])<dy/2 : y = cooy[1]
            coo2[i]=(x,y)
        #print 'mpw l408', self.domn,cooy,coolist,coo2
        return coo2
        
    def getTimeList(self):
        tlist = self.core.getTlist2()
        s = str(len(tlist))+'\n'
        for i in range(len(tlist)):
            if mod(i,4)==0: s+='\n'
            s += str(tlist[i])+' '
        return s
        
    def getFileKwd(self,line,i):
        return self.fileKeys[line][i]
        

    def writeDistributedFile(self,line,extension):
        """ to write a file with values distributed on the grid, 
        get the object from zone and call write function for external file"""
        f1=open(self.fullPath+'.'+extension,'w')
        s = 'title 1 \n title 2 \n title 3 \n'
        if line[:4] in ['trans','engp']: model = 'Min3pTrans'
        else : model = 'Min3pFlow'
        tblZones = self.getTableZones(model,line)
        nz,nv = shape(tblZones)
        if self.min3p.nodes != None: #unstructured
            nodvalues = zone2mesh(self.core,model,line)
            mat = c_[self.min3p.nodes[:,1:],tblZones[list(nodvalues),:]]
            s += arr2string1(mat)
        else : # structured
            coords =  getMesh3D(self.core) # in x,y,z order
            mat = [self.core.getValueLong(model,line,0)]
            m0 = self.core.getValueLong(model,line,0) # matrix contianing zone numbers
            mat = []
            for iv in range(nv):
                m1 = m0*0+tblZones[0,iv]
                for j in range(1,nz): m1[m0==j] = tblZones[j,iv]
                mat.append(m1)
            s += self.writeExtFile(coords,mat)
        f1.write(s)
        f1.close()
        
    def getTableZones(self,model,line):
        '''reads all zones and transforms the values in a table having nzones lines
        and nvalues columns'''
        dicz = self.core.diczone[model].dic[line]
        nz = len(dicz['value'])
        a = dicz['value'][0].split('$')[1].split('\n');a.remove('')
        nval = len(a)
        tbl = zeros((nz,nval))
        for iz in range(nz):
            a = dicz['value'][iz].split('$')[1].split('\n');a.remove('');#print a
            tbl[iz,:]=[float(x) for x in a]
        #print tbl
        return tbl
        
    def writeExtFile(self,coords,mat):
        # formats and write the object
        nvar,s = len(mat),'';   #print shape(mat),shape(coords[0])
        nz,ny,nx = shape(mat[0])
        X,Y,Z = coords
        if self.core.addin.getDim() in ['3D','2D']:
            for iz in range(nz):
                for iy in range(ny):
                    for ix in range(nx): 
                        s += str(X[iz,iy,ix])+' '+str(Y[iz,iy,ix])+' '+str(Z[iz,iy,ix])
                        for iv in range(nvar): 
                            s += ' '+str(mat[iv][iz,iy,ix])
                        s += '\n'
        else :# Xsection or radial
            for iz in range(nz):
                for ix in range(nx): 
                    s += str(X[iz,0,ix])+' 0.0 '+str(Z[iz,0,ix])
                    for iv in range(nvar): 
                        s += ' '+str(mat[iv][iz,0,ix])
                    s += '\n'
        return s
        
    def writeGeochem(self):
        #creates the geochemical system block from chem database
        # ! redox are only considered as kinetics up to now
        lgrp = ['complex','redox','gases','sorption','mineral']
        lnames = ['secondary aqueous species','intra-aqueous kinetic reactions','gases',
                  'sorbed species','minerals']
        s = '!-----------------------------\n\'geochemical system\'\n'
        s += '\'use new database format\'\n\n\'database directory\'\n'
        s += '\''+self.fDir+'\'\n\n'
        # add components (modified for immobile compounds)
        base = self.chem.Base['MChemistry']['comp']
        s1,s2 = [],[]
        for i,n in enumerate(base['rows']):
            if base['data'][i][0] and base['data'][i][1] : s2.append(n)
            elif base['data'][i][0] and not base['data'][i][1] : s1.append(n)
        s += '\'components\'\n' +str(len(s1))+'\n'+'\n'.join(s1)+'\n\n'
        if len(s2)>0: 
            s += '\'biomass components\'\n' +str(len(s2))+'\n'+'\n'.join(s2)+'\n\n'
        #others
        for a in zip(lgrp,lnames):
            group,name = a;#print a,group,name
            l = self.chem.getListSpecies(group)
            if len(l)>0 :
                s += '\''+name + '\'\n' # name of group
                s += str(len(l)) +'\n\'' # nb of solutes
                s += '\'\n\''.join(l)
                s += '\'\n\n'
        # add linear sorption
        l = self.chem.getListSpecies('linear sorption')
        if len(l)>0:
            s += '\'linear sorption\'\n' +str(len(l))+'\n'
            for n in l:
                s += '\''+n+'\' \n'
            s += '\n'
        # add ion exchange
        l,s1 = self.chem.getListSpecies('exchange'),''
        if len(l)>0:
            for n in l: 
                if n != '-x': s += '\''+n+'\' \n'
            s += '\'sorbed species\'\n' +str(len(l)-1)+'\n' +s1 +'\n'
        # add non aqueous conc for surface (it's name)
        l = self.chem.getListSpecies('sorption')
        s1,nb = '',0
        if len(l)>0:
            s+='\'define sorption type\' \n\'surface-complex\'\n\n'
            dChem = self.chem.Base['MChemistry']['sorption']
            for n in l:
                if n in dChem['text']:
                    indx = dChem['rows'].index(n)
                    nb += 1
                    s1 += '\''+dChem['data'][indx][1]+'\' \'surface\'\n'  
            s += '\'non-aqueous components\'\n' +str(nb)+'\n'+ s1 + '\n'
        # add scaling for kinetics
        l = self.chem.getListSpecies('redox')
        if len(l)>0:
            s += '\'scaling for intra-aqueous kinetic reactions\'\n'
            dChem = self.chem.Base['MChemistry']['redox']
            for n in l:
                indx = dChem['rows'].index(n)
                s += dChem['data'][indx][1] +'\n'
            s += '\n'
        # add gas decay
        l = self.chem.getListSpecies('gases')
        if len(l)>0:
            s += '\'gases decay rate constant\'\n'
            dChem = self.chem.Base['MChemistry']['gases']
            for n in l:
                indx = dChem['rows'].index(n)
                s += dChem['data'][indx][1] +'\n'
        return s+'\n\'done\'\n'
        
    def getDiffusion(self):
        diffChoice = self.core.getValueFromName('Min3pTrans','Diff_choice')
        if diffChoice==0: # no dusty gas
            return ''
        elif diffChoice==1 : s = '\n\'binary gas diffusion coefficients\' \n'
        elif diffChoice==2 : s = '\n\'dusty gas model\' \n'
        base = self.chem.Base['MChemistry']['gases']
        glist = base['rows']
        ng = len(glist)
        if diffChoice==1:
            for ig in range(ng):
                for i2 in range(ig+1,ng):
                    diff = self.calcDij(base['data'][ig][3:],base['data'][i2][3:])
                    suff = '  ;'+glist[ig]+'-'+glist[i2] # just the names for comments
                    s += str(diff).replace('e','d')+suff+'\n'
        # writes gas couples for dusty gas, data: wilke visco(2), LJ sigma and LJ e/K
        elif diffChoice==2:
            for ig in range(ng):
                for i2 in range(ig+1,ng):
                    diff = base['data'][ig][1]
                    if diff =='LJ': diff = 1e-7
                    suff = '  ;'+glist[ig]+'-'+glist[i2] # just the names for comments
                    s += str(diff).replace('e','d')+suff+'\n'
            # writes wilke viscosisty
            for ig in range(ng):
                visco = base['data'][ig][2]
                s += str(visco).replace('e','d')+'  ;'+glist[ig]+'\n'
            # writes lennard and jones parameters
            s += '\n\'lennard-jones\'\n'
            for ig in range(ng):
                LJ = base['data'][ig][3:]
                s += str(LJ[0])+'  '
                s += str(LJ[1])+'  ;'+glist[ig]+'\n'
        return s + '\n'
        
    def calcDij(self,data1,data2):
        """cal coeff diffusion in cm2/s for two species from chapman enskog theroy"""
        sig1,eK1,m1 = [float(a) for a in data1]
        sig2,eK2,m2 =  [float(a) for a in data2]
        sig_12 = (sig1+sig2)/2 # collision diameter in angstrom
        eps12_K = sqrt(eK1*eK2) # it is eps/kb for each species
        T = 293 # temperature K
        p = 1 # pressure atm
        kbT_e = T/eps12_K # kboltz*Temp/epsilon
        Omega = 0.9 + 0.5/kbT_e-log10(kbT_e)/4.8 # approcimation from table 2.1
        D12 = 0.00186*T**1.5*sqrt(1/m1+1/m2)/(p*sig_12**2*Omega)*1e-4
        return D12

    def writeDatabases(self):
        """writes the databases locally using only species used"""
        for name in ['comp','gases','complex','sorption','redox','mineral']:
            fil1 = self.core.gui.mainDir+os.sep+'utils'+os.sep+name+'.dbs'
            fil2 = self.fDir+os.sep+name+'.dbs'
            if name+'.dbs' not in os.listdir(self.fDir): # OA added 23/5/17
                os.system('copy '+fil1+' '+fil2)
#        for name in ['redox','mineral']: # 23/5 OA removed ,not usefull and print bad data
#            if len(self.chem.Base['MChemistry'][name]['rows'])>0:
#                f1 = open(self.fDir+os.sep+name+'.dbs','w')
#                f1.write(self.chem.writeDb(name))
#                f1.close()
                    
class min3pReader:
    
    def __init__(self,core,fDir, fName):
        self.fDir,self.fName = fDir,fName
        self.fullPath = fDir+os.sep+fName;#print self.fullPath
        self.data = {'gsp':None,'gst':None,'gsg':None,'gsv':None,'vel':None,'gsm':None,'gsc':None}
        self.struct = core.getValueFromName('Min3pFlow','P_Uns')==0
        self.names, self.core = {},core         
        
    def readOutput(self,ext,iper,varname):
        """reads an output file from Min3p, with three head lines
        and retrieves a matrix of correct shape (nvar,z,y,x), reads all tieme step and stores them
        gst : total conc, gsm : general (pH, Alk...), gsc : ion conc (incl complex)
        gsp : pressure.. wcontent..temp"""
        nper = self.core.nper 
        if self.data[ext] == None:
            if self.struct: # structured
                nx,ny,nz=self.core.getValueFromName('Min3pFlow','NX'),self.core.getValueFromName('Min3pFlow','NY'),self.core.getValueFromName('Min3pFlow','NZ')
                if ext=='vel': # velocity have diff shape
                    nx,ny,nz = max(nx-1,1),max(ny-1,1),max(nz-1,1)
                f1=open(self.fullPath+'_0.'+ext,'r');f1.readline()
                nam = f1.readline().replace('"','').replace(' ','').replace('fh_w','h_w').replace('\n','')
                f1.close()
                self.names[ext] = nam.split(',')[3:];#print nx,ny,nz,self.names
                for iper in range(nper-1):
                    mat = loadtxt(self.fullPath+'_'+str(iper+1)+'.'+ext,skiprows = 3) # OA 14/6 pb read _0
                    nr,nc = shape(mat);#print nr,nc
                    if iper==0: self.data[ext] = zeros((nper,nc-3,nz,ny,nx))
                    mat2 = zeros((nc-3,nz,ny,nx));#print ext,iper,shape(self.data[ext]),shape(mat2)
                    for iv in range(nc-3):
                        mat2[iv] = reshape(mat[:,iv+3],(nz,ny,nx))
                    self.data[ext][iper] = mat2
            else : # usntructured
                self.names[ext],self.data[ext] = self.readOutputVtk(ext,nper)
        ivar = self.names[ext].index(varname)
        return self.data[ext][iper,ivar]
        
    def readOutputVtk(self,ext,nper):
        '''reads a series of files defined by the extension name ext in vtk format'''
        names = []
        for iper in range(nper):
            f1=open(self.fullPath+'_'+str(iper)+'.'+ext+'.vtk')
            a = f1.read(); f1.close()
            b = a.split('SCALARS')[1:]
            nvar = len(b)
            if iper==0: 
                npts = len(b[0].split('\n'))-3
                data=zeros((nper,nvar,npts))
            n=-1
            for ivar in range(nvar):
                b1 = b[ivar].split('\n');#print b1[0]
                names.append(b1[0].split()[0])
                data[iper,ivar]=array(b1[2:n]).astype('float')
        return names,data
            
    def readHeadFile(self,core,iper):
        self.core = core
        hd = self.readOutput('gsp',iper,'h_w') # first variable is head
        return hd

    def readWcontent(self,core,iper):
        self.core = core
        wc = self.readOutput('gsp',iper,'s_w') # 3rd variable is theta_w
        return wc

    def readFloFile(self,core,iper):
        self.core = core
        vel = self.readOutput('vel',iper,rnge(2))
        nv,nz,ny,nx = shape(vel);#print 'redflo',shape(vel)
        # velocity have diff shape
        vel = concatenate((vel[:,:,:,:1],vel,vel[:,:,:,-1:]),axis=3);#print 'redflo',shape(vel)
        vel = vel[:,:,:,1:]/2+vel[:,:,:,:-1]/2;#print 'redflo',shape(vel)
        if ny>1:
            vel = concatenate((vel[:,:,:1,:],vel,vel[:,:,-1:,:]),axis=2);#print 'redflo',shape(vel)
            vel = vel[:,:,1:,:]/2+vel[:,:,:-1,:]/2;#print 'redflo',shape(vel)
        if nz>1:
            vel = concatenate((vel[:,:1,:,:],vel,vel[:,-1:,:,:]),axis=1);#print 'redflo',shape(vel)
            vel = vel[:,1:,:,:]/2+vel[:,:-1,:,:]/2;#print 'redflo',shape(vel)
        return vel

    def readUCN(self,core,option,iper,ispec,specname):
        self.core = core
        specname = specname.lower();#print option,iper,specname
        eng=core.getValueFromName('Min3pFlow','energy_balance');
        if eng==1:
            return self.readOutput('gsp',iper,'temp_n') 
        if specname =='ph':
            return self.readOutput('gsm',iper,'pH') 
        if specname =='pe':
            return self.readOutput('gsm',iper,2) #!!! wrong but i don't know where pe is
        if specname =='tracer':
            trc = self.readOutput('gsc',iper,'tracer'); #print shape(trc)
            return trc
        lgrp = ['comp','gases','mineral']
        suff = ['t','g','v'] # totals, gases, mineral vol fraction
        for i,grp in enumerate(lgrp): 
            slist = self.core.addin.chem.getListSpecies(grp)
            if specname in slist:
                cnc = self.readOutput('gs'+suff[i],iper,specname) # OA 24/5 name are now read directly
        return cnc
        
    def getPtObs(self,core,irow,icol,ilay,iper,option,ispec=0,specname=''):
        """get an observation point or a list of obs points for one period
        up to now, no different periods"""
        #print irow,icol,ilay,iper
        if core.addin.getDim() in ['Xsection','Radial']:
            ilay=irow*1;irow=[0]*len(ilay)
        specname = specname.lower()
        if option in ['Tracer','Chemistry']: # transport or chemistry
            a = self.readUCN(core,option,iper,ispec,specname);#print shape(a)
            obs = a[:,ilay,irow,icol]
        elif option == 'head':
            obs = self.readHeadFile(core,iper)[ilay,irow,icol]
        elif option == 'Wcontent':
            obs = self.readWcontent(core,iper)[ilay,irow,icol]
        elif option == 'flux': 
            obs = self.readFloFile(core,iper)[:,ilay,irow,icol]
        return obs
        
