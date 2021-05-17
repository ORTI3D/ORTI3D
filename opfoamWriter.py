# -*- coding: utf-8 -*-
"""
Created on Sun May 17 10:55:23 2020

@author: olivier
"""
# a file with point coordinates nbofpts ((x0 y0 z0) (x1 y1 z1)...)
# faces file nboffaces (nbpts (pt0 pt1 pt2..) nbpts (pt0 pt1...) ...)
# owner file  al ist of the cell that owns the considered face nb (0 1 2 0)
# neighbour 
from scipy import zeros,ones,array,arange,r_,c_,around,argsort,unique,cumsum,where,shape,\
    amin,amax,mod,ravel
from ilibq.geometry import *

import os

class opfoamWriter:
    def __init__(self,md,opfoam):
        '''md is the model, opf the openoam object'''
        self.mesh = md.addin.mesh
        self.md,self.opf = md,opfoam
            
    def writeFiles(self,fDir,dicBC,options,wriGeom=True):
        self.fDir, self.dicBC,self.options  = fDir, dicBC, options
        self.group = options['group']
        if 'system' not in os.listdir(fDir): os.mkdir(fDir+'system')
        if 'constant' not in os.listdir(fDir): os.mkdir(fDir+'constant')
        if '0' not in os.listdir(fDir): os.mkdir(fDir+'0')
        self.ttable = self.md.makeTtable()
        self.tlist = self.ttable['tlist']
        if 'polyMesh' not in os.listdir(fDir+'constant'): os.mkdir(fDir+'constant\\polyMesh')
        fDir1 = self.fDir +'constant\\polyMesh\\'
        if wriGeom:
            if self.md.getValueFromName('Modflow','MshType') == 0:
                points,faces,bfaces,fcup = self.opf.opfMeshReg(self.md)
            else :
                points,faces,bfaces,fcup = self.opf.opfMeshFromUsg(self.md);print('got geom from usg')  
        self.nlay = getNlayers(self.md)
        self.ncell_lay = self.opf.ncell_lay
        self.area,self.nbc = self.opf.area,self.opf.nbc
        if wriGeom:self.writeGeom(fDir1,points,faces,bfaces,fcup);print("geom written")
        self.writeCtrlDict()
        self.writeFvSchemes()
        self.writeFvSolutions()
        # write variables
        self.bcD0={'top':{'type':'zeroGradient'},'bottom':{'type':'zeroGradient'}}
        for i in range(self.nbc):
            self.bcD0['bc'+str(i)] = {'type':'zeroGradient'}
        self.writeConstantFields();print('constant field written')
        self.writeInitFields();print('init field written')
        if self.group=='Transport': self.writeTransport()
        if self.group=='Chemistry': self.writeChemistry()
        if 'options' not in os.listdir(fDir+'constant'): os.mkdir(fDir+'constant\\options')
        if 'sets' not in os.listdir(fDir1): os.mkdir(fDir1+'sets')
        self.writeOptions()
        
    def writeGeom(self,fDir,points,faces,bfaces,fcup):
        '''
        md is the model,fDir the major openfoam folder
        '''
        md,mesh = self.md,self.mesh
        nf0=shape(faces)[0];nfb=shape(bfaces)[0];nf=nf0+nfb
        nplay = shape(points)[0] # nb of points for each layer
        ncell = len(fcup) # nb of cells per layer
        nlay = getNlayers(md);self.nlay = nlay
        # below approximative, just for flat layers
        if md.addin.getDim()=='3D':
            if 'importArray' in md.dictype['Modflow']['disu.7']:
                zlist = self.getZfromPoints(points)
            else :
                zlist = [float(a) for a in md.addin.get3D()['topMedia']]
                zlist.append(float(md.addin.get3D()['zmin']))
        else : 
            zlist = [md.dicval['Modflow']['disu.7'][0],md.dicval['Modflow']['disu.8'][0]]
        zlist = zlist[-1::-1] # layers are reversed, start from 0 in OpF
        # points now same height for layer top
        sp = 'FoamFile \n{ version 2.0;\n format ascii;\n class vectorField;\n'
        sp += ' location  \"constant/polyMesh\";\n object points;\n}\n'
        sp += str(shape(points)[0]*(nlay+1))+'\n(\n'
        for z in zlist:
            if type(z)==type(0.): z = [z]*nplay
            for i in range(nplay):
                sp += '(%g %g %g)\n'%(points[i,0],points[i,1],z[i]) 
            #sp += '\n'
        sp += ')'       
        f1=open(fDir+os.sep+'points','w');f1.write(sp);f1.close()
        
        ############# write face list containing all faces
        sp = 'FoamFile \n{ version 2.0;\n format ascii;\n class faceList;\n'
        sp += ' location  \"constant/polyMesh\";\n object faces;\n}\n'
        sp += str(nf*nlay+ncell*(nlay+1))+'\n(\n'
        # internal faces for 1st nlay-1 # This part is quite slow!! because top face added there
        for il in range(nlay-1): 
            for i in range(ncell):
                fc = faces[faces[:,2]==i,:2].astype('int')
                for f in fc:
                    sp += '4('+str(f[0]+nplay*(il))+' '+str(f[1]+nplay*(il))+' '+str(f[1]+nplay*(il+1))+' '+str(f[0]+nplay*(il+1))+')\n'
                sp+=str(len(fcup[i])-1)+'('+' '.join([str(ip+nplay*(il+1)) for ip in fcup[i][:-1]])+ ')\n'
        # internal for top layer
        for i in range(ncell):
            fc = faces[faces[:,2]==i,:2].astype('int')
            for f in fc:
                sp += '4('+str(f[0]+nplay*(nlay-1))+' '+str(f[1]+nplay*(nlay-1))+' '+str(f[1]+nplay*(nlay))+' '+str(f[0]+nplay*(nlay))+')\n'
        print('internal faces written')
        # boundary faces
        for ib in range(self.nbc):
            for i in range(ncell):
                for il in range(nlay):
                    fc = bfaces[(bfaces[:,2]==i)&(bfaces[:,3]==-ib-1),:2].astype('int')
                    for f in fc:
                        sp += '4('+str(f[1]+nplay*(il+1))+' '+str(f[0]+nplay*(il+1))+' '+str(f[0]+nplay*(il))+' '+str(f[1]+nplay*(il))+')\n'
        print('bdy faces written')
        # top faces
        for i in range(ncell):
            sp+=str(len(fcup[i])-1)+'('+' '.join([str(ip+nplay*nlay) for ip in fcup[i][:-1]])+ ')\n'
        # bottom faces
        for i in range(ncell):
            sp+=str(len(fcup[i])-1)+'('+' '.join([str(ip) for ip in fcup[i][-1:0:-1]])+ ')\n'
        sp += '\n)'        
        f1=open(fDir+os.sep+'faces','w');f1.write(sp);f1.close()
        
        ######## owner (cell nb) contains all faces
        sp = 'FoamFile \n{ version 2.0;\n format ascii;\n class labelList;\n'
        sp += ' location  \"constant/polyMesh\";\n object owner;\n}\n'
        sp += str(nf*nlay+ncell*(nlay+1))+'\n(\n'
        # internal + top (internal)
        for il in range(nlay-1):
            for i in range(ncell):
                fc = faces[faces[:,2]==i,:2].astype('int')
                sp += '\n'.join([str(a+ncell*il) for a in [i]*(len(fc)+1)])+'\n'
        sp += '\n'.join([str(a+ncell*(nlay-1)) for a in sort(faces[:,2].astype('int'))])+'\n' # internal
        # boundaries
        for ib in range(self.nbc):
            for i in range(ncell):
                for il in range(nlay):
                    fc = bfaces[(bfaces[:,2]==i)&(bfaces[:,3]==-ib-1)].astype('int')
                    for f in fc: sp += str(i+ncell*(il))+'\n'
        # top and bottom
        sp += '\n'.join([str(f+ncell*(nlay-1)) for f in range(ncell)])+'\n' 
        sp += '\n'.join([str(f) for f in range(ncell)])+'\n' 
        sp += ')'
        f1=open(fDir+os.sep+'owner','w');f1.write(sp);f1.close()
        
        ########## neighbour (cell nb) contains only inside faces 
        sp = 'FoamFile \n{ version 2.0;\n format ascii;\n class labelList;\n'
        sp += ' location  \"constant/polyMesh\";\n object neighbour;\n}\n'
        sp += str(nf0*nlay+ncell*(nlay-1))+'\n(\n'
        # internal
        for il in range(nlay):
            for i in range(ncell):
                fc = faces[faces[:,2]==i,:].astype('int')
                if len(fc)>0: sp += '\n'.join([str(a+ncell*il) for a in fc[:,3]])+'\n'
                if il<nlay-1 : sp += str(i+ncell*(il+1))+'\n'
        sp += ')'
        f1=open(fDir+os.sep+'neighbour','w');f1.write(sp);f1.close()
        
        # boundary file
        sp = 'FoamFile \n{ version 2.0;\n format ascii;\n class polyBoundaryMesh;\n'
        sp += ' location  \"constant/polyMesh\";\n object boundary;\n}\n'
        sp += str(self.nbc+2)+'\n (\n';i0=nf0*nlay+ncell*(nlay-1)
        for i in range(self.nbc):
            sp += 'bc'+str(i)+' {type patch;'
            i1 = len(where(bfaces[:,3]==-i-1)[0])*nlay
            sp += 'nFaces '+str(i1)+';startFace '+str(i0)+';}\n'
            i0 += i1
        # top and bottom 
        sp += 'top {type patch;startFace '+str(i0)+';nFaces  '+str(ncell)+';}\n'
        i0 += ncell
        sp += 'bottom {type patch;startFace '+str(i0)+';nFaces  '+str(ncell)+';}\n'
        i0 += ncell
        #sp += 'defaultFaces {type empty;inGroups 1(empty);startFace '+str(i0)+';nFaces  '+str(ncell*2)+';}\n'
        f1=open(fDir+os.sep+'boundary','w');f1.write(sp+')');f1.close()
        ####### il faut faire renumberMesh apres ??? pas sur
        
    def getZfromPoints(self,points):
        '''
        get the z coords of the points grid, cannot be done through getValueLong or zblock
        as the points are not at the center of the cell
        this is calle donly when arrays are present (interpolated case not treated)
        as z will be reversed, here the list of points for each layer starts at the top
        '''
        md = self.md;lzout=[]
        md.lcellInterp = [] # to reset the values where to search (default cell centers)
        grd = md.addin.getFullGrid();intp,ysign,zdx,zdy=False,0,None,None
        fName0 = md.dicarray['Modflow']['disu.7'][0] #top 
        xx,yy,intp = points[:,0],points[:,1],False
        for i in range(self.nlay): 
            fNameExt = md.dicarray['Modflow']['disu.7'][i] #tops
            if fNameExt[-3:] == 'var' : ysign,zdx,zdy,zgrd = md.importGridVar(self.md.fileDir,fNameExt) # OA 13/6/20 add ysign
            else : zgrd = loadtxt(self.md.fileDir+fNameExt)
            if ysign == -1 : 
                zgrd = zgrd[-1::-1];zdy = zdy[-1::-1]
            lzout.append(linIntpFromGrid(md,grd,zgrd,xx,yy,intp,zdx,zdy))# removed [::-1]
        fNameExt = md.dicarray['Modflow']['disu.8'][-1] #bottom
        if fNameExt[-3:] == 'var' : ysign,zdx,zdy,zgrd = md.importGridVar(self.md.fileDir,fNameExt) # OA 13/6/20 add ysign
        else : zgrd = loadtxt(self.md.fileDir+fNameExt)
        if ysign == -1 : 
            zgrd = zgrd[-1::-1];zdy = zdy[-1::-1]
        lzout.append(linIntpFromGrid(md,grd,zgrd,xx,yy,intp,zdx,zdy))# removed [::-1]
        md.lcellInterp = [] # to reset the values where to search (default cell centers)
        return lzout
    
    def writeCtrlDict(self):
        ''' we still have simple definition of maxDeltaT'''
        tunit = self.md.dicval['Modflow']['dis.2'][4]
        tunit=4 # not done corectly on modlfow files!!
        if tunit== 5 : dt = 86400*365
        if tunit== 4 : dt = 86400
        if tunit== 3: dt = 3600
        if tunit== 2: dt = 60
        if tunit== 1 : dt = 1
        self.maxT = int(float(self.md.dicaddin['Time']['final'][0])*dt)
        intv = int(float(self.md.dicaddin['Time']['steps'][0])*dt)
        nstp = 10; #self.md.dicval['Modflow']['dis.8'][1]
        ctrlDict={'startTime':0,'endTime': int(self.maxT),'deltaT':str(int(intv/nstp/1000)),'writeInterval':str(intv),
                  'maxDeltaT':str(int(intv/nstp)),'maxCo':0.75,'dCmax':1e-2,'dCresidual':1e-3}
        solv = 'gwaterFoam'
        s = 'FoamFile\n{\n version 2.0;\n format ascii;\n class dictionary;\n'
        s += ' location \"system\";\n object controlDict;\n}\n'
        s += 'application '+solv+';\n'
        s += ' startFrom startTime;\n stopAt endTime;\n writeControl adjustableRunTime;\n'
        s += ' purgeWrite 0;\n writeFormat ascii;\n writePrecision 6;\n'
        s += ' writeCompression uncompressed;\n timeFormat fixed;\n timePrecision 0;\n'
        s += ' runTimeModifiable yes;\n\n'
        for k in ctrlDict.keys():
            s += k+' '+str(ctrlDict[k])+';\n'
        f1=open(self.fDir+os.sep+'system'+os.sep+'controlDict','w');f1.write(s);f1.close()

    def writeFvSchemes(self):
        schemeDict = {'ddtSchemes':{'default':'Euler'},
        'gradSchemes':{'default':'Gauss linear'},
        'divSchemes':{'default':'none','div(phit,Cw)':'Gauss vanLeer','div(phit,Cwi)':'Gauss vanLeer'}, #SuperBee
        'laplacianSchemes':{'default':'Gauss linear corrected'},
        'interpolationSchemes':{'default':'linear','M': 'vanLeer phiw'},
        'snGradSchemes':{'default':'corrected'},
        'fluxRequired':{'default':'no','h':' '}
        }
        s = 'FoamFile\n{\n version 2.0;\n format ascii;\n class dictionary;\n'
        s += ' location \"system\";\n object fvSchemes;\n}\n'
        for k in schemeDict.keys():
            s += '\n'+k+'\n{\n'
            for k1 in schemeDict[k].keys():
                s += '  '+k1+' '+str(schemeDict[k][k1])+';\n'
            s += '}'
        f1=open(self.fDir+os.sep+'system'+os.sep+'fvSchemes','w');f1.write(s);f1.close()
        
    def writeFvSolutions(self):# writing fvsolution
        s = 'FoamFile{version 2.0;format ascii;class  dictionary;location \"system\";object fvSolution;}'
        s += 'solvers \n{ \n'
        s += 'h{solver PCG;preconditioner  DIC;tolerance 1e-5;relTol 1e-8;}\n'
        s += '\"C.*\"{solver PBiCG;preconditioner  DILU;tolerance 1e-12;relTol 0;}\n}'
        s += '\nPicard {tolerance 0.01;maxIter 10;minIter 3;nIterStability  5;}'
        f1=open(self.fDir+'system\\fvSolution','w');f1.write(s);f1.close()
        
    def getVariable(self,modName,line):
        return ravel(self.md.getValueLong(modName,line,0)[-1::-1])
    
    def writeConstantFields(self):
        '''
        !!!now eps comes from mt3dms and K from modflow (and Kh/Kv)
        '''
        md = self.md
        # gravity field g
        s = 'FoamFile{version 2.0;format ascii;class uniformDimensionedVectorField;'
        s += 'location \"constant\"; object g;}\n'
        s += ' dimensions [0 1 -2 0 0 0 0];\n value  ( 0 0 -9.81 );'
        f1=open(self.fDir+'constant/g','w');f1.write(s+'\n}');f1.close()
        # porosity
        self.eps = self.getVariable('Mt3dms','btn.11')
        self.writeScalField('constant','eps',self.eps,self.bcD0)
        # permeability
        K = self.getVariable('Modflow','lpf.8')/86400/9.81e6;self.K=K
        self.writeScalField('constant','Kh',K,self.bcD0,dim='[0 2 0 0 0 0 0]')
        vrt = self.getVariable('Modflow','lpf.9')
        if self.options['Kvert'] == 'Kz':
            kv = vrt/86400/9.81e6
            #kv = maximum(kv,K/50) #♣ openfoam doe snot like too small ratio
        elif self.options['Kvert'] == 'Ratio':
            kv = K/vrt
        self.writeScalField('constant','Kv',kv,self.bcD0,dim='[0 2 0 0 0 0 0]')
        #if amax(K)==amin(K): kt = (K[0],K[0],kv[0]) ;#/khkv[0]
        #else : kt = c_[K,K,kv]#khkv]
        #self.writeTensorField('constant','K',kt,self.bcD0,dim='[0 2 0 0 0 0 0]')
        # thickness and bottom
        zb = md.Zblock
        thk = zb[:-1]-zb[1:];thk=thk[-1::-1];thk=ravel(thk);self.thk = thk  #-1::-1 for modflow
        self.writeScalField('constant','thk',thk,self.bcD0,dim='[0 1 0 0 0 0 0]')
        zbot = ravel(zb[1:][-1::-1])
        self.writeScalField('constant','zbot',zbot,self.bcD0,dim='[0 1 0 0 0 0 0]')
        # transport properties
        s = 'FoamFile{version 2.0;format ascii;class dictionary;location "constant";object transportProperties;}\n'
        if self.group == 'Flow' :s += 'activateTransport 0.;\n'
        elif self.group == 'Transport' :s += 'activateTransport 1.;\n'
        elif self.group == 'Chemistry' :s += 'activateTransport 3.;\n'
        uncf = '0'
        if md.dicaddin['Model']['type'] == 'Unconfined': uncf = '1'
        s += 'activateUnconfined '+uncf+';\n'
        std = '0'
        #☺if md.dicval['Modflow']['disu.9'][3] == 0: std = '1'
        s += 'flowStartSteady ' + std+';\n'
        s += 'phase.w{rho	rho [1 -3 0 0 0 0 0] 1e3;mu mu [1 -1 -1 0 0 0 0] 1e-3;}\n'
        if self.group in ['Chemistry','Transport']:
            s += 'alphaL alphaL [0 1 0 0 0 0 0] '+str(md.dicval['MfUsgTrans']['bct.6'][0])+';\n'
            s += 'alphaT alphaT [0 1 0 0 0 0 0] '+str(md.dicval['MfUsgTrans']['bct.7'][0])+';\n'
        f1=open(self.fDir+'constant\\transportProperties','w');f1.write(s);f1.close()
    
    def writeInitFields(self):
        # Uw
        self.writeVectField('0','Uw',[0,0,0],self.bcD0,'[0 1 -1 0 0 0 0]')
        # head h
        h0 = self.getVariable('Modflow','bas.5')
        self.writeScalField('0','h',h0,self.dicBC,dim='[0 1 0 0 0 0 0]')

    def writeTransport(self):
        C0 = self.getVariable('MfUsgTrans','bct.20')
        self.writeScalField('0','Cw',C0,self.bcD0,dim='[1 -3 0 0 0 0 0]')
        
    def writeChemistry(self):
        self.ncomp,self.lcomp,self.lspec = self.findSpecies(self.md)
        self.writePhreeqc(self.md,self.fDir)
        self.writePhqFoam(self.md,self.fDir)
        
    def writeOptions(self):
        '''
        '''
        # writing the fvOption file
        ncl = self.ncell_lay
        so ='FoamFile{version 2.0;format ascii;class dictionary;location \"constant\";object fvOptions;}\n'
        if 'bas.5' in self.md.diczone['Modflow'].dic.keys():
            so += 'hFix \n{type scalarmyFixedValueConstraint;\nactive true;\n'
            so += 'selectionMode cellSet;\ncellSet hfix;\nvolumeMode absolute;\n'
            so += 'fieldValues {h 1;}\n}\n'
            # writing the hfix file in cellSet (cells) and data in options
            lcell,ncell,zcell = self.writeOptionCells('bas.5','hfix')
            hmat = self.ttable['bas.5'];nr,nz = shape(hmat)
            mult = []
            for iz in range(len(lcell)) : 
                mult.append([1]*ncell[iz])
            self.writeOptionData(hmat,lcell,ncell,'hfix',mult=mult)
        if 'wel.1' in self.md.diczone['Modflow'].dic.keys():
            so += 'hSouwel \n{type scalarmySemiImplicitSource;\nactive true;\n'
            so += 'selectionMode cellSet;\ncellSet hwel;\nvolumeMode absolute;\n'
            so += 'injectionRateSuSp {kr (1 0);h (1 0);}\n}\n'
            # writing the hWel file in cellSet (cells) and data in options
            lcell,ncell,zcell = self.writeOptionCells('wel.1','hwel')
            qmat = self.ttable['wel.1'];nr,nz = shape(qmat)
            mult = [];eps=0.25
            for iz in range(nz) : 
                kr = self.getPermScaled(lcell[iz]);#print(kr)#/self.thk[lcell[iz]]
                area = self.area[mod(lcell[iz],ncl)]
                mult.append(kr/area/86400/self.thk[lcell[iz]])
            self.writeOptionData(qmat,lcell,ncell,'hwel',mult=mult)
            
        if ('rch.2' in self.md.diczone['Modflow'].dic.keys()) or (self.md.dicval['Modflow']['rch.2'][0] != 0):
            so += 'hSourch \n{type scalarmySemiImplicitSource;\nactive true;\n'
            so += 'selectionMode cellSet;\ncellSet hrch;\nvolumeMode absolute;\n'
            so += 'injectionRateSuSp {sw (1 0);h (1 0);}\n}\n'
            lcell,ncell,zcell = self.writeOptionCells('rch.2','hrch')
            nz = len(lcell)
            rmat = self.ttable['rch.2'].astype(float);nr,nc = shape(rmat)
            vbase = self.md.dicval['Modflow']['rch.2'][0]
            rmat  = c_[rmat,ones((nr,1))*vbase]# for the background
            mult = []
            for iz in range(len(lcell)):
                mult.append(ones(ncell[iz])/86400)# thickness will be calc in fvOptions
            self.writeOptionData(rmat,lcell,ncell,'hrch',mult=mult)
            # write thk in constant/options
            sthk ='\n'.join(['%g'%x for x in self.thk])
            f1=open(self.fDir+'constant\\options\\thk','w');f1.write(sthk);f1.close()

        for n in ['riv','drn','ghb']:
            if n+'.1' in self.ttable.keys():
                so += 'hSou'+n+' \n{type scalarmySemiImplicitSource;\nactive true;\n'
                so += 'selectionMode cellSet;\ncellSet h'+n+';\nvolumeMode absolute;\n'
                so += 'injectionRateSuSp {sw (1 0);h (1 1);}\n}\n'
                lcell,ncell,zcell = self.writeOptionCells(n+'.1','h'+n)
                mat = self.ttable[n+'.1'];nr,nz = shape(mat);mat = ones((nr,nz)) # mat sueless here, put 1
                mult, mult2 = [],[]
                dicz = self.md.diczone['Modflow'].dic[n+'.1']
                for iz in range(len(lcell)) : 
                    if len(lcell[iz])==0:
                        mult.append([]);mult2.append([]);continue
                    v1 = self.area[mod(lcell[iz],ncl)]*self.thk[lcell[iz]]
                    a = dicz['value'][iz].split('$')[1].split('\n')
                    zz,hcond = float(a[0]),float(a[1])
                    if zcell[iz]==0: zcell[iz]=[zz]*ncell[iz]
                    mult.append(-hcond*array(zcell[iz])/v1/86400)
                    mult2.append(hcond/v1/86400)
                self.writeOptionData(mat,lcell,ncell,'h'+n,mult=mult,mult2=mult2)
        if self.group == 'Chemistry':
            ## need to know if ph.4 is at a flow bc or at a well
            if 'c' not in self.dicBC.keys():
                if 'ph.4' in self.ttable.keys():
                    for i in range(self.ncomp):
                        so += 'cfix'+str(i)+' \n{type scalarmyFixedValueConstraint;\nactive true;\n'
                        so += 'selectionMode cellSet;\ncellSet cfix;\nvolumeMode absolute;\n'
                        so += 'fieldValues {Cw'+str(i)+' 1;}\n}\n'
                    lcell,ncell,zcell = self.writeOptionCells('ph.4','cfix')
                    mat = self.ttable['ph.4']
                    self.writeOptionData(mat,lcell,ncell,'cfix',formt='int')
                
                    
        f1=open(self.fDir+'constant\\fvOptions','w');f1.write(so);f1.close()

    def writeOptionCells(self,line,fname):
        '''
        writes the cell numbers in constant/polyMesh/sets, lcell is the list of cells
        for each considered zone, ncell the nb of these cells
        '''
        modName = 'Modflow'
        if line[:2]=='ph': modName = 'Pht3d'
        dicz = self.md.diczone[modName].dic[line]
        nclay,nlay = self.ncell_lay,self.nlay
        if self.md.getValueFromName('Modflow','MshType')==0:
            lcell,ncell,zcell = [],[],[]
            nx,ny,xv,yv = getXYvects(self.md)
            dx,dy = xv[1:]-xv[:-1],yv[1:]-yv[:-1]
            dxm,dym = meshgrid(dx,dy)
            self.area = ravel(dxm*dym)
            for media in range(self.nlay):
                m = zone2grid(self.md,modName,line,media)[-1::-1,:] # modflow inversion
                m1 = zone2grid(self.md,modName,line,media,opt='zon')[-1::-1,:] # modflow inversion
                lzon = unique(m1).nonzero()[0]
                for iz in lzon:
                    if iz>len(lcell): 
                        lcell.append([]);zcell.append([]);ncell.append(0)
                    idx = where(m1==iz)
                    lcell[iz-1].extend(list(idx[0]*nx+idx[1]))
                    zcell[iz-1].extend(list(m[idx[0],idx[1]]))
                    ncell[iz-1] += len(idx[0])
        else :
            lcell,ncell,zcell = [],[],[]
            for iz in range(len(dicz['name'])):
                lmedia = dicz['media'][iz]
                if type(lmedia) != type([5,6]) : lmedia = [lmedia]
                lc0,zval = zmesh(self.md,dicz,lmedia[0],iz)
                if len(lc0)==1: lc0 = lc0[0]
                nc,lc1 = 0,[]
                for im in lmedia:
                    lc1.extend(lc0+nclay*(self.nlay-im-1));nc += 1  # nlay-im because inversion vs modflow
                lcell.append(lc1);zcell.append(zval);ncell.append(len(lc0)*nc)
        nlay = self.nlay
        if fname=='hrch': # to add the base value
            lc_toplay = set(list(range(nclay*(nlay-1),nclay*nlay)))
            lctot = []
            for lc in lcell: lctot.extend(lc)
            lc_toplay.difference_update(lctot)
            lcell.append(list(lc_toplay))
            ncell.append(len(lc_toplay))
            zcell.append(0)
        s ='FoamFile{version 2.0;format ascii;class cellSet;location \"constant/polyMest/sets\";object hWel;}\n'
        s += '\n'+str(sum(ncell))+'\n(\n'
        s += '\n'.join([str(a) for lc1 in lcell for a in lc1])+'\n)'
        f1=open(self.fDir+'constant\\polyMesh\\sets\\'+fname,'w');f1.write(s);f1.close()
        return lcell,ncell,zcell
    
    def writeOptionData(self,mat,lcell,ncell,fname,mult=None,mult2=None,formt='float'): 
        '''
        write the option semiImplcit for one variable, mult is a list of vector (for well)
        '''
        nr,a = shape(mat);nzo=len(lcell);dt = 86400
        if mult==None: mult = [[1]*nc for nc in ncell]
        if mult2==None: mult2 = [[0]*nc for nc in ncell]
        if formt == 'float': fmt0 = '0 %9i %e %e\n';fmt1='%9i %9i %e %e\n'
        elif formt == 'int': fmt0 = '0 %9i %9i %9i\n';fmt1='%9i %9i %9i %9i\n'
        s = '%5i %8i \n'%(self.nlay,self.ncell_lay)
        val = zeros((sum(ncell),3))
        sc1 = cumsum(ncell);sc0 = r_[0,sc1[:-1]]
        for iz in range(nzo):  # 1st time step
            if len(lcell[iz])==0: continue
            val[sc0[iz]:sc1[iz],0] = lcell[iz]
            val[sc0[iz]:sc1[iz],1] = float(mat[0,iz])*array(mult[iz])
            val[sc0[iz]:sc1[iz],2] = float(mat[0,iz])*array(mult2[iz])
        for v in val:
            s += fmt0%(int(v[0]),v[1],v[2])
        for it in range(1,nr):
            if all(mat[it] == mat[it-1]): continue
            for iz in range(nzo):
                val[sc0[iz]:sc1[iz],1] = float(mat[it,iz])*array(mult[iz])
                val[sc0[iz]:sc1[iz],2] = float(mat[it,iz])*array(mult2[iz])
            for v in val:
                s += fmt1%(self.tlist[it]*dt,int(v[0]),v[1],v[2])
        for v in val:
            s += '%9i %9i 0 0\n'%(self.tlist[-1]*dt+10,int(v[0]))
        f1=open(self.fDir+'constant\\options\\'+fname,'w');f1.write(s);f1.close()
        
    def getPermScaled(self,lcell):
        """return the permeability for a list of layer, col rows scaled by the
        sum of permeability for this list"""
        ka=ones(len(lcell))*0.;#print 'mfi permsc',shape(ka),shape(K),ilay,irow,icol
        for i in range(len(lcell)):
            ka[i] = self.K[lcell[i]]*self.thk[lcell[i]]
        return ka/sum(ka)

    def writeScalField(self,loc,name,values,dicBC,dim='[0 0 0 0 0 0 0]'):
        '''
        values is a list of values of the scalar
        '''
        s ='FoamFile{version 2.0;format ascii;class volScalarField;location '+loc+';object '
        s += name+';}\n\n'
        s += 'dimensions '+dim+';\n\n'
        if amax(values)==amin(values):
            s += 'internalField uniform '+str(amin(values))+';\n'
        else :
            s += 'internalField nonuniform List<scalar>\n'+str(len(values))+'('
            s += '\n'.join(['{:.3e}'.format(a) for a in values])+');\n'
        s += '\nboundaryField\n{'
        bc1 = self.bcD0.copy()
        if name in dicBC.keys():
            for k in self.bcD0.keys():
                if k in dicBC[name] : bc1[k] = dicBC[name][k]
        for k in bc1.keys():
            s += '\n  '+k+'\n  {\n'
            for k1 in bc1[k].keys():
                s += '    '+k1+' '+str(bc1[k][k1])+';\n'
            s += '  }'
        f1=open(self.fDir+os.sep+loc+'/'+name,'w');f1.write(s+'\n}');f1.close()
        
    def writeVectField(self,loc,name,values,bcDict,dim='[0 0 0 0 0 0 0]'):
        '''
        values is a list of three values
        '''
        s ='FoamFile{version 2.0;format ascii;class volVectorField;location '+loc+';object '
        s += name+';}\n\n'
        s += 'dimensions '+dim+';\n\n'
        s += 'internalField uniform ('+' '.join(['%g'%a for a in values])+');\n'
        s += '\nboundaryField\n{'
        for k in bcDict.keys():
            s += '\n  '+k+'\n  {\n'
            for k1 in bcDict[k].keys():
                s += '    '+k1+' '+str(bcDict[k][k1])+';\n'
            s += '  }'
        f1=open(self.fDir+os.sep+loc+'/'+name,'w');f1.write(s+'\n}');f1.close()

    def writeTensorField(self,loc,name,values,bcDict,dim='[0 0 0 0 0 0 0]'):
        '''
        values is an array 3 colum that are transofrmed in a symmetric tensor
        '''
        ncell = len(values)
        s ='FoamFile{version 2.0;format ascii;class volTensorField;location '+loc+';object '
        s += name+';}\n\n'
        s += 'dimensions '+dim+';\n\n'
        if len(values)==3:
            s += 'internalField uniform '
            v1,v3 = values[0],values[2]
            s +='('+str(v1)+' 0 0 0 '+str(v1)+' 0 0 0 '+str(v3)+');\n'
        else : #vairable values in space
            s += 'internalField nonuniform List<tensor>\n'+str(ncell)+'\n('
            for ic in range(ncell):
                v1,v3 = values[ic,0],values[ic,2]
                s +='(%g 0 0 0 %g 0 0 0 %g)\n'%(v1,v1,v3)
            s += ');\n'
        s += '\nboundaryField\n{'
        for k in bcDict.keys():
            s += '\n  '+k+'\n  {\n'
            for k1 in bcDict[k].keys():
                s += '    '+k1+' '+str(bcDict[k][k1])+';\n'
            s += '  }'
        f1=open(self.fDir+os.sep+loc+'/'+name,'w');f1.write(s+'\n}');f1.close()
                
    def findSpecies(self,md):
        '''find the components and species for phreeqc'''
        listE = md.addin.pht3d.getDictSpecies()
        listS = listE['i'];listS.extend(listE['k']);listS.extend(listE['kim'])
        listS.sort()
        nbs,ncomp,lcomp,lspec,i = len(listS),0,[],[],0
        while i<nbs:
            if listE['i'][i] == 'O(0)': lspec.append('O(0)')
            elif '(' in listE['i'][i]: 
                lspec.append(listE['i'][i])
                lcomp.append(listE['i'][i]);ncomp +=1
                if '(' in listE['i'][i+1]: 
                    lspec.append(listE['i'][i+1]);i +=1
            else :
                lcomp.append(listE['i'][i]);ncomp += 1
            i += 1
        ncomp += 2 # there is pH and pe but we need to have H20,H,0,charge and not pH,pe
        if len(lspec)==0:
            lspec = [listE['i'][0]] # we need at least one species for sel out
        return ncomp,lcomp,lspec
    
    def writePhreeqc(self,md,fDir):
        """this routine writes a phreeqc file where all solutions are written in
        phreqc format to be able to test their equilibrium before running pht3d
        1. tabke background, then cycle through pht3d zones
        2. get the solution number, phase number...it does not take rates
        3 write them in phreeqc format"""
        s = 'Database '+fDir+'\pht3d_datab.dat \n'
        s += 'Selected_output \n  -totals '+' '.join(self.lspec)+'\n\n'
        listE = md.addin.pht3d.getDictSpecies()
        chem = md.addin.pht3d.Base['Chemistry']
        ncell = self.ncell_lay
        solu = chem['Solutions'];
        dicz = md.diczone['Pht3d']
        if 'ph.4' in dicz.dic.keys():
            nbsol = max([int(a) for a in dicz.dic['ph.4']['value']])+1
        else :
            nbsol = 1
        listS = listE['i'];listS.extend(listE['k']);listS.extend(listE['kim'])
        listS.sort()
        for isol in range(nbsol):
            s += '\nSolution '+str(1-isol)+' \n units mol/L \n'
            for esp in listS: # go through phase list
                ie=solu['rows'].index(esp);#print esp,phases['rows'],ip,phases['data'][ip] # index of the phase
                conc=solu['data'][ie][isol+1] #backgr concentration of species
                s += esp+' '+str(conc)+'\n'
        pht = md.getValueLong('Pht3d','ph.3',0) #this is a vector ncell
        dInd={'Solutions':pht/1000.,
              'Phases':mod(pht,1000)/100,'Gases':mod(pht,1000)/100,
              'Exchange':mod(pht,100)/10,
              'Surface':mod(pht,10)}
        # mineral phases
        phases=chem['Phases'];
        nbph = len(unique(dInd['Phases']))
        for ip in range(nbph):
            if len(listE['p'])>0 : s += '\nEquilibrium_Phases '+str(ip)+'\n'
            for esp in listE['p']: # go through phase list
                ie = phases['rows'].index(esp);#print esp,phases['rows'],ip,phases['data'][ip] # index of the phase
                IS,conc = phases['data'][ie][ip+1:ip+3] #backgr SI and concentration of phase
                s += esp+' '+str(IS)+' '+str(float(conc)/unique(self.eps)[0])+'\n'
        # exchanger
        exc=chem['Exchange'];
        nbex = len(unique(dInd['Exchange']))
        for ix in range(nbex):
            if len(listE['e'])>0 :s += '\nExchange '+str(ix+1)+'\n'
            for esp in listE['e']: # go through phase list
                ie = exc['rows'].index(esp);#print esp,phases['rows'],ip,phases['data'][ip] # index of the phase
                conc = exc['data'][ie][ix+1] #backgr SI and concentration of phase
                s += esp+' '+str(float(conc)/unique(self.eps)[0])+'\n'
            if len(listE['e'])>0 :s += '-equilibrate with solution 1\n'
        # surface
        surf=chem['Surface'];
        nbsu = len(unique(dInd['Surface']))
        for isu in range(nbsu):
            if len(listE['s'])>0 :s += '\nSurface '+str(isu+1)+'\n'
            for esp in listE['s']: # go through phase list
                ie = surf['rows'].index(esp);#print esp,phases['rows'],ip,phases['data'][ip] # index of the phase
                conc = surf['data'][ie][isu+1] #backgr SI and concentration of phase
                s += esp+' '+str(float(conc)/unique(self.eps)[0])+' '
                s += ' '.join(surf['data'][ie][6:8])+' \n'
            if len(listE['s'])>0 :s += '-equilibrate with solution 1\n-no_edl\n'
        # kinetics (up to now just one)
        rates = chem['Rates'];
        parmk = md.addin.pht3d.calcNbParm()
        lkin = listE['k'];lkin.extend(listE['kim'])
        if len(lkin)>0 :s += '\nKinetics 1-'+str(ncell)+'\n'
        for nom in lkin:
            iek = rates['rows'].index(nom);
            s += '\n'+nom+'\n -parms '  #param k
            for ip in range(parmk[nom]):
                s += str(rates['data'][iek][ip+2])+' '
            s += '\n-formula '+rates['data'][iek][-1] +'\n' # formula
        rates = chem['Kinetic_Minerals'];
        for nom in listE['kp']:
            iek = rates['rows'].index(nom);
            iep = phases['rows'].index(nom);
            s += '\n'+nom+'\n -m0 '+ str(phases['data'][iep][2])+'\n'  
            s += ' -parms '
            for ip in range(parmk[nom]): s += str(rates['data'][iek][ip+1])+' '
        f1 = open(fDir+os.sep+'initChem.pqi','w')
        f1.write(s);f1.close()

    def writePhqFoam(self,md,fDir):
        '''phqfoam contains the nb of cell, nb of mobile components and species
        below the lines are 
        Solution,Equilibrium phases,Exchange,Surface,Gas phase,Solid solutions,Kinetics
        dInd contain vectors of solu/assembl nb in the domain, as a dict
        '''
        listE = md.addin.pht3d.getDictSpecies();
        pht = md.getValueLong('Pht3d','ph.3',0)[-1::-1] # -1::-1 as layer are inverted
        dInd={'Solutions':pht/1000.,
              'Phases':mod(pht,1000)/100,
              'Gases':mod(pht,1000)/100,
              'Exchange':mod(pht,100)/10,
              'Surface':mod(pht,10)}
        shortn=['k','i','kim','g','p','e','s','kp']
        longn=['Solutions','Solutions','Solutions','Gases','Phases','Exchange','Surface',
               'Phases']
        s = str(self.ncell_lay*self.nlay)+' '+str(self.ncomp)+' 1 \n' # solid-units=1 here
        for sp in ['i','p','e','s','g']:
            indx = shortn.index(sp)
            val = ravel(dInd[longn[indx]])
            if len(listE[sp])>0: # there is at least one speces in this group
                if len(unique(val))==1:  # just one solution or phase...
                    s += '0  1 \n'  
                else :
                    s += '-1 '+' '.join([str(int(a)) for a in val]) + '\n'
            else:
                s += '0  -1 \n'
        s += '0  -1 \n' # solid solution not used
        flgK = 0
        for sp in ['k','kim','kp']: # special as k and kp are in the same category
            if len(listE[sp]) > 0: flgK = 1
        if flgK : s += '0  1 \n'
        else : s += '0  -1 \n'
        f1 = open(fDir+os.sep+'phqfoam.txt','w')
        f1.write(s);f1.close()

class opfoamReader:
    def __init__(self,md,opfoam):
        '''md is the model, opf the openfoam object'''
        self.ncell_lay = opfoam.ncell_lay
        self.md,self.opf = md,opfoam
        self.ttable = md.makeTtable()
        self.tlist = self.ttable['tlist']
        self.nlay = getNlayers(md)
        self.ncell = self.ncell_lay*self.nlay
        
    def readScalar(self,name):
        ''' returns a scalar for all the times reshapes in layers'''
        def rd1(dt,i,n):
            f1=open(self.md.fileDir+os.sep+str(i*dt)+os.sep+n,'r')
            s = f1.read();f1.close()
            s1 = s.split('<scalar>')[1].split('(')[1].split(')')[0].split('\n')
            data = np.array(s1[1:-1]).astype('float')
            return data
        dt = int((self.tlist[1]-self.tlist[0])*86400);
        ntime=len(self.tlist)
        data = np.zeros((ntime,self.ncell))
        for i in range(1,ntime):
            data[i,:]= rd1(dt,i,name)
        V = reshape(data,(ntime,self.nlay,self.ncell_lay))
        return V


