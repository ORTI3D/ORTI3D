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
    def __init__(self,core,opfoam,orientation):
        '''core is the model, opf the openoam object'''
        self.mesh = core.addin.mesh
        self.core,self.opf = core,opfoam
        self.orientation = orientation # 'z' or 'y' for xsect or 'r' for radial
        
    def writeFiles(self,fDir,dicBC,options,wriGeom=True):
        self.fDir, self.dicBC,self.options  = fDir, dicBC, options
        self.group = options['group']
        if 'system' not in os.listdir(fDir): os.mkdir(fDir+'system')
        if 'constant' not in os.listdir(fDir): os.mkdir(fDir+'constant')
        if '0' not in os.listdir(fDir): os.mkdir(fDir+'0')
        self.ttable = self.core.makeTtable()
        self.tlist = self.ttable['tlist']
        self.nlay = getNlayers(self.core)
        self.ncell_lay = self.opf.ncell_lay
        self.ncell = self.nlay*self.ncell_lay
        self.carea,self.nbc = self.opf.carea,self.opf.nbc
        if 'polyMesh' not in os.listdir(fDir+'constant'): os.mkdir(fDir+'constant\\polyMesh')
        fDir1 = self.fDir +'constant\\polyMesh\\'
        self.MshType = self.core.getValueFromName('OpenFlow','MshType')
        if wriGeom:
            if self.MshType == 0:
                points,faces,bfaces,fcup = self.opf.opfRect()
            else :
                points,faces,bfaces,fcup = self.opf.opfMesh2Faces();print('got geom from usg')  
            self.writeGeom(fDir1,points,faces,bfaces,fcup);print("geom written")
        self.writeCtrlDict()
        self.writeFvSchemes()
        self.writeFvSolutions()
        # write variables
        self.bcD0={'top':{'type':'zeroGradient'},'bottom':{'type':'zeroGradient'}}
        for i in range(self.nbc):
            self.bcD0['bc'+str(i)] = {'type':'zeroGradient'}
        self.writeConstantFields();print('constant field written')
        self.writeInitFields();print('init field written')
        if 'options' not in os.listdir(fDir+'constant'): os.mkdir(fDir+'constant\\options')
        if 'sets' not in os.listdir(fDir1): os.mkdir(fDir1+'sets')
        if self.group in ['Chemistry']: #'Transport',
            self.ncomp,self.lcomp,self.lspec = self.findSpecies(self.core)
        self.writeOptions()
        if self.group=='Transport': 
            self.writeTransport()
        if self.group=='Chemistry': 
            self.writeChemistry()
        
    def writeGeom(self,fDir,points,faces,bfaces,fcup):
        '''
        core is the model,fDir the major openfoam folder
        '''
        core,mesh = self.core,self.mesh
        nplay = shape(points)[0] # nb of points for each layer
        ncell = len(fcup) # nb of cells per layer
        nlay = getNlayers(core);self.nlay = nlay
        # below approximative, just for flat layers
        if core.addin.getDim()=='3D':
            if 'importArray' in core.dictype['OpenFlow']['dis.6']:
                zlist = self.getZfromPoints(points)
            else :
                zlist = [float(a) for a in core.addin.get3D()['topMedia']]
                zlist.append(float(core.addin.get3D()['zmin']))
            zlist = zlist[-1::-1] # layers are reversed, start from 0 in OpF
        elif self.orientation == 'r': # radial case
            zlist = [-points[:,0]/100,points[:,0]/100] # z_radial = x/100
            zlist[0][points[:,0]==0]= -points[:,0][1]/100
            zlist[1][points[:,0]==0]= points[:,0][1]/100
            #bfaces = bfaces[bfaces[:,-1]==-2] # the last face is the one at the well
        else : 
            zlist = [core.dicval['OpenFlow']['dis.7'][0],core.dicval['OpenFlow']['dis.6'][0]]
        nf0=shape(faces)[0];nfb=shape(bfaces)[0];nf=nf0+nfb
        # points now same height for layer top
        sp = 'FoamFile \n{ version 2.0;\n format ascii;\n class vectorField;\n'
        sp += ' location  \"constant/polyMesh\";\n object points;\n}\n'
        nbp = shape(points)[0]*(nlay+1)
        sp += str(nbp)+'\n(\n'
        x0,y0 = amin(points[:,0]),amin(points[:,1])
        strt = 0
        for i,z in enumerate(zlist):
            sp+='('
            if type(z)==type(0.): z = [z]*nplay
            coo = zeros((shape(points)[0],3))
            coo[:,0],coo[:,1],coo[:,2] = points[:,0]-x0,points[:,1]-y0,z
            sp += ')\n('.join([' '.join(x) for x in coo.astype('str')]) 
            sp += ')\n'
        sp += ')'       
        f1=open(fDir+os.sep+'points','w');f1.write(sp);f1.close()
        
        ############# write face list containing all faces
        sp = 'FoamFile \n{ version 2.0;\n format ascii;\n class faceList;\n'
        sp += ' location  \"constant/polyMesh\";\n object faces;\n}\n'
        nbf = nf*nlay+ncell*(nlay+1)
        sp += str(nbf)+'\n(\n'
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
        nbo = nf*nlay+ncell*(nlay+1)
        sp += str(nbo)+'\n(\n'
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
        
        ########################### boundary file
        sp = 'FoamFile \n{ version 2.0;\n format ascii;\n class polyBoundaryMesh;\n'
        sp += ' location  \"constant/polyMesh\";\n object boundary;\n}\n'
        i0=nf0*nlay+ncell*(nlay-1);sp1='';nbb=0
        for i in range(self.nbc):
            i1 = len(where(bfaces[:,3]==-i-1)[0])*nlay
            if i1==0:continue
            nbb+=1
            sp1 += 'bc'+str(i)+' {type patch;'
            sp1 += 'nFaces '+str(i1)+';startFace '+str(i0)+';}\n'
            i0 += i1
        sp += str(nbb+2)+'\n (\n'+sp1;
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
        this is called only when arrays are present (interpolated case not treated)
        as z will be reversed, here the list of points for each layer starts at the top
        '''
        core = self.core;lzout=[];zb=core.Zblock;dzmin=(amax(zb)-amin(zb))/500
        core.lcellInterp = [] # to reset the values where to search (default cell centers)
        grd = core.addin.getFullGrid();intp,ysign,zdx,zdy=False,0,None,None
        fName0 = core.dicarray['OpenFlow']['dis.6'][0] #top 
        xx,yy,intp,z0 = points[:,0],points[:,1],False,1e6
        for i in range(self.nlay): 
            fNameExt = core.dicarray['OpenFlow']['dis.6'][i] #tops
            if fNameExt[-3:] == 'var' : ysign,zdx,zdy,zgrd = core.importGridVar(self.core.fileDir,fNameExt) # OA 13/6/20 add ysign
            else : zgrd = loadtxt(self.core.fileDir+fNameExt)
            if ysign == -1 : 
                zgrd = zgrd[-1::-1];zdy = zdy[-1::-1]
            z1 = linIntpFromGrid(core,grd,zgrd,xx,yy,intp,zdx,zdy)
            lzout.append(minimum(z1,z0-dzmin)) # to avoid negative thickness
            z0 = z1*1
        fNameExt = core.dicarray['OpenFlow']['dis.7'][-1] #bottom
        if fNameExt[-3:] == 'var' : ysign,zdx,zdy,zgrd = core.importGridVar(self.core.fileDir,fNameExt) # OA 13/6/20 add ysign
        else : zgrd = loadtxt(self.core.fileDir+fNameExt)
        if ysign == -1 : 
            zgrd = zgrd[-1::-1];zdy = zdy[-1::-1]
        z1 =linIntpFromGrid(core,grd,zgrd,xx,yy,intp,zdx,zdy)
        lzout.append(minimum(z1,z0-dzmin))# removed [::-1]
        core.lcellInterp = [] # to reset the values where to search (default cell centers)
        return lzout
    
    def writeCtrlDict(self):
        ''' we still have simple definition of maxDeltaT'''
        tunit = self.core.dicval['Modflow']['dis.2'][4]
        tunit=4 # not done corectly on modlfow files!!
        if tunit== 5 : dt = 86400*365
        if tunit== 4 : dt = 86400
        if tunit== 3: dt = 3600
        if tunit== 2: dt = 60
        if tunit== 1 : dt = 1
        self.maxT = int(float(self.core.dicaddin['Time']['final'][0])*dt)
        self.ttable = self.core.makeTtable()
        self.tlist = self.ttable['tlist']
        intv = int(float(self.core.dicaddin['Time']['steps'][0])*dt)
        nstp = 10; #self.core.dicval['Modflow']['dis.8'][1]
        ctrlDict={'startTime':0,'endTime': int(self.maxT),'deltaT':str(int(intv/nstp/100)),
                  'maxDeltaT':str(int(intv/nstp)),'maxCo':5,'dCmax':1e-2,'dCresidual':1e-3}
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
        # write the function for variable times
        s += '#include "$FOAM_CASE/system/writeInterval"\n'
        s += 'functions \n { \n    fileUpdate1\n    {\n'
        s += '     type timeActivatedFileUpdate;\n	 libs ("libutilityFunctionObjects.so");\n'
        s += '     writeControl timeStep;\n     writeInterval 1;\n'
        s += '     fileToUpdate  "$FOAM_CASE/system/writeInterval";\n'
        s += '     timeVsFile \n    (\n'
        t_old,dt_old = -1,-1
        for i,t in enumerate(self.tlist):
            dt = t-t_old
            #if dt != dt_old:
            s1 = str(int(t_old*86400))
            s += '    ('+s1+'  "$FOAM_CASE/system/writeInterval.%0*i")\n'%(2,i)
            f1=open(self.fDir+os.sep+'system'+os.sep+'writeInterval.%0*i'%(2,i),'w')
            s1 = str(int(t*86400))  #dt*86400
            f1.write('writeInterval '+s1+';')
            f1.close()
            t_old,dt_old = t*1,dt*1
        s += '    );\n    }\n }'
        f1=open(self.fDir+os.sep+'system'+os.sep+'controlDict','w');f1.write(s);f1.close()
        f1=open(self.fDir+os.sep+'system'+os.sep+'writeInterval','w')
        f1.write('writeInterval 100;')
        f1.close()

    def writeFvSchemes(self):
        schemeDict = {'ddtSchemes':{'default':'Euler'},
        'gradSchemes':{'default':'Gauss linear'},
        'divSchemes':{'default':'none','div(phiw,Cw)':'Gauss vanLeer','div(phiw,Cwi)':'Gauss vanLeer',
                     'div(phig,Cg)':'Gauss vanLeer','div(phig,Cgi)':'Gauss vanLeer' }, #vanLeer,SuperBee
        'laplacianSchemes':{'default':'Gauss linear corrected'},
        'interpolationSchemes':{'default':'linear','M': 'vanLeer phiw','krg':'upwind phig','krw':'upwind phiw'},
        'snGradSchemes':{'default':'corrected'},
        'fluxRequired':{'default':'no','h':' ','p':' '}
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
        s += 'h{solver PBiCGStab;preconditioner  DIC;tolerance 1e-12;relTol 0;}\n'
        s += 'hp{solver PBiCGStab;preconditioner  DIC;tolerance 1e-12;relTol 0;}\n'
        s += 'p{solver PBiCGStab;preconditioner  DIC;tolerance 1e-12;relTol 0;}\n'
        s += 'sw{solver PBiCGStab;preconditioner  DIC;tolerance 1e-12;relTol 0;}\n'
        s += 'Cw{solver PBiCG;preconditioner  DILU;tolerance 1e-12;relTol 0;}\n'
        s += '\"Cw.*\"{solver PBiCG;preconditioner  DILU;tolerance 1e-12;relTol 0;}\n'
        s += 'Cg{solver PBiCG;preconditioner  DILU;tolerance 1e-12;relTol 0;}\n'
        s += '\"Cg.*\"{solver PBiCG;preconditioner  DILU;tolerance 1e-12;relTol 0;}\n}'
        s += '\nPicard {tolerance 0.01;maxIter 10;minIter 3;nIterStability  5;}'
        s += '\nSIMPLE {residualControl {h 1e-7;} }\n'
        s += 'relaxationFactors { fields { h 0.5;} }\n'
        f1=open(self.fDir+'system\\fvSolution','w');f1.write(s);f1.close()
        
    def getVariable(self,modName,line):
        return ravel(self.core.getValueLong(modName,line,0)[-1::-1])
    
    def writeConstantFields(self):
        '''
        !!!now eps comes from mt3dms and K from modflow (and Kh/Kv)
        '''
        core = self.core
        # gravity field g
        s = 'FoamFile{version 2.0;format ascii;class uniformDimensionedVectorField;'
        s += 'location \"constant\"; object g;}\n'
        s += ' dimensions [0 1 -2 0 0 0 0];\n value  '
        if self.orientation == 'z':s += '( 0 0 -9.81 );'
        elif self.orientation in ['y','r']: s += '( 0 -9.81 0 );'
        f1=open(self.fDir+'constant/g','w');f1.write(s+'\n}');f1.close()
        # porosity
        self.eps = self.getVariable('OpenTrans','poro.1')
        self.writeScalField('constant','eps',self.eps,self.bcD0)
        # permeability
        K = self.getVariable('OpenFlow','khy.2')/86400/9.81e6;self.K=K
        self.writeScalField('constant','Kh',K,self.bcD0,dim='[0 2 0 0 0 0 0]')
        vrt = self.getVariable('OpenFlow','khy.3');kv=vrt*1
        for il in range(self.nlay):
            r = range((self.nlay-il-1)*self.ncell_lay,(self.nlay-il)*self.ncell_lay)
            if self.core.dicval['OpenFlow']['khy.1'][0] == 0:  # ratio
                kv[r] = vrt[r]/86400/9.81e6
            else :
                kv[r] = K[r]/vrt[r]
        self.writeScalField('constant','Kv',kv,self.bcD0,dim='[0 2 0 0 0 0 0]')
        # thickness and bottom
        zb = core.Zblock
        thk = zb[:-1]-zb[1:];thk=thk[-1::-1];thk=ravel(thk);self.thk = thk  #-1::-1 for modflow
        self.writeScalField('constant','thk',thk,self.bcD0,dim='[0 1 0 0 0 0 0]')
        zbot = ravel(zb[1:][-1::-1])
        self.writeScalField('constant','zbot',zbot,self.bcD0,dim='[0 1 0 0 0 0 0]')
        # transport properties
        s = 'FoamFile{version 2.0;format ascii;class dictionary;location "constant";object transportProperties;}\n'
        if self.group == 'Flow' :s += 'activateTransport 0.;\n'
        elif self.group == 'Transport' :s += 'activateTransport 1.;\n'
        elif self.group == 'Chemistry' :s += 'activateTransport 3.;\n'
        if core.dicval['OpenFlow']['dis.8'][3] == 0: s += 'flowStartSteady 1;\n'
        else :s += 'flowStartSteady 0;\n'
        if core.dicaddin['Model']['type'] == 'Unconfined': s += 'activateUnconfined 1;\n'
        if core.dicaddin['Model']['type'] == 'Unsaturated': s += 'activateUnsat 1;\n'
        if core.dicaddin['Model']['type'] == '2phases': s += 'activate2phase 1;\n'
        if core.dicaddin['Model']['type'] in ['Unsaturated','2phases']:
            s += 'activateCapillarity 1;\n'
            s += 'sw_min '+str(core.dicval['OpenFlow']['uns.2'][0])+';\n'
            s += 'alpha_vg alpha_vg [0 -1 0 0 0 0 0] '+str(core.dicval['OpenFlow']['uns.3'][0])+';\n'
            s += 'n_vg '+str(core.dicval['OpenFlow']['uns.4'][0])+';\n'
        s += 'phase.w{rho	rho [1 -3 0 0 0 0 0] 1e3;mu mu [1 -1 -1 0 0 0 0] 1e-3;}\n'
        s += 'phase.g{rho	rho [1 -3 0 0 0 0 0] 1.2;mu mu [1 -1 -1 0 0 0 0] 1e-5;}\n'
        if self.group in ['Chemistry','Transport']:
            s += 'alphaL alphaL [0 1 0 0 0 0 0] '+str(core.dicval['MfUsgTrans']['bct.6'][0])+';\n'
            s += 'alphaT alphaT [0 1 0 0 0 0 0] '+str(core.dicval['MfUsgTrans']['bct.7'][0])+';\n'
        s += 'nlay '+str(self.nlay)+'; ncell_lay '+str(self.ncell_lay)+';'
        f1=open(self.fDir+'constant\\transportProperties','w');f1.write(s);f1.close()
    
    def writeInitFields(self):
        # Uw 
        # considering recharge
        '''
        if ('rch.2' in self.core.diczone['Modflow'].dic.keys()) or (self.core.dicval['Modflow']['rch.2'][0] != 0):
            lcell,ncell,zcell = self.writeOptionCells('rch.2','hrch')
            vbase = self.core.dicval['Modflow']['rch.2'][0]
            if 'rch.2' in self.ttable.keys():
                rmat = self.ttable['rch.2'].astype(float);nt,nc = shape(rmat)
                rmat  = c_[rmat,ones((nt,1))*vbase]# for the background
            else :
                rmat  = ones((nt,1))*vbase# for the background  
            uz = zeros(sum(ncell))
            for i in range(len(ncell)):
                nc = (array(lcell[i])-self.ncell_lay*(self.nlay-1)).astype('int')
                uz[nc] = -float(rmat[0,i])/86400
            self.writeVectField('0','Uw',[0,0,0],self.bcD0,'[0 1 -1 0 0 0 0]',uz=uz)
        else :
        ''' 
        self.writeVectField('0','Uw',[0,0,0],self.bcD0,'[0 1 -1 0 0 0 0]')
        self.writeScalField('0','sw',1,self.bcD0,'[0 0 0 0 0 0 0]')
        self.writeVectField('0','Ug',[0,0,0],self.bcD0,'[0 1 -1 0 0 0 0]')
           
# head h
        h0 = self.getVariable('OpenFlow','head.1')
        name = 'h'
        if self.core.dicaddin['Model']['type'] == 'Unsaturated': name ='hp'
        self.writeScalField('0',name,h0,self.dicBC,dim='[0 1 0 0 0 0 0]')

    def writeTransport(self):
        C0 = self.getVariable('OpenTrans','conc.2') # initial concentrations
        self.writeScalField('0','Cw',C0,self.bcD0,dim='[1 -3 0 0 0 0 0]')
        cactiv  = arange(self.ncell) #cell number but -1 for inactive cells
        if 'conc.1' in self.core.diczone['OpenTrans'].dic.keys(): # active zone
            c_bc = maximum(self.getVariable('OpenTrans','conc.1'),0)
            self.cactiv = cactiv[c_bc==1];
        else :
            self.cactiv = cactiv.astype('int')
        s = '\n'.join(self.cactiv.astype('str'))
        f1=open(self.fDir+'constant\\options\\cactive','w');f1.write(s);f1.close()
        
    def writeChemistry(self):
        self.writePhreeqc(self.core,self.fDir)
        ractiv  = arange(self.ncell) #cell number but -1 for inactive cells
        if 'conc.1' in self.core.diczone['OpenTrans'].dic.keys():
            c_bc = maximum(self.getVariable('OpenTrans','conc.1'),0)
            self.ractiv = ractiv[c_bc==1];
        else :
            self.ractiv = ractiv.astype('int')
        s = '\n'.join(self.ractiv.astype('str'))
        f1=open(self.fDir+'constant\\options\\ractive','w');f1.write(s);f1.close()
        self.writePhqFoam()
        
    def writeOptions(self):
        '''
        '''
        # writing the fvOption file
        ncl = self.ncell_lay
        so ='FoamFile{version 2.0;format ascii;class dictionary;location \"constant\";object fvOptions;}\n'
        if self.core.dicaddin['Model']['type'] == '2phases': vr = 'p'
        else : vr = 'h'
        if 'head.2' in self.core.diczone['OpenFlow'].dic.keys():
            so += vr+'Fix \n{type scalarmyFixedValueConstraint;\nactive true;\n'
            so += 'selectionMode cellSet;\ncellSet '+vr+'fix;\nvolumeMode absolute;\n'
            so += 'fieldValues {'+vr+' 1;}\n}\n'
            # writing the hfix file in cellSet (cells) and data in options
            lcell,ncell,zcell = self.writeOptionCells('head.2',vr+'fix')
            hmat = self.ttable['head.2'];nr,nz = shape(hmat)
            mult = []
            for iz in range(len(lcell)) : 
                mult.append([1]*ncell[iz])
            self.writeOptionData(hmat,lcell,ncell,zcell,'hfix',mult=mult)
        if 'wel.1' in self.core.diczone['OpenFlow'].dic.keys():
            so += vr+'Souwel \n{type scalarmySemiImplicitSource;\nactive true;\n'
            so += 'selectionMode cellSet;\ncellSet '+vr+'wel;\nvolumeMode absolute;\n'
            so += 'injectionRateSuSp {sw (1 0);'+vr+' (1 0);}\n}\n'
            # writing the hWel file in cellSet (cells) and data in options
            lcell,ncell,zcell = self.writeOptionCells('wel.1',vr+'wel')
            qmat = self.ttable['wel.1'];nr,nz = shape(qmat)
            mult = [];eps=0.25
            for iz in range(nz) : 
                kr = self.getPermScaled(lcell[iz]);#print(kr)#/self.thk[lcell[iz]]
                #area = self.area[mod(lcell[iz],ncl)]
                mult.append(kr) #/area/self.thk[lcell[iz]])
            self.writeOptionData(qmat,lcell,ncell,zcell,'hwel',mult=mult)
            
            # can be replaced by change of top bc in Uw
        if ('rch.1' in self.core.diczone['OpenFlow'].dic.keys()) or (self.core.dicval['OpenFlow']['rch.1'][0] != 0):
            so += 'hSourch \n{type scalarmySemiImplicitSource;\nactive true;\n'
            so += 'selectionMode cellSet;\ncellSet hrch;\nvolumeMode absolute;\n'
            so += 'injectionRateSuSp {sw (1 0);h (1 0);}\n}\n'
            lcell,ncell,zcell = self.writeOptionCells('rch.1','hrch')
            vbase = self.core.dicval['OpenFlow']['rch.1'][0]
            if 'rch.2' in self.ttable.keys():
                rmat = self.ttable['rch.2'].astype(float);nr,nc = shape(rmat)
                rmat  = c_[rmat,ones((nr,1))*vbase]# for the background
            else :
                rmat  = ones((nr,1))*vbase# for the background                
            mult = [0]*len(lcell)
            for iz in range(len(lcell)):
                if len(lcell[iz])>0:
                    mult[iz] = self.carea[mod(lcell[iz],ncl)]# thickness will be calc in fvOptions
            self.writeOptionData(rmat,lcell,ncell,zcell,'hrch',mult=mult)
            # write thk in constant/options

        for n in ['riv','drn','ghb']:
            flgAr = False # classical case, conductance given by cell
            if 'conductance' in self.core.dicaddin['Model'].keys():
                if self.core.dicaddin['Model']['conductance'] == 'byZone': flgAr = True
            if n+'.1' in self.ttable.keys():
                so += 'hSou'+n+' \n{type scalarmySemiImplicitSource;\nactive true;\n'
                so += 'selectionMode cellSet;\ncellSet h'+n+';volumeMode absolute;\n'
                so += 'injectionRateSuSp {sw (1 0);h (1 1);}\n}\n'
                lcell,ncell,zcell = self.writeOptionCells(n+'.1','h'+n)
                mat = self.ttable[n+'.1'];nr,nz = shape(mat);
                mult, mult2 = [],[]
                dicz = self.core.diczone['OpenFlow'].dic[n+'.1']
                for iz in range(len(lcell)) : 
                    if len(lcell[iz])==0:
                        mult.append([]);mult2.append([]);continue
                    #zarea = sum(self.area[mod(lcell[iz],ncl)])
                    #if flgAr : ar = self.area[mod(lcell[iz],ncl)]/zarea# conducatnace already given by zone
                    #else : ar=1
                    a = dicz['value'][iz].split('$')[1].split('\n')
                    zz,hcond = float(a[0]),float(a[1])
                    #if zcell[iz]==0: zcell[iz]=[zz]*ncell[iz]
                    mult.append(-hcond)
                    mult2.append(hcond) #
                self.writeOptionData(mat,lcell,ncell,zcell,'h'+n,mult=mult,mult2=mult2)
        if self.group == 'Transport':
            if 'conc.3' in self.ttable.keys():
                so += 'cfix {type scalarmyFixedValueConstraint; active true;'
                so += 'selectionMode cellSet; cellSet cfix; volumeMode absolute;'
                so += 'fieldValues {Cw 1;}}\n'
                lcell,ncell,zcell = self.writeOptionCells('conc.3','cfix')
                mat = self.ttable['conc.3']
                self.writeOptionData(mat,lcell,ncell,zcell,'cfix',formt='float')
            if 'wel.1' in self.ttable.keys():# pumping or injection
                so += 'cwel {type scalarmySemiImplicitSource; active true;'
                so += 'selectionMode cellSet; cellSet cwel; volumeMode absolute;'
                so += 'injectionRateSuSp {sw (1 0); Cw (1 1);}}\n'
                lcell,ncell,zcell = self.writeOptionCells('wel.1','cwel')
                mat = self.ttable['wel.1']
                mat = zeros(shape(mat)) # A dummy way to set all solutions to 0
                if 'cwell.1' in self.ttable.keys(): # find the injecting wells with conc
                    mat1 = self.ttable['cwell.1']
                    lcell1,ncell1,zcell1 = self.writeOptionCells('cwell.1','dum')
                    for iw,w in enumerate(lcell1):
                        if w in lcell: 
                            mat[:,lcell.index(w)] = mat1[:,iw]
                self.writeOptionData(mat,lcell,ncell,zcell,'cwel',formt='int')
            if ('crch.1' in self.core.diczone['MfUsgTrans'].dic.keys()) or (self.core.dicval['MfUsgTrans']['crch.1'][0] != 0):
                so += 'crch {type scalarmySemiImplicitSource; active true;'
                so += 'selectionMode cellSet; cellSet crch; volumeMode absolute;'
                so += 'injectionRateSuSp {sw (1 0); Cw (1 0);}}\n'
                lcell,ncell,zcell = self.writeOptionCells('crch.1','crch')
                vbase = self.core.dicval['MfUsgTrans']['crch.1'][0]
                if 'crch.1' in self.ttable.keys():
                    rmat = self.ttable['crch.1'].astype(float);nr,nc = shape(rmat)
                    rmat  = c_[rmat,ones((nr,1))*vbase]# for the background
                else :
                    rmat  = ones((nr,1))*vbase# for the background                
                self.writeOptionData(rmat,lcell,ncell,zcell,'crch',formt='float')
            for n in ['ghb.1','drn.1','head.2','riv.1']:
                n0 = n[:3]
                if n in self.ttable.keys():# if flow ghb cells put 0 at conc
                    so += 'c'+n0+' {type scalarmyFixedValueConstraint; active true;'
                    so += 'selectionMode cellSet; cellSet c'+n0+'; volumeMode absolute;'
                    so += 'fieldValues {Cw 1;}}\n'
                    lcell,ncell,zcell = self.writeOptionCells(n,'c'+n0)
                    mat = self.ttable[n]
                    mult = [0]*len(lcell)
                    for iz in range(len(lcell)): mult[iz] = [0]*len(lcell[iz])
                    self.writeOptionData(mat,lcell,ncell,zcell,'c'+n0,mult=mult,formt='float')

        if self.group == 'Chemistry':
            ## need to know if ph.4 is at a flow bc or at a well
            self.solucell=[] # solucell conaint the list of fixed chem cells
            if 'ph.4' in self.ttable.keys():
                for i in range(self.ncomp):
                    so += 'cfix'+str(i)+' {type scalarmyFixedValueConstraint; active true;'
                    so += 'selectionMode cellSet; cellSet cfix; volumeMode absolute;'
                    so += 'fieldValues {Cw'+str(i)+' 1;}}\n'
                lcell,ncell,zcell = self.writeOptionCells('ph.4','cfix')
                self.solucell,self.solunb = lcell,zcell
                mat = self.ttable['ph.4']
                self.writeOptionData(mat,lcell,ncell,zcell,'cfix',formt='int')
            if 'ph.5' in self.ttable.keys():# ph recharge
                for i in range(self.ncomp):
                    so += 'crch'+str(i)+' {type scalarmySemiImplicitSource; active true;'
                    so += 'selectionMode cellSet; cellSet crch; volumeMode absolute;'
                    so += 'injectionRateSuSp {sw (1 0); Cw'+str(i)+' (1 0);}}\n'
                lcell,ncell,zcell = self.writeOptionCells('ph.5','crch')
                nz = len(lcell)
                vbase = self.core.dicval['Pht3d']['ph.5'][0]
                rmat = self.ttable['ph.5'].astype(float);nr,nc = shape(rmat)
                rmat  = c_[rmat,ones((nr,1))*vbase]# for the background
                self.writeOptionData(rmat,lcell,ncell,'crch',formt='int')
            if 'wel.1' in self.ttable.keys():# pumping or injection
                for i in range(self.ncomp):
                    so += 'cwel'+str(i)+' {type scalarmySemiImplicitSource; active true;'
                    so += 'selectionMode cellSet; cellSet cwel; volumeMode absolute;'
                    #so += 'fieldValues {Cw'+str(i)+' 1;}}\n'
                    so += 'injectionRateSuSp {sw (1 0); Cw'+str(i)+' (1 1);}}\n'
                lcell,ncell,zcell = self.writeOptionCells('wel.1','cwel')
                mat = self.ttable['wel.1']
                mat = zeros(shape(mat)) # A dummy way to set all solutions to 0
                if 'ph.4' in self.ttable.keys(): # find the injecting wells with conc
                    mat1 = self.ttable['ph.4']
                    lcell1,ncell1,zcell1 = self.writeOptionCells('ph.4','dum')
                    for iw,w in enumerate(lcell1):
                        if w in lcell: 
                            mat[:,lcell.index(w)] = mat1[:,iw]
                self.writeOptionData(mat,lcell,ncell,zcell,'cwel',formt='int')

            for n in ['ghb.1','drn.1','head.2','riv.1']:
                n0 = n[:3]
                if n in self.ttable.keys():# if flow ghb cells put 0 at conc
                    for i in range(self.ncomp):
                        so += 'c'+n0+str(i)+' {type scalarmyFixedValueConstraint; active true;'
                        so += 'selectionMode cellSet; cellSet c'+n0+'; volumeMode absolute;'
                        so += 'fieldValues {Cw'+str(i)+' 1;}}\n'
                    lcell,ncell,zcell = self.writeOptionCells(n,'c'+n0)
                    mat = self.ttable[n]
                    mult = [0]*len(lcell)
                    for iz in range(len(lcell)): mult[iz] = [0]*len(lcell[iz])
                    self.writeOptionData(mat,lcell,ncell,zcell,'c'+n0,mult=mult,formt='float')
                
        f1=open(self.fDir+'constant\\fvOptions','w');f1.write(so);f1.close()

    def writeOptionCells(self,line,fname):
        '''
        writes the cell numbers in constant/polyMesh/sets, lcell is the list of cells
        for each considered zone, ncell the nb of these cells
        '''
        #if fname[:1]=='h': modName = 'Modflow'
        #if fname[:1]=='c' : modName = 'MfUsgTrans'
        #if line[:3] in ['bas','wel','rch','ghb','riv','drn']: modName = 'Modflow'
        #if line[:3] in ['bct','pcb','cwe','crc','cgh','cri']: modName = 'MfUsgTrans'
        #if line[:2]=='ph': modName = 'Pht3d'
        if line[:3] in ['hea','wel','rch','ghb','riv','drn']: modName = 'OpenFlow'
        if line[:3] in ['bct','pcb','con','por','crc']: modName = 'OpenTrans'
        if line[:2]=='ph': modName = 'Pht3d'
        d0 = self.core.diczone[modName].dic
        if line in d0.keys(): dicz = d0[line]
        else : dicz={'name':[]}
        nclay,nlay = self.ncell_lay,self.nlay
        if self.MshType ==0:
            lcell,ncell,zcell = [],[],[]
            nx,ny,xv,yv = getXYvects(self.core)
            dx,dy = xv[1:]-xv[:-1],yv[1:]-yv[:-1]
            dxm,dym = meshgrid(dx,dy)
            self.carea = ravel(dxm*dym)
            for media in range(self.nlay):
                m = zone2grid(self.core,modName,line,media)[-1::-1,:] # modflow inversion
                m1 = zone2grid(self.core,modName,line,media,opt='zon')[-1::-1,:] # modflow inversion
                lzon = unique(m1).nonzero()[0]
                for iz in lzon:
                    if iz>len(lcell): 
                        lcell.append([]);zcell.append([]);ncell.append(0)
                    idx = where(m1==iz) # id0 layer, id1 cell
                    lcell[iz-1].extend(list(idx[0]*nx+idx[1]))
                    zcell[iz-1].extend(list(m[idx[0],idx[1]]))
                    ncell[iz-1] += len(idx[0])
        else :#unstructured
            lcell,ncell,zcell = [],[],[]
            for iz in range(len(dicz['name'])):
                lmedia = dicz['media'][iz]
                if type(lmedia) != type([5,6]) : lmedia = [lmedia]
                lc0,zval = zmesh(self.core,dicz,lmedia[0],iz)
                if len(lc0)==1: lc0 = lc0[0]
                lc1 = []
                for im in lmedia:
                    lc1.extend(lc0+nclay*(self.nlay-im-1));  # nlay-im because inversion vs modflow
                lcell.append(lc1);zcell.append(zval);ncell.append(len(lc0)*len(lmedia))
        # remove cells that can be in two zones, the last one is to be kept
        nl = len(lcell)
        for i in range(nl-1,-1,-1):
            l0 = lcell[i]
            for j in range(i-1,-1,-1):
                a = set(lcell[j])
                a.difference_update(l0)
                lcell[j] = list(a)
                '''
                b = set(zcell[j])
                b.difference_update(l0)
                zcell[j] = list(b)
                ncell[j] = len(a)'''
        if fname in ['hrch','crch']: # to add the base value
            lc_toplay = set(list(range(nclay*(nlay-1),nclay*nlay)))
            lctot = []
            for lc in lcell: lctot.extend(lc)
            lc_toplay.difference_update(lctot)
            lcell.append(list(lc_toplay))
            ncell.append(len(lc_toplay))
            zcell.append(0)
        s ='FoamFile{version 2.0;format ascii;class cellSet;location \"constant/polyMest/sets\";object '+fname+';}\n'
        s += '\n'+str(sum(ncell))+'\n(\n'
        s += '\n'.join([str(a) for lc1 in lcell for a in lc1])+'\n)'
        fD1 = self.fDir+'constant\\polyMesh\\'
        if 'sets' not in os.listdir(fD1): os.mkdir(fD1+'sets')
        f1=open(self.fDir+'constant\\polyMesh\\sets\\'+fname,'w');f1.write(s);f1.close()
        return lcell,ncell,zcell
    
    def writeOptionData(self,mat,lcell,ncell,zcell,fname,mult=None,mult2=None,formt='float'): 
        '''
        write the option semiImplcit for one variable, mult is a list of vector (for well)
        '''
        nr,a = shape(mat);nzo=len(lcell);#dt = 86400
        if mult==None: mult = [[1]*nc for nc in ncell]
        if mult2==None: mult2 = [[0]*nc for nc in ncell]
        if len(zcell)==0: zcell= [[0]*nc for nc in ncell]
        #if formt == 'float': fmt0 = '0 %9i %e %e\n';fmt1='%9i %9i %e %e\n'
        #elif formt == 'int': fmt0 = '0 %9i %9i %9i\n';fmt1='%9i %9i %9i %9i\n'
        s = '%5i %8i \n'%(self.nlay,self.ncell_lay)
        val = zeros((sum(ncell),4))
        #float_formatter = "{:.3e}".format
        #out ='\n'.join(map('{:.3e}'.format, val[:5]))
        #a = ['         ']*sum(ncell)*3
        sc1 = cumsum(ncell);sc0 = r_[0,sc1[:-1]]
        for iz in range(nzo):  # 1st time step
            if len(lcell[iz])==0: continue
            if isinstance(mat[0,iz],float): zz = mat[0,iz]#+array(zcell[iz])
            else : zz = float(mat[0,iz].split()[0])#+array(zcell[iz])
            val[sc0[iz]:sc1[iz],1] = lcell[iz]
            val[sc0[iz]:sc1[iz],2] = zz*array(mult[iz])
            val[sc0[iz]:sc1[iz],3] = array(mult2[iz])
        ag = argsort(val[:,1]) # to sort by cell number
        val1 =val[ag]*1
        s += '\n'.join([' '.join(x) for x in val1.astype('str')])
        s += '\n'
        for it in range(1,nr):
            if all(mat[it] == mat[it-1]): continue
            for iz in range(nzo):
                if isinstance(mat[0,iz],float): zz = mat[it,iz]
                else : zz = float(mat[it,iz].split()[0])
                val[sc0[iz]:sc1[iz],2] = zz*array(mult[iz])
                val[sc0[iz]:sc1[iz],3] = array(mult2[iz])
            val[:,0] = self.tlist[it]+0.001
            val1 =val[ag]*1
            s += '\n'.join([' '.join(x) for x in val1.astype('str')])
            s+='\n'
        val[:,0] = self.tlist[-1]+0.01
        val[:,1:3] = 0
        s += '\n'.join([' '.join(x) for x in val.astype('str')])
        fD1 = self.fDir+'constant\\'
        if 'options' not in os.listdir(fD1): os.mkdir(fD1+'options')
        f1=open(self.fDir+'constant\\options\\'+fname,'w');f1.write(s);f1.close()
        print('option '+fname+' written')
        
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
        
    def writeVectField(self,loc,name,values,bcDict,dim='[0 0 0 0 0 0 0]',uz=[]):
        '''
        values is a list of three values
        case of uz !=None is for recharge on the top boundary
        '''
        s ='FoamFile{version 2.0;format ascii;class volVectorField;location '+loc+';object '
        s += name+';}\n\n'
        s += 'dimensions '+dim+';\n\n'
        s += 'internalField uniform ('+' '.join(['%g'%a for a in values])+');\n'
        s += '\nboundaryField\n{'
        for k in bcDict.keys():
            if (k == 'top')&(len(uz) > 0): # for recharge
                s += '\n  top\n  {\n'
                s += 'type fixedValue;\nvalue nonuniform List<vector> '
                s += str(len(uz)) +'( ('
                u = c_[uz*0,uz*0,uz]
                s1 = ')\n('.join([' '.join(x) for x in u.astype('str')])
                s += s1+') ); }'
            else : # all other cases
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
                
    def findSpecies(self,core):
        '''find the components and species for phreeqc'''
        listE = core.addin.pht3d.getDictSpecies()
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
    
    def writePhreeqc(self,core,fDir):
        """this routine writes a phreeqc file where all solutions are written in
        phreqc format to be able to test their equilibrium before running pht3d
        1. tabke background, then cycle through pht3d zones
        2. get the solution number, phase number...it does not take rates
        3 write them in phreeqc format"""
        s = 'Database '+fDir+'\pht3d_datab.dat \n'
        s += 'Selected_output \n  -totals '+' '.join(self.lspec)+'\n\n'
        listE = core.addin.pht3d.getDictSpecies()
        chem = core.addin.pht3d.Base['Chemistry']
        ncell = self.ncell_lay
        solu = chem['Solutions'];
        dicz = core.diczone['Pht3d']
        nsolu = 0
        if 'ph.3' in dicz.dic.keys():
            nsolu = max([int(a.strip()[0]) for a in dicz.dic['ph.3']['value']])+1
        if 'ph.4' in dicz.dic.keys():
            nsolu = max(nsolu,max([int(a) for a in dicz.dic['ph.4']['value']])+1)
        else :
            nsolu = 2
        self.nsolu = nsolu
        listS = listE['i'];listS.extend(listE['k']);listS.extend(listE['kim'])
        #listS.sort()
        for isol in range(nsolu):
            s += '\nSolution '+str(isol)+' \n units mol/L \n'
            for esp in listS: # go through phase list
                ie=solu['rows'].index(esp);#print esp,phases['rows'],ip,phases['data'][ip] # index of the phase
                conc=solu['data'][ie][isol+1] #backgr concentration of species
                s += esp+' '+str(conc)+'\n'
        pht = self.getVariable('Pht3d','ph.3') #this is a vector ncell
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
                s += esp+' '+str(IS)+' '+str(float(conc)/unique(self.eps)[0])+'\n' #
        # exchanger
        exc=chem['Exchange'];
        nbex = len(unique(dInd['Exchange']))
        for ix in range(nbex):
            if len(listE['e'])>0 :s += '\nExchange '+str(ix)+'\n'
            for esp in listE['e']: # go through phase list
                ie = exc['rows'].index(esp);#print esp,phases['rows'],ip,phases['data'][ip] # index of the phase
                conc = exc['data'][ie][ix+1] #backgr SI and concentration of phase
                s += esp+' '+str(float(conc)/unique(self.eps)[0])+'\n' 
            if len(listE['e'])>0 :s += '-equilibrate with solution '+str(ix)+'\n'
        # surface
        surf=chem['Surface'];
        nbsu = len(unique(dInd['Surface']))
        for isu in range(nbsu):
            if len(listE['s'])>0 :s += '\nSurface '+str(isu)+'\n'
            for esp in listE['s']: # go through phase list
                ie = surf['rows'].index(esp);#print esp,phases['rows'],ip,phases['data'][ip] # index of the phase
                conc = surf['data'][ie][isu+1] #backgr SI and concentration of phase
                s += esp+' '+str(float(conc)/unique(self.eps)[0])+' ' #
                s += ' '.join(surf['data'][ie][6:8])+' \n'
            if len(listE['s'])>0 :s += '-equilibrate with solution '+str(isu)+'\n-no_edl\n'
        # gas phase
        gases=chem['Gases'];
        #nbg = len(unique(dInd['Gases']))
        for ig in range(self.nsolu):
            if len(listE['g'])>0 : s += '\nGas_Phase '+str(ig)+'\n'
            for esp in listE['g']: # go through phase list
                ie = gases['rows'].index(esp);#print esp,phases['rows'],ip,phases['data'][ip] # index of the phase
                IS,conc = gases['data'][ie][ig+1:ig+3] #backgr SI and concentration of phase
                s += esp+' '+str(IS)+' '+str(float(conc))+'\n' #
        # kinetics (up to now just one)
        rates = chem['Rates'];
        parmk = core.addin.pht3d.calcNbParm()
        lkin = listE['k'];lkin.extend(listE['kim'])
        if len(lkin)>0 :s += '\nKinetics 1-'+str(ncell)+'\n'
        for nom in lkin:
            iek = rates['rows'].index(nom);
            s += '\n'+nom+'\n -parms '  #param k
            for ip in range(parmk[nom]):
                s += str(rates['data'][iek][ip+2])+' '
            s += '\n-formula '+rates['data'][iek][-1] +'\n' # formula
        rates = chem['Kinetic_Minerals'];
        if len(listE['kp'])>0 :s += '\nKinetics 1\n'
        for nom in listE['kp']:
            iek = rates['rows'].index(nom);
            iep = phases['rows'].index(nom);
            s += '\n'+nom+'\n -m0 '+ str(phases['data'][iep][2])+'\n'  
            s += ' -parms '
            for ip in range(parmk[nom]): s += str(rates['data'][iek][ip+1])+' '
        f1 = open(fDir+os.sep+'initChem.pqi','w')
        f1.write(s);f1.close()

    def writePhqFoam(self):
        '''phqfoam contains the nb of cell, nb of mobile components and species
        below the lines are 
        Solution,Equilibrium phases,Exchange,Surface,Gas phase,Solid solutions,Kinetics
        dInd contain vectors of solu/assembl nb in the domain, as a dict
        here we also produce the solu_cell file
        '''
        listE = self.core.addin.pht3d.getDictSpecies();
        pht = self.getVariable('Pht3d','ph.3')
        pht = ravel(pht)
        pht = pht[self.ractiv]
        dInd={'Solutions':(pht/1000).astype('int'),
              'Phases':mod(pht,1000)/100,
              'Gases':mod(pht,1000)/100,
              'Exchange':mod(pht,100)/10,
              'Surface':mod(pht,10)}
        shortn=['k','i','kim','g','p','e','s','kp']
        longn=['Solutions','Solutions','Solutions','Gases','Phases','Exchange','Surface',
               'Phases']
        nreac = len(self.ractiv)
        # now write phqfoam        
        # initial solutions
        val0 = ravel(dInd[longn[0]]) 
        nsolu = len(unique(val0))
        s = str(nreac)+' '+str(self.ncomp)+' 0 '+str(nsolu)+' 1\n' # solid-units=1 mol/l here
        if len(unique(val0))==1:  # just one solution or phase...
            s += '0  0 \n'  
        else :
            s += '-1 '+' '.join([str(int(a)) for a in val0]) + '\n'
        #phase, exchange, surface, gases
        for sp in ['p','e','s','g']:
            indx = shortn.index(sp)
            val = ravel(dInd[longn[indx]]);
            if len(listE[sp])>0: # there is at least one speces in this group
                if len(unique(val))==1:  # just one solution or phase...
                    s += '0  0 \n'  
                else :
                    s += '-1 '+' '.join([str(int(a)) for a in val]) + '\n'
            else:
                s += '0  -1 \n'
        s += '0  -1 \n' # solid solution not used
        #kinetics
        flgK = 0
        for sp in ['k','kim','kp']: # special as k and kp are in the same category
            if len(listE[sp]) > 0: flgK = 1
        if flgK : s += '0  1 \n'
        else : s += '0  -1 \n'
        f1 = open(self.fDir+os.sep+'phqfoam.txt','w')
        f1.write(s);f1.close()
        '''now write the phinit.txt file which makes one cell for each solution'''
        s = str(self.nsolu)+' '+str(self.ncomp)+' 0 '+str(self.nsolu)+' 1\n' # solid-units=1 mol/l here
        s += '-1 '+' '.join([str(int(a)) for a in range(self.nsolu)])+'\n'
        for sp in ['p','e','s','g']:
            if len(listE[sp])>0:
                s += '-1 '+' '.join([str(int(a)) for a in range(self.nsolu)])+'\n'
            else :
                s += '0  -1 \n'
        s += '0  -1 \n'  # solid solutions
        if flgK : s += '0  1\n'
        else : s += '0  -1 \n'
                
        f1 = open(self.fDir+os.sep+'phqinit.txt','w')
        f1.write(s);f1.close()

    def writeGSF(self):
        '''
        GSF file
        "UNSTRUCTURED GWF"
        ncell nlay 1 1
        npoints
        x y z (of points)
        icell x y z (cell coord) lay npts (nb of pts) ip1 ip2 ...(points index)
        '''
        core,ncell = self.core,self.mesh.ncell
        #get geom
        if self.MshType == 0:
        	points,faces,bfaces,fcup = self.opf.opfMeshReg(self.core)
        else :
        	points,faces,bfaces,fcup = self.opf.opfMeshFromUsg(self.core);print('got geom from usg')  
        # To get z
        if core.addin.getDim()=='3D':
        	if 'importArray' in core.dictype['OpenFlow']['dis.6']:
        		zlist = self.getZfromPoints(points)
        	else :
        		zlist = [float(a) for a in core.addin.get3D()['topMedia']]
        		zlist.append(float(core.addin.get3D()['zmin']))
        else : 
        	zlist = [core.dicval['OpenFlow']['dis.6'][0],core.dicval['OpenFlow']['dis.7'][0]]
        zlist = zlist[-1::-1] # layers are reversed, start from 0 in OpF
        #write first line
        nlay=1
        s = '"UNSTRUCTURED GWF"\n'+str(ncell)+' '+str(nlay)+' 1 1 \n'
        s += str(len(points)*2)+'\n'
        #write points
        x0,y0 = amin(points[:,0]),amin(points[:,1])
        nplay=len(points)
        for z in zlist:
        	if type(z)==type(0.): z = [z]*nplay
        	coo = zeros((nplay,3))
        	coo[:,0],coo[:,1],coo[:,2] = points[:,0]-x0,points[:,1]-y0,z
        	s += '\n'.join([' %5g %5g %5g ' %(x[0],x[1],x[2]) for x in coo]) 
        	s += '\n'
        # write the cells uses fcup
        for ic in range(ncell):
            x,y = self.mesh.elcenters[ic];z=5; #to be corrected
            s += str(ic+1)+' %5g %g5 %5g 1 ' %(x,y,z)
            s += str((len(fcup[ic])-1)*2)+' '
            s += ' '.join([str(a+1) for a in fcup[ic][:-1]])+' '
            s += ' '.join([str(a+nplay+1) for a in fcup[ic][:-1]])+'\n'
        # writing the file
        f1=open(core.fileDir+os.sep+core.fileName+'.gsf','w');f1.write(s);f1.close()

class opfoamReader:
    def __init__(self,core,opfoam):
        '''core is the model, opf the openfoam object'''
        self.ncell_lay = opfoam.ncell_lay
        self.core,self.opf = core,opfoam
        self.ttable = core.makeTtable()
        self.tlist = self.ttable['tlist']
        self.nlay = getNlayers(core)
        self.ncell = self.ncell_lay*self.nlay
        self.MshType = self.core.getValueFromName('OpenFlow','MshType')
        
    def readHeadFile(self,core,iper):
        '''reads h as the head'''
        return self.readScalar('h',iper)
        
    def readScalar(self,name,iper=None):
        ''' returns a scalar for all the times reshapes in layers
        if iper =None if one value just one period'''
        def rd1(t,n):
            f1=open(self.core.fileDir+os.sep+str(int(t))+os.sep+n,'r')
            s = f1.read();f1.close()
            if 'nonuniform' in s:
                s1 = s.split('<scalar>')[1].split('(')[1].split(')')[0].split('\n')
                data = np.array(s1[1:-1]).astype('float')
            else :
                data = 0
            return data
        #dt = int((self.tlist[-1]-self.tlist[-2])*86400);
        ntime=len(self.tlist)
        if self.MshType==0:
            grd = self.core.addin.getFullGrid()
            nx,ny = grd['nx'],grd['ny']
            shp0 = (ntime,self.nlay,ny,nx)
            shp1 = (self.nlay,ny,nx)
        else :
            shp0 = (ntime,self.nlay,self.ncell_lay)
            shp1 = (self.nlay,self.ncell_lay)
        if iper==None:
            data = np.zeros((ntime,self.ncell))
            for i,t in enumerate(self.tlist[1:]):
                data[i+1,:]= rd1(t*86400,name)
            V = reshape(data,shp0)
        else:
            data = rd1(self.tlist[iper]*86400,name)
            V = reshape(data,shp1)            
        return V


