#
import os,sys,time,base64,types,io # OA 25/10/18 add types
import subprocess as sbp
from numpy import frombuffer,float64
'''juste below : to remove the . before each module'''
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
#from scipy.interpolate import griddata
from modflowWriter import *
from mtphtWriter import *
from mtUsgWriter import *
#from min3pWriter import *
import xml.dom.minidom as xdom
from geometry import *
from importExport import *
from addin import *
from timeperiod import *
from numpy import loadtxt,size #OA 25/1/21
from numpy import savez_compressed
from numpy import load as npload
from config import *
from modflowKeywords import Mf
from mtPhtKeywords import Mt
from mtUsgKeywords import Mtu  # OA 27/7/19
from pht3dKeywords import Ph
#from min3pFlowKeywords import m3F
#from min3pTransKeywords import m3T
#from min3pChemKeywords import m3C
from opfoamKeywords import OpF
from opfoamKeywords import OpT
from opfoamKeywords import OpC
from obsKeywords import Obs
from pestKeywords import Pst
from plugins import *
from collections import deque

class Core:
    """this is the central class that makes the link between all objects
    it can be used on graphic mode or batche mode
    """
    def __init__(self,gui=None):
        self.modelList = ['Modflow','Mt3dms','MfUsgTrans','Pht3d', # OA 27/7/19 added MfUsgTrans
            'OpenFlow','OpenTrans','OpenChem','Observation','Pest']
        self.gui = gui
        self.baseDir = os.getcwd(); # OA 11/9/18 2 lines below modfiied
        if gui!=None:
            if gui.gtyp=='qgis': self.baseDir= self.gui.plugin_dir
        self.dicval = {}
        self.dictype = {}
        self.diczone = {}
        self.dickword = {}
        self.dicarray = {}
        self.dicformula = {}
        self.dicinterp = {}
        self.dicplugins = {}
        for mod in self.modelList:
            self.dicval[mod] = {}
            self.dictype[mod] = {}
            self.dicarray[mod] = {}
            self.dicformula[mod] = {}
            self.dicinterp[mod] = {} # EV 20/02/20
        # OA 30/07 Mf... are now defined as classes : more correct and usefull for python3
        self.dickword['Modflow'] = Mf()
        self.dickword['Mt3dms'] = Mt()
        self.dickword['MfUsgTrans'] = Mtu()
        self.dickword['Pht3d'] = Ph()
        #self.dickword['Min3pFlow'] = m3F()
        #self.dickword['Min3pTrans'] = m3T()
        #self.dickword['Min3pChem'] = m3C()
        self.dickword['Observation'] = Obs()
        self.dickword['OpenFlow'] = OpF()
        self.dickword['OpenTrans'] = OpT()
        self.dickword['OpenChem'] = OpC()
        self.dickword['Pest'] = Pst()
        self.plugins = plugins(self)
        self.initAll()
        self.mfUnstruct = False # OA 17/9/17 made to work with modflow USG
        self.dicaddin = {}
        self.addin = addin(self)
        self.addin.initAddin()
        self.fileDir,self.fileName = None,None
        self.createKwList()
        
#********************* initialisation ****************************
    def initAll(self):
        self.radfact,self.Zblock,self.grd = 1.,None,None
        self.mfUnstruct,self.ttable,self.lcellInterp = False,None,[] # OA 9/1/21
        self.MshType = 0
        for mod in self.modelList:
            self.dicval[mod] = self.initVal(mod)
            self.dictype[mod] = self.initArray(mod) # ,self.dicarray[mod] #EV 06/02/20
            self.diczone[mod] = dicZone(self,mod)

    def initVal(self,modName):
        val = {}
        lines = self.dickword[modName].lines
        for n in list(lines.keys()):
            if 'default' in lines[n]:
                val[n] = lines[n]['default']
            else : val[n] = [0]*len(lines[n]['kw'])
        return val
    
    def initArray(self,modName):
        array,atype = {},{}
        lines = self.dickword[modName].lines
        for n in list(lines.keys()):
            if lines[n]['type'][0][:3] == 'arr': 
                array[n],atype[n] = None,['one_value']
            else : 
                atype[n]=['one_value']
        return atype#,array #EV 06/02/20
    
    def getFormula(self,modName,line,media): #EV 06/02/20
        """if the formula exists it returns it, if not, it
        creates a void one in the dic and returns it"""
        nmedia = getNmedia(self)
        if line in self.dicformula[modName]: 
            try : self.dicformula[modName][line][media]
            except IndexError :
                form=len(self.dicformula[modName][line])
                nform=len(form)
                if nform<nmedia : form.extend(['None']*(nmedia-nform))
                self.dicformula[modName][line][media] = 'value ='
                return 'value ='
            else : 
                return self.dicformula[modName][line][media]
        else : 
            self.dicformula[modName][line] = ['None']*nmedia
            self.dicformula[modName][line][media] = 'value ='
            return 'value ='
            
    def updateDicts(self):
        """this methods copies the values in one model to another one if 
        the key is the same or if the correspondance is in the list"""
        def copies(nam1,nam2,kM1,kM2):
            for k in list(kM2.keys()):
                if k in list(kM1.keys()):
                    ll1,nb1 = kM1[k]
                    ll2,nb2 = kM2[k]
                    self.dicval[nam2][ll2][nb2] = self.dicval[nam1][ll1][nb1]
        nam1, nam2 = 'Modflow','Mt3dms'
        kM1 = self.KwList[nam1]
        kM2 = self.KwList[nam2]
        copies(nam1,nam2,kM1,kM2)
#        nam2 = 'Pht3d'
#        kM2 = self.KwList[nam2]
#        copies(nam1,nam2,kM1,kM2)
        
    def getUsedModulesList(self,modName):
        return self.addin.getUsedModulesList(modName)
    
    def setUsedModulesList(self,modName,grpList):
        self.addin.setUsedModulesList(modName,grpList)
        
    def makeTtable(self):
        self.ttable = makeTransientTable(self)
        return self.ttable
    
    def getTlist2(self):
        tlist = array(self.ttable['tlist'])
        if self.dicaddin['Model']['group'][:4]=='Open': 
            tlist=array(self.ttable['wtimes']);return tlist
        else : return tlist[1:] #(tlist[:-1]+tlist[1:])/2.
        
#*************************** load, save,run models ********************    
    def openModel(self,fDir,fName):
        """opens an orti file, stored in xml"""
        self.fileDir,self.fileName = fDir,fName
        self.initAll()
        fullName = fDir+os.sep+fName; #print(fullName) #EV 06/01/20
        #if fName == '' : return
        if fName+'.orti' in os.listdir(self.fileDir): fullName +='.orti'
        else : fullName +='.iqpht'
        #flgArr = False
        f1 = open(fullName, 'r');doc = f1.read();f1.close()
        #if 'compressdata.npz' in os.listdir(fDir):
        #    darr = npload(fDir+os.sep+'compressdata.npz');flgArr = True
        dom = xdom.parseString(doc)    
        dicts = dom.getElementsByTagName("dict")
        for d in dicts:
            dname = d.getElementsByTagName("name")[0].childNodes[0].data
            model,typ = dname.split('_')
            if (model not in self.modelList) & (model!='dic'): continue
            keys = d.getElementsByTagName("key");dict1 = {}
            for k in keys:
                kname = str(k.getElementsByTagName("name")[0].childNodes[0].data)
                dict1[kname] = eval(k.getElementsByTagName("content")[0].childNodes[0].data)
                #exec('dict1[kname] = '+kdata);
            if typ=='val': self.dicval[model].update(dict1);#print self.dicval[model]
            elif typ=='type': self.dictype[model].update(dict1)
            elif typ=='zone': self.diczone[model].setDic(dict1)
            elif typ=='array' : 
                self.dicarray[model].update(dict1) #EV 07/02/20
            #elif typ=='array' and flgArr: #EV 07/02/20
                #for k in keys : 
                    #kname = str(k.getElementsByTagName("name")[0].childNodes[0].data)
                    #self.dicarray[model][kname] = darr[kname.replace('.','')]
            elif typ=='formula': 
                self.dicformula[model].update(dict1)
            elif typ=='interp': 
                self.dicinterp[model].update(dict1) # EV 20/02/20
            elif typ=='addin': 
                self.addin.update1(dict1)
            elif typ=='plugins': 
                self.dicplugins=dict1
            #print self.dicaddin
        #self.addin.initAddin() seems to make trouble
        if self.gui!=None: self.gui.onInitGui(self)
        #print('in open')
        self.addin.grd = makeGrid(self,self.dicaddin['Grid']);#print 'core 152',self.addin.grd
        group = self.dicaddin['Model']['group'];self.addin.group = group;print(group)
        self.makeUsedModelList(group)
        self.makeTtable()
        self.MshType = 0                    
        self.mfUnstruct = False  # OA 13/3/21
        if self.gui!=None: self.gui.currentModel = 'Modflow'
        if group == 'Modflow USG': # OA 02/20
            self.mfUnstruct = True
            self.MshType = self.getValueFromName('Modflow','MshType')
            self.addin.setMfUnstruct();
            self.addin.setGridInModel('old')
            #if self.mfUnstruct and self.getValueFromName('Modflow','MshType')>0:#OA 4/3/20   
            #    self.flgMesh = 1 #18/12/20                      
        if group == 'Openfoam': # OA 02/20
            if self.gui!=None: self.gui.currentModel = 'OpenFlow'
            self.MshType = self.getValueFromName('OpenFlow','MshType')
            self.addin.MshType = self.MshType
            self.addin.setGridInModel('old')
            self.dtu,self.lu=self.addin.opfoam.getUnits()
            self.flowReader = opFlowReader(self,self.addin.mesh)
            self.transReader = opTransReader(self,self.addin.mesh)
        if group[:5] == 'Modfl': #OA 02/20
            self.flowReader = modflowReader(fDir,fName)
            if 'USG' in group: self.transReader = mtUsgReader(fDir,fName)
            else : self.transReader = mtphtReader(fDir,fName)
        elif group[:5] == 'Min3p' :
            self.addin.min3p.buildMesh(opt='read')
            self.flowReader = min3pReader(self,fDir,fName)
            self.transReader = min3pReader(self,fDir,fName)
        if self.gui != None and self.MshType>0: self.gui.onGridMesh('Mesh') #EV 30/09/19 # OA removed on 8/2/20
        if type(self.Zblock)!=type(ones(3)): self.Zblock = makeZblock(self)
        self.addin.setChemType()
        #self.usePostfix()
        
    def makeUsedModelList(self,group):
        '''crate a list of models that are used in the opened file'''
        if group == 'Modflow series':
           self.usedModelList=['Modflow','Mt3dms','Pht3d']
        if group == 'Modflow USG':
           self.usedModelList=['Modflow','MfUsgTrans','Pht3d']
        if group[:4] == 'Open':
           self.usedModelList=['OpenFlow','OpenTrans','OpenChem']
        self.usedModelList.append('Observation')
            
    def saveModel(self,fDir = None,fName = None):
        """save the model"""
        if fDir!= None:
            self.fileDir, self.fileName = fDir, fName
        filename = self.fileDir+os.sep+self.fileName + '.orti'
        darray = {}
        farrname = self.fileDir+os.sep+'compressdata.npz'
        f1 = open(filename,'w');str1 = '<ORTi3Ddoc>\n'
        for md in self.modelList:
            for t in ['val','type','zone','formula','array','interp']: # EV 20/02/20
                dic = eval('self.dic'+t+'[md]'); # OA 1/8/17 exec to eval for python3
                str1+= '<dict>\n<name>'+md+'_'+t+'</name>\n'
                if t=='zone' :
                    for k in list(dic.dic.keys()):
                        str1 += '<key><name>'+k+'</name><content>'+str(dic.dic[k])+\
                            '</content></key>\n'
                else :
                    for k in list(dic.keys()):
                        str1 += '<key><name>'+k+'</name><content>'+str(dic[k])+\
                            '</content></key>\n'                    
                str1+= '</dict>\n'             
        str1+= '<dict>\n<name>dic_addin</name>\n'
        for k in list(self.dicaddin.keys()):
            str1 += '<key><name>'+k+'</name><content>'+str(self.dicaddin[k])+\
                '</content></key>\n'
        str1+= '</dict>\n'
        str1+= '<dict>\n<name>dic_plugins</name>\n'
        for k in list(self.dicplugins.keys()):
            str1 += '<key><name>'+k+'</name><content>'+str(self.dicplugins[k])+\
                '</content></key>\n'
        str1+= '</dict>\n'
        str1+= '</ORTi3Ddoc>'
        f1.write(str1)
        f1.close()
        s = ' ' # will be a string giving the name of the keys for saving
        for k in list(darray.keys()): s += k.replace('.','')+'=darray[\''+k+'\'], '
        #print 'core 213',filename,farrname,s
        if s != ' ':exec('savez_compressed(r\''+farrname+'\','+s[:-1]+')') #r to get raw string for \ problem
        
    def writeModel(self,modName,info=True):
        """ writes the ascii file for modflow series, does nothing for fipy"""
        mtype = self.dicaddin['Model']['group']  # OA 10/2/2020
        #print('in core',modName)
        if modName  in ['Modflow','Modflow_USG']:
            self.mfWriter = modflowWriter(self,self.fileDir,self.fileName)
            self.mfWriter.writeModflowFiles(self)
            self.flowReader = modflowReader(self.fileDir,self.fileName)
        if modName in ['Mt3dms','MfUsgTrans','Pht3d']: # OA 28/7/19
            if 'USG' in mtype: #28/7/19 this and 2 below # oa modif 10/2/20
                self.mtWriter = mtUsgWriter(self,self.fileDir,self.fileName)
                self.transReader = mtUsgReader(self.fileDir,self.fileName)    # OA 21/08/19
            else : 
                self.mtWriter = mtphtWriter(self,self.fileDir,self.fileName)
                self.transReader = mtphtReader(self.fileDir,self.fileName)   
            parmk = None
            if modName =='Pht3d':
                dicSpec = self.addin.pht3d.getDictSpecies()
                parmk = self.addin.pht3d.calcNbParm()
            else : 
                dicSpec ={'mcomp':1,'ncomp':1,'gcomp':1,'kim':[]}
            info = self.mtWriter.writeMtphtFiles(dicSpec,modName,parmk) # OA 18/12/21 adde info
        if modName[:5]  == 'Min3p':
            self.m3pWriter = min3pWriter(self,self.fileDir,self.fileName)
            self.m3pWriter.writeMin3pFiles(self,modName[5:])
            self.transReader=self.flowReader = min3pReader(self,self.fileDir,self.fileName)
        if modName  in ['OpenFlow','OpenTrans','OpenChem']:
            self.opfWriter = opfoamWriter(self,self.addin.mesh)
            dicBC= {};options = {'group':modName[4:]}
            self.opfWriter.writeFiles(self.fileDir,dicBC,options)
            self.flowReader = opFlowReader(self,self.addin.mesh)
            self.transReader = opTransReader(self,self.addin.mesh)
        if modName == 'Pest':
            info=self.addin.pest.writeFiles() #EV 11/12/19
        self.writeBatchFile(modName)
        if info ==True :return 'Files written'
        else : return info  #EV 11/12/19   
        
    def runModel(self,modName,info=False):
        if modName in ['Modflow','Modflow_USG']: #self.runProgress('runMflow.bat')
            sbp.run(['runMflow.bat'],creationflags=sbp.CREATE_NEW_CONSOLE)
        if modName == 'Mt3dms': #self.runProgress('runMt3d.bat')
            sbp.run(['runMt3d.bat'],creationflags=sbp.CREATE_NEW_CONSOLE)
        if modName == 'MfUsgTrans': #self.runProgress('runUtrp.bat')
            sbp.run(['runUtrp.bat'],creationflags=sbp.CREATE_NEW_CONSOLE)
        if modName == 'Pht3d': #self.runProgress('runPht3d.bat')
            sbp.run(['runPht3d.bat'],creationflags=sbp.CREATE_NEW_CONSOLE)
        if modName in ['OpenFlow','OpenTrans','OpenChem']: #self.runProgress('runOpf.bat')
            sbp.run(['runOpf.bat'],creationflags=sbp.CREATE_NEW_CONSOLE)
        a = self.returnState(modName,info);
        return a

    def writeBatchFile(self,modName):
        tabRes, sep = [],os.sep
        lf = os.listdir(self.fileDir) 
        cfg = Config(self);#print cfg.typInstall
        #if self.gui != None and 'dist' not in self.baseDir: #OA 22/6 modified
            #if cfg.typInstall=='exe': self.baseDir += os.sep+'dist'
        if modName[:4] not in ['Opge','Min3','Open'] and self.fileName+'.nam' not in lf: 
            return 'Files not written'
        try : 
            b=str(self.baseDir).encode("utf-8")
            b=str(self.fileDir).encode("utf-8")
        except UnicodeEncodeError: return 'Bad caracters in folder name'
        if modName in ['Modflow','Modflow_USG']:
            mod,lastline = 'mf2k_PMwin',3 # OA 22/8/19 added lastline for search line for usg too
            if 'USG' in modName: 
                mod,lastline = 'PHT_USG',7 #EV 19/03/21
            if 'NWT' in self.getUsedModulesList('Modflow'): mod = 'mfNWT_dev'
            if os.name == 'nt':
                exec_name = '"'+self.baseDir+sep+'bin'+sep+mod+'.exe"'
            else :
                exec_name = self.baseDir+sep+'mf2k '
            s=exec_name+' '+self.fileName+'.nam'
            os.chdir(self.fileDir)
            f1 = open('runMflow.bat','w');f1.write(s);f1.close()

        if modName == 'Mt3dms':
            mod1,mod2 = 'mt3dms5b','Mt3dms'
            if 'VDF' in self.getUsedModulesList('Mt3dms'): 
                mod1,mod2 ='swt_v4','Mt3dms'
            s=self.baseDir+sep+'bin'+sep+mod1+'.exe '+mod2+'.nam'
            os.chdir(self.fileDir)
            f1 = open('runMt3d.bat','w');f1.write(s);f1.close()
            
        if modName == 'MfUsgTrans': # OA 19/8/19
            s=self.baseDir+sep+'bin'+sep+'PHT_USG.exe '+self.fileName #EV 19/03/21
            os.chdir(self.fileDir)
            f1 = open('runUtrp.bat','w');f1.write(s);f1.close()
            
        if modName == 'Pht3d':
            if self.dicaddin['Model']['group'] == 'Modflow USG': # OA 03/20
                s = self.baseDir+sep+'bin'+sep+'PHT_USG.exe '+self.fileName
            else :
                if self.dicval['Pht3d']['ph.6'][5]== 1 : #EV 2/7/21
                    N=self.dicval['Pht3d']['ph.6'][6]
                    s='mpiexec -n '+ str(N) +' '+ self.baseDir+sep+'bin'+sep+'pht3dv217_mpi_fett.exe Pht3d.nam'
                else : s=self.baseDir+sep+'bin'+sep+'Pht3dv217.exe Pht3d.nam'
            os.chdir(self.fileDir)
            f1 = open('runPht3d.bat','w');f1.write(s);f1.close()
                    
        if modName in ['OpenFlow','OpenTrans','OpenChem']:
            sbin=self.baseDir+'\\bin'
            s= '@echo off\n call '+sbin+'\\opflib\\mySetvars_OF8.bat\n'
            s += 'cd '+self.fileDir+'\n'
            s += 'call '+sbin+'\muFlowRT.exe >log.txt \n echo.'
            os.chdir(self.fileDir)
            f1 = open('runOpf.bat','w');f1.write(s);f1.close()
            #process = sbp.Popen('runOpf.bat', shell=True) #, stdout=sbp.PIPE)
            
        if modName[:5] =='Min3p':
            #print('name',self.fileName)
            s=self.baseDir+sep+'bin'+sep+'min3p.exe '+self.fileName ; # OA 19/3/19
            os.chdir(self.fileDir)
            p = sbp.Popen(s).wait() #),creationflags=CREATE_NEW_CONSOLE).wait(); #OA 8/6/19
            if info !=False :
                return self.getTxtFileLastLine(self.fileName+'.log',5)
            
        if modName[:5] == 'Pest': #EV 07/11
            if self.dicval['Pest']['ctd.1'][1]==0:
                #s=self.baseDir+sep+'bin'+sep+'pest.exe '+self.fileName #; print(s)
                s=self.baseDir+sep+'bin'+sep+'pestpp-glm.exe '+self.fileName #EV 15/7/21
            else : 
                s=self.baseDir+sep+'bin'+sep+'pestpp-glm.exe '+self.fileName
                #s=self.baseDir+sep+'bin'+sep+'pest.exe '+self.fileName+'r' #; print(s)
            os.chdir(self.fileDir)
            p = sbp.Popen(s).wait()#),creationflags=CREATE_NEW_CONSOLE).wait(); #OA 8/6/19
            #subprocess.call('start /wait '+s, shell=True)

    def runProgress(self,cmd):
        """Returns the nth before last line of a file (n=1 gives last line)"""
        t1=self.dicaddin['Time'];
        lfin,lstp=t1['final'],t1['steps']
        print (self.gui,self.gui.notify.text())
        self.setRunTime('0')
        wtimes=[]
        for i in range(len(lfin)):
            stp,fin=float(lstp[i])*self.dtu,float(lfin[i])*self.dtu
            wtimes.extend(list(linspace(stp,fin,int(fin/stp)).astype('int')))
        s=''
        if 'log.txt' in os.listdir(self.fileDir): os.system("del log.txt")
        f1=open(self.fileDir+'log.txt','a')
        with sbp.Popen(cmd, stdout=sbp.PIPE, bufsize=100,universal_newlines=True) as p:
            for line in p.stdout:
                s +=line
                if 'time =' in line: 
                    if line.split()[2] in wtimes : 
                        self.setRunTime(line.split()[2])
                        f1.write(s);s=''
        rc = p.returncode
        f1.close()
        
    def setRunTime(self,time):
        if self.gui !=None: self.gui.onNotify('run time : '+str(time))
    
    def returnState(self,modName,info):
        if info !=False and modName[:6]=='Modflo':
            try :  # EV 13/11 "show model fail to converge"
                if modName == 'Modflow':time_model=self.dicval['Modflow']['dis.2'][4]
                else : time_model=self.dicval['Modflow']['disu.4'][5] # EV 3/12/21
                if time_model==0 : 
                    time_out=self.getTxtFileLastLine(self.fileName+'.lst',3).split()[4]
                else : 
                    time_out=self.getTxtFileLastLine(self.fileName+'.lst',3).split()[(time_model+1)]
                time_last=self.getTlist2()[-1]
                if float(time_out)==float(time_last):
                    if modName == 'Modflow':
                        return ('Normal termination of MODFLOW-2000')
                    else : return ('Normal termination of Modflow USG') # EV 3/12/21
                else: return('Model fail to converge')
            except :# IndexError:
                return('Model fail to converge')

        if info !=False and modName[:5]=='Mt3dms':
            try : # EV 13/11 "show model fail to converge"
                line_out=self.getTxtFileLastLine('Mt3dms.out',3).split()[4]
                if line_out=='END':
                    return ('Normal termination of MT3DMS')
                else: return('Model fail to converge') #return self.getTxtFileLastLine('Mt3dms.out',3)#+'\n Mt3dms run done'
            except IndexError:
                return('Model fail to converge')

            
        if info !=False and modName[:5]=='Pht3d':
            if self.dicaddin['Model']['group'] != 'Modflow USG': 
                try : # EV 13/11 "show model fail to converge"
                    line_out=self.getTxtFileLastLine('Pht3d.out',3).split()[4]
                    if line_out=='END':
                        return ('Normal termination of PHT3D')
                    else: return('Model fail to converge') 
                except IndexError:
                    return('Model fail to converge')
            else : 
                try :  # EV 3/12/21
                    time_model=self.dicval['Modflow']['disu.4'][5]
                    if time_model==0 : 
                        time_out=self.getTxtFileLastLine(self.fileName+'.lst',7).split()[4]
                    else : 
                        time_out=self.getTxtFileLastLine(self.fileName+'.lst',7).split()[(time_model+1)]
                    time_last=self.getTlist2()[-1]
                    if float(time_out)==float(time_last):
                        return ('Normal termination of PHT-USG')
                    else: return('Model fail to converge')
                except :
                    return('Model fail to converge')
                
        if info !=False and modName[:5]=='MtUsg':
            try :  # EV 3/12/21
                time_model=self.dicval['Modflow']['disu.4'][5]
                if time_model==0 : 
                    time_out=self.getTxtFileLastLine(self.fileName+'.lst',7).split()[4]
                else : 
                    time_out=self.getTxtFileLastLine(self.fileName+'.lst',7).split()[(time_model+1)]
                time_last=self.getTlist2()[-1]
                if float(time_out)==float(time_last):
                    return ('Normal termination of Modflow USG Transport')
                else: return('Model fail to converge')
            except :
                return('Model fail to converge')
                
        if info!=False and modName in ['OpenFlow','OpenTrans','OpenChem']:
            #time_model=int(float(self.dicaddin['Time']['final'][-1])*self.dtu)
            lines= self.getTxtFileLastNLines(self.fileDir+os.sep+'log.txt',3)
            #print(time_model);print(lines)
            for i in range(2,-1,-1):
                if len(lines[i])>6:
                    if lines[i][:6]=='Normal': 
                        return 'Normal termination of OpenFoam'
            return 'Model fail to converge'
            
    def getTxtFileLastLine(self,fname,line):
        f1 = open(fname,'r')
        a= f1.read().split('\n')
        f1.close()
        return a[-line]
   
    def getTxtFileLastNLines(self,fname,nline):
        f1 = open(fname,'r')
        a= f1.read().split('\n')
        f1.close()
        return a[-nline:]
    
    def runZonebud(self,modName): # EV 04/03/20
        s = self.baseDir+os.sep+'bin'+os.sep
        if 'USG' in modName: s += 'zonbudusg.exe'
        else : s += 'zonbud.exe'
        myinput = open(self.fileDir+os.sep+'zonbud.in')
        os.chdir(self.fileDir)
        p = sbp.Popen(s).wait() #),stdin=myinput,creationflags=CREATE_NEW_CONSOLE).wait()
        
#********************** import and export functions *****************
    def importData(self,fileDir,fileName):
        importer = impFile(self.gui,self)
        dicData = importer.impTabFile(fileDir+os.sep+fileName+'.txt')
        self.data = dicData

    #def importSolutions(self,fileDir,fileName):                #EV 14/11/19
    #    importer = impFile(self.gui,self)
    #    dicData = importer.impTabFile(fileDir+os.sep+fileName+'.txt')
    #    self.addin.pht3d.setImportedSolutions(dicData)
        
    def importZones(self,fileDir,fileName,modName,line):
        importer = impFile(self.gui,self)
        if fileName=='': return #EV 27/11
        importer.impZones(fileDir,fileName,modName,line)
        
    def importLayerValues(self,fileName,spname):
        """import from existing txt files the concentrations of one given species 
        (this is for the restart option)
        """
        m0 = loadtxt(self.fileDir+os.sep+fileName) # first line layer numb from top to bott, 2nd :a value
        nlay,a = shape(m0)
        g = self.addin.getFullGrid();nx=g['nx'];dx=g['dx']
        m1 = tile(m0[-1::-1,1:],(1,nx));#print m1
        dim = self.addin.getDim()
        dictE = self.addin.pht3d.getDictSpecies();#print spname,dictE
        for kw in list(dictE.keys()):
            if iterable(dictE[kw])==0: continue
            if spname in dictE[kw] : groupname = kw
        if dim =='Radial' and groupname in ['p','g','e','s']:
            for l in range(nlay): 
                m1[l] = m1[l]*(cumsum(dx)-dx/2.)*6.28;
        return m1
    
    def importAscii(self,fileDir,fileName): #EV 04/02/19
        importer = impFile(self.gui,self)
        if fileName=='': return 
        m0=importer.impAsciiGrid(fileDir,fileName)
        return m0   
    
    def importGridVar(self,fileDir,fileName): #EV 02/04/20
        importer = impFile(self.gui,self)
        if fileName=='': return 
        ysign,zdx,zdy,m0=importer.impGridVar(fileDir,fileName)
        return ysign,zdx,zdy,m0 
            
#********************* working with keywords and values***************            
    def createKwList(self):
        """creates a list of all keywords per model as a dict
        with the line and number for each kw
        """
        self.KwList={}
        for modName in self.modelList:
            self.KwList[modName]={}
            lines = self.dickword[modName].lines
            for ll in list(lines.keys()):
                kw0=lines[ll]['kw']
                for ik in range(len(kw0)):
                    kw=kw0[ik].split('(')[0] # removes the dimension in ()
                    self.KwList[modName][kw]=(ll,ik)
        #print self.KwList
                    
    def setValue(self,modName,line,ik,value):
        """sets a value to a dicvalue place"""
        self.dicval[modName][line][ik] = value
        #print line,value
        if line in ['dis.6','dis.7'] : 
            self.Zblock = makeZblock(self);#print self.Zblock #top and bottom

    def getValueFromName(self,modName,vName,option=-1):
        """returns a value from the name of the keyword"""
        #print modName,vName,self.KwList[modName]
        if (vName =='NLAY') and (self.addin.getDim() in ['Radial','Xsection']):
            return getNlayers(self)
        if vName in list(self.KwList[modName].keys()): 
            line,ik = self.KwList[modName][vName]  # OA 23/7/20 3 lines below modified
            if len(self.dickword[modName].lines[line]['kw'])==1 and option !=-1 : # case where we want a va lue for a given layer
                val = self.dicval[modName][line][option]
            else : 
                if ik<len(self.dicval[modName][line]): val = self.dicval[modName][line][ik];#classical case get the true keyword value
                else : val =None
            return val
            
    def getSingleValueFromName(self,modName,vName,nb): # OA 22/7/20 added nb
        #if type(val)==type([]): return float(val)
        print('hello')
        try : 
            val = self.getValueFromName(modName,vName,nb)
            return float(val)
        except : 
            return 0
        
    def setValueFromName(self,modName,vName,value):
        if vName in self.KwList[modName]:
            line,ik=self.KwList[modName][vName]
            if type(value)==type([5]):
                self.dicval[modName][line] = value;#print vName,value
            else :
                dicv = self.dicval[modName][line]
                if ik<len(self.dicval[modName][line]):
                    self.dicval[modName][line][ik] = value
                else : 
                    self.dicval[modName][line].extend([0]*(ik-len(dicv)+1))
                    self.dicval[modName][line][ik] = value
        #print 'core setvn',vName,type(value),value,self.dicval[modName][line]

    def testCondition(self,modName,cond0,option=0): # OA 22/7/20 adde any option
        """tests if a condition from the dictionnary is validated; 
        it allows to use several conditions in the same line
        option can be set to the layer nb or 'any' """
        a=True;
        if option == 'any' : llay = range(getNlayers(self))
        else: llay=[option]
        kwl=self.KwList[modName];#print kwl
        s1='self.getSingleValueFromName(modName,';
        if cond0=='': return True
        for ilay in llay :
            cond = cond0*1
            cond=cond.replace(' and ',') and (')
            cond=cond.replace(' or ',') or (')
            for k in kwl:
                if k in cond:
                    cond=cond.replace(k,s1+'\''+k+'\','+str(ilay)+')')
            a=eval('('+cond+')');#print(kwl,a) # OA 1/8/17 for python 3 
            if a : return True
        return False  

    def getSizeKw(self,modName,kw):
        """get the size of the object designed by a keyword (vector or array)"""
        if modName == 'Pht3d' : modName = 'Mt3dms' #the variables are located in mt3d
        size = [];
        a = kw.split('(')
        if len(a)==1 : return size
        b = a[1][:-1].split(',') # list of dimensions
        for s in b:
            n = self.getValueFromName(modName,s)
            size.append(int(n))
        return size

    def getValueLong(self,modName,line,ik,iper=0):
        """get the vector or array of a keyword using the size of a vector, 
        an array or a formula"""
        valIn = self.dicval[modName][line][ik]
        kw = self.dickword[modName].lines[line]['kw'][ik]
        #print('vtype', vtype, 'valIn', valIn, 'kw', kw)
        cond = self.dickword[modName].lines[line]['cond']
        if self.testCondition(modName,cond,'any')==False: return None # OA 22/7/20 added any
        #size = self.getSizeKw(modName,kw)
        kw = kw.split('(')[0]
        #print('size',size, 'kw', kw)
        numtype = self.dickword[modName].lines[line]['type'][ik];#int, float or lay
        #print('numtype', numtype)
        if line in ['dis.6','disu.7'] : return self.Zblock[:-1] # several top, OA 25/7/20 added disu
        if line =='btn.9' : return self.Zblock[0] # only the top        
        if line == 'dis.7': return self.Zblock[-1:] # one bottom for mflow 2000
        if line == 'disu.8': return self.Zblock[1:] # all bottom for usg added OA 25/7/20
        if line == 'btn.10': return abs(self.Zblock[1:]-self.Zblock[:-1]) #[-1::-1]
        ### generic writing
        vtype = self.dictype[modName][line]
        nmedia = getNmedia(self) #OA 21/10/20
        if len(vtype)<nmedia: vtype=[vtype[0]]*nmedia #OA 21/10/20
        intp=[]
        for i in range(nmedia):
            if vtype[i]=='importArray': intp.append(4) 
            elif vtype[i] in ['one_value','zone'] : intp.append(3) 
            elif vtype[i]=='interpolate' : intp.append(1)  
            else: intp.append(0)

        value = block(self,modName,line,intp,iper=iper)
        if 'formula' in vtype : # case of a formula  
            for i in range(nmedia):
                f = self.dicformula[modName][line][i]
                if f != 'None':
                    value[i] = array(self.formExec(f)) ; # modif OA 3/10/18
        
        value = value.astype(numtype[3:]) # OA re-added line, it is necessary for int
        return value

    def formExec(self,s):
        s1 = s.replace('\n','\n\t')
        s1 = 'def yoyo(self):\n\t'+s1+'\n\treturn value'
        dct={}
        exec(s1,globals(),dct)
        b = types.MethodType(dct['yoyo'], self)
        return b()
    
    def runInterp(self,model,line,media,allOpt): # function added EV 19/02/20
        value,mess=zone2interp(self,model,line,media,allOpt,iper=0)
        if self.addin.mesh == None or self.getValueFromName(model,'MshType')==0: # OA 29/2/20 added mstType rect
            value=value[::-1] 
            nx,ny,x,y = getXYvects(self)
            extent = np.min(x), np.max(x), np.min(y), np.max(y)
        else : extent = None
        return value,mess,extent
    
    def save2array(self,model,line,media): # function added EV 19/02/20
        cols=self.dicval['Modflow']['dis.4'] # EV 30/04/20
        rows=self.dicval['Modflow']['dis.5']
        if line in 'dis.6': value=makeZblock(self)[:-1][media]
        elif line in 'dis.7': value=makeZblock(self)[-1:]
        else : value= self.getValueLong(model,line,media,iper=0)[media]
        value=value[::-1] ;#print('val1',value[0])
        return cols,rows,value
        
    def getUnits(self,modName,line,ik):   
        '''returns the units for a given line and keyword index'''
        s = ''
        tlist = ['-','sec','min','hours','days','years']
        llist = ['-','cm','m','km','ft']
        if modName in ['Modflow','Mt3dms']:
            tunit = tlist[self.dicval['Modflow']['dis.2'][4]]
            lunit = llist[self.dicval['Modflow']['dis.2'][5]]
        if modName[:4]=='Open':
            tunit = tlist[self.dicval['OpenFlow']['dis.3'][4]]
            lunit = llist[self.dicval['OpenFlow']['dis.3'][5]]            
        d0 = self.dickword[modName].lines[line]#;print ('d0',d0)
        if 'units' in d0: 
            s = d0['units'][ik];#print s
            s = s.replace('T',tunit)
            s = s.replace('L',lunit)
        return s
    
    def getPorosity(self,modName,line,iy,ix,iz): # OA 11/4/20 modified to consider Xsection
        ''' to get the porosity for a list of cell indices
        getValueLong returns data oriented in the x,y way (not modflow y)
        don't transform here -> NOT TRUE, transformation is needed here'''
        grd  = self.addin.getFullGrid()
        ncol, nrow = grd['nx'], grd['ny']
        val=self.getValueLong(modName,line,ik=0,iper=0)
        if self.addin.getDim() in ['Xsection','Radial']:
            return val[iy,0,ix] # OA 13/5/20 to get correct list directly
        else :
            return val[iz,nrow-iy-1,ix] #EV 20/04/20
    
    def getZcoord(self,iy,ix,iz): #EV 20/04/20
        ''' returns Z coordinates of the cell centers in 3D for a list of cell 
        indices  '''
        x2,y2,z2=getMesh3Dcenters(self)
        grd  = self.addin.getFullGrid()
        nrow = grd['ny']
        if self.addin.getDim() in ['Xsection','Radial']:
            iz=iy*1;nrow=1;iy=zeros(len(iz)).astype('int')
        return z2[iz,nrow-iy-1,ix]  
            
    def isObsFile(self,group,zname,esp):
        '''returns a list of bool stating if for the given zone and given species there
        is a file written by opf
        for species if solutes is 1 (first in the list) then all species have been written
        '''
        lFlow=['Head','Wcontent','Darcy V']
        lTrans=['Tracer','Temperature']
        lout=[]
        if 'Obspts' not in self.dicaddin.keys(): return [False]*len(esp)
        if len(self.dicaddin['Obspts'])<2: return [False]*len(esp)
        if len(self.dicaddin['Obspts'][1])==0 or self.dicaddin['Obspts'][0]==0: return [False]*len(esp)
        if zname in self.dicaddin['Obspts'][1]: 
            if group=='Flow':
                for e in esp: lout.append(e in lFlow)
            if group=='Transport':
                for e in esp: lout.append(e in lTrans)
            if group=='Chemistry':
                if 'Solutes' in self.dicaddin['Obspts'][4] : return [True]*len(esp)
            return lout
        else : return [False]*len(esp)
        
    def onPtObs(self,typ,iper,group,zname,esp,layers_in=0,ss=''): #EV 23/03/20
        """ get the values at observation zones, esp can be a list of species names
        typ[0]: B: breakthrough, P: profile, X : XYplot, V: vertical profile
        typ[1]: 0:head/conc, 1:weighted conc, 2:total discharge, 3:average flux
        iper is a period number
        group : Flow, Transport, Chemistry
        zname : zone name
        esp : list of species (F:Head,Wcontent,Flux;Tracer, C:species)
        layers_in : list of layer or 'all' for all layers of one zone
        ss : solute ss='' or sorbed species ss='S'
        """
        print('in core obspt', group,zname)
        ofile=self.isObsFile(group,zname,esp) # returns a list for the given variables
        #if self.addin.getDim() in ['Xsection','Radial']: layers_in='all'
    ### Get some parameters
        zlist=self.diczone['Observation'].dic['obs.1'] ## list of zone observation
        nx,ny,xvect,yvect = getXYvects(self) 
        grd = self.addin.getFullGrid()
        mtype = self.dicaddin['Model']['group'][:3]; ## model type, modflow or other
        modName='Modflow'
        if self.dicaddin['Model']['group'][:4]=='Open': modName='OpenFlow'
        MshType = self.getValueFromName(modName,'MshType') # OA 18/12/20
        ### Get a list of icol, irow, ilay
        if typ[0]=='X': ## XYplot
            ix=[];iy=[];typ='X0'
            #for i in zname: # zname is a list
               # ind = zlist['name'].index(zname)
               # x,y = list(zip(*zlist['coords'][ind]))
            for izon,xy in enumerate(zlist['coords']): # OA 18/12/20
                x,y = list(zip(*xy))
                if flgMesh==0: # OA 18/12/20
                    a,b,c = zone2index(self,x,y,x*1)
                    ix.append(a[0]);iy.append(b[0])
                else : # OA 18/12/20 fror Usg
                    ix = cellsUnderPolyOrd(self,zlist,0,izon);iy=[0]
            ix2=array(ix);iy2=array(iy);iz2=ix2*0.
        else:    ## TimeSerie & Profile
            izon = zlist['name'].index(zname)
            x,y = list(zip(*zlist['coords'][izon]))
            if MshType==0: # OA 18/12/20
                ix,iy,a,asin,acos = zone2index(self,x,y,x*1,'angle') ## ix,iy are orti indices (not modflow)
                if isclosed(self,x,y): # polygon
                    iy,ix = where(fillZone(nx,ny,ix,iy,a));
            else : # OA 18/12/20
                ix = cellsUnderPolyOrd(self,zlist,0,izon); # OA 8/6/21
                if len(ix) != 1 : ix = ix[0] # EV 21/6/21
                a = [0]*len(ix);iy=a*1;asin=a*1;acos=a*1 
                # should we add closed
            ix2,iy2,iz2,asin2,acos2=[],[],[],[],[]
            zlayers = media2layers(self,zlist['media'][izon]) # OA 19/3/20 OA added
    ### Get lists : ix2,iy2,iz2 are list of cell position in 3d
        try: 
            layers=[int(layers_in)]
            ix2,iy2,iz2,asin2,acos2=ix*1,iy*1,layers*len(ix),asin,acos
        except:
            if '-' in layers_in :
                    l1 = layers_in.split('-')
                    layers = list(range(int(l1[0]),int(l1[1])+1))
            elif ',' in layers_in:
                a = layers_in.split(',');layers = [int(x) for x in a]
                #print('lay :',layers)
            elif  layers_in == 'all' : layers = zlayers
            else : layers = layers_in # OA 13/5/20 layers->layers_in already a list
            if self.addin.getDim() in ['Xsection','Radial']:
                ix2,iy2,iz2,asin2,acos2=ix*1,iy*1,layers,asin,acos
            else :
                for il in layers: 
                    ix2.extend(ix);iy2.extend(iy);iz2.extend([il]*len(ix))
                    asin2.extend(asin);acos2.extend(acos)
        #if layers_in == 'all': layers=[-1]*len(zlayers) # OA 19/3/20 to make later the avergae on layers
        llay=iz2  #EV 23/03/20 
    ### Transform icol for modflow
        if mtype=='Mod' and MshType==0: iym = [ny-y-1 for y in iy2] # transform to modflow coords NOT for Usg #OA 6/11/21
        else : iym =iy2
        ix2,iym,iz2 = array(ix2),array(iym),array(iz2) # OA 11/4/20 added to index later
    ### Get list of time and period 
        t2 = self.getTlist2()
        if typ[0]=='B': iper=list(range(len(t2))) ## time-serie graph
        else : iper = [iper]
    ### Get the value
        pt=[] 
        '''pt is a list of tables (iper,irows) for each species/layers for the given zone'''
        labels=[''] 
       ## For Head and Wcontent (pt index is the layer)
        if group=='Flow': 
            m = self.flowReader.getPtObs(self,iym,ix2,iz2,iper,esp[0],ofile=ofile[0],zname=zname)
            if ofile[0]: tlong=m[:,0];m=m[:,1:2]
            if layers_in == 'all': # +below OA 11/4/2 to consider all
                pt.append(m)
                labels.append('all layers')
            else :
                for il in layers:
                    pt.append(m[:,iz2==il]); ## irow,icol,ilay #EV 23/03/20 
                    if typ[0]!='V':labels.append('lay'+str(il)) #EV 20/04/20
                if typ[0]=='V':
                    str1 = ','.join(str(l) for l in layers)
                    labels.append('lay'+str1)
       ## For Tracer
        elif group=='Transport':
            if mtype == 'Mod': 
                opt ='Mt3dms';ss=''
                m = self.transReader.getPtObs(self,iym,ix2,iz2,iper,opt,-1,esp[0],ss)
            else: 
                opt = esp[0];ss='' # OA 19/3/19
                m = self.transReader.getPtObs(self,iym,ix2,iz2,iper,opt,-1,esp[0],ss,ofile=ofile[0],zname=zname) #tracer -1
                if ofile[0]: tlong=m[:,0];m=m[:,1:2]
            if layers_in == 'all':  # +below OA 11/4/2 to consider all
                pt.append(m)
                labels.append('all layers')
            else :
                for il in layers:
                    pt.append(m[:,iz2==il]); ## irow,icol,ilay #EV 23/03/20 
                    if typ[0]!='V':labels.append('lay'+str(il)) #EV 20/04/20
                if typ[0]=='V':
                    str1 = ','.join(str(l) for l in layers)
                    labels.append('lay'+str1)
       ## For chemistry
        elif group=='Chemistry': 
            iesp,lesp = 0,self.addin.chem.getListSpecies();print(lesp)
            if mtype == 'Mod': opt ='Pht3d' # OA added 25/5
            else : opt = 'Chemistry'
            for i,e in enumerate(esp):
                if e in lesp: iesp = lesp.index(e) ;print('in core ptobs',e,iesp)
                m = self.transReader.getPtObs(self,iym,ix2,iz2,iper,opt,iesp,e,ss=ss,ofile=ofile[i],zname=zname)
                if ofile[i]: tlong=m[:,0];m=m[:,1:2]
                if layers_in == 'all':  # +below OA 11/4/2 to consider all
                    pt.append(m)
                    labels.append(str(e)+'_all layers')
                else :
                    for il in layers:
                        pt.append(m[:,iz2==il]);
                        if typ[0]!='V':labels.append(e+'_lay'+str(il)) #EV 20/04/20
                    if typ[0]=='V':
                        str1 = ','.join(str(l) for l in layers)
                        labels.append(str(e)+'_lay'+str1)
        else : pt=[1.] ## for flux
       ## For Darcy flux
        if esp[0]=='Flux': ## to get the flux, shall be a matrix (nper,nrow)
            #print('iper',iper)
            disx,disy = self.flowReader.getPtObs(self,iym,ix2,llay,iper,'flux'); ## provides the total flux from each cell [L3.T-1]  #EV 19/03/20 llay and not iz2
            disch = disx+disy # OA 19/3/20 modified
            thick = self.flowReader.getThicknessZone(self,iper,llay,ix2,iym) #EV 19/03/20 llay and not iz2
            dx,dy = array(grd['dx'])[ix2],array(grd['dy'])[iy2]
            if len(ix)==1: asin2 = 1; acos2 = 1 # OA 19/3/20 just one point
            f1,f2 = disx*asin2/dy/thick,disy*acos2/dx/thick ## this is the flux [L3.T-1.M-2]
            flux = sqrt(f1**2+f2**2) ## flux shall be a vector
            flux[flux<1e-12]=1e-12
         ## Por distributed
            if mtype == 'Mod' : 
                opt='Mt3dms' ; line='btn.11'
            else : opt='Min3pFlow' ; line='poro.1'
            lpor=self.getPorosity(opt,line,iym,ix2,llay) #print('por',lpor)
         ## Cells pore volume
            vol=thick*dx*dy*lpor #EV 23/03/20 
            #print('vol',vol)
            if layers_in == 'all': # OA 11/4/20 added below to have flux for all
                disch1,flux1,vol1 = [disch]*len(esp),[flux]*len(esp),[vol]*len(esp) #EV 20/04/20
            else :
                disch1,flux1,vol1 = [],[],[]
                for il in layers: 
                    disch1.append(disch[:,iz2==il])
                    flux1.append(flux[:,iz2==il])
                    vol1.append(vol[:,iz2==il])
                disch1=disch1*len(esp) #EV 20/04/20
                flux1=flux1*len(esp)
                vol1=vol1*len(esp)
            #print('vol1',vol1,'\n') ; print('disch1',disch1,'\n') ; print('flux1 1',(flux1),'\n')
        #print('pt',pt,'\n')
    ### Transform values for different type of graph and return them
        if esp[0] in ['Head','Wcontent']: typ=typ[0]+'0'
        #if esp[0] == 'Flux': pt=flux1 #EV 20/04/20 # OA 13/5/20 no it is 1 mutliplied later by flux or discharge
       ## Time-serie graph
        if typ[0]=='B': 
            labels[0]='time';
#            if ofile[0]:
#                t2,p1=pt[0][:,0],pt[0][:,1:2]
#                return t2,p1,labels
            #p1=[]; ## p1 : to make a table of (ntimes,nspecies)
            tlst2 = self.getTlist2()
            if len(ofile)>0:
                if ofile[0]: tlst2=tlong #*86400/self.dtu
            p1=zeros((len(tlst2),len(pt)))
            #tlst2=[]
            for i in range(len(pt)): # OA 11/4/20 modified flux to flux1 below
                #else: tlst2.append(t2)
                if typ[1]=='0': p1[:,i]=mean(pt[i],axis=1);## conc, pt[i] is a table (nper,nrow) 
                elif typ[1]=='1': p1[:,i]=sum(pt[i]*flux1[i],axis=1)/sum(flux,axis=1) ## weighted conc
                elif typ[1]=='2': p1[:,i]=sum(pt[i]*disch1[i],axis=1); ## total discharge [mol.T-1]
                elif typ[1]=='3': p1[:,i]=mean(pt[i]*flux1[i],axis=1); ## average flux [mol.T-1.M-2]
                elif typ[1]=='4': p1[:,i]=sum(pt[i]*vol1[i],axis=1); ## mass [mol] #EV 23/03/20 
            return tlst2,p1,labels
       ## Profile
        elif typ[0]=='P':  
            if MshType==0:
                xzon=xvect[ix];yzon=yvect[iy];#print(xzon,yzon) #EV 20/04/20 ix2 -> ix & iy2 -> iy
                d0=sqrt((xzon[1:]-xzon[:-1])**2.+(yzon[1:]-yzon[:-1])**2.)
            else :
                elc = self.addin.mesh.elcenters
                xc,yc = elc[ix,0],elc[ix,1]
                d0 = sqrt((xc[1:]-xc[:-1])**2.+(yc[1:]-yc[:-1])**2.)
            dist = concatenate((array(0.,ndmin=1),cumsum(d0)))
            p1=zeros((len(dist),len(pt)))
            labels[0]='distance'
            for i in range(len(pt)):
                if layers_in == 'all': # EV 30/04/20 
                    pt[i]=pt[i].reshape(len(layers),len(dist)).T
                    pt[i]=np.mean(pt[i],axis=1)
                if typ[1]=='0': p1[:,i]=pt[i]
                elif typ[1] in ['1','3']: p1[:,i]=pt[i]*flux1[i] ## weigthed conc= darcy flux #EV 20/04/20
                elif typ[1]=='2': p1[:,i]=pt[i]*disch1[i] ## discharge or mass discharge #EV 20/04/20
            return dist,p1,labels
        
        elif typ[0]=='V': #EV 20/04/20
            labels[0]='distance'
            distz=self.getZcoord(iym,ix2,llay) #print('distz',distz)
            #distz=list(dict.fromkeys(distz))
            p1=zeros((len(distz),len(pt)))
            for i in range(len(pt)): # i layer or species
#                if layers_in == 'all': # EV 30/04/20 
#                    pt[i]=pt[i].reshape(len(distz),len(pt)).T
#                    pt[i]=np.mean(pt[i],axis=1)
                if typ[1]=='0': p1[:,i]=pt[i][0] # there is jus tone species
                elif typ[1] in ['1','3']: p1[:,i]=pt[i]*flux1[i] ## weigthed conc= darcy flux #EV 20/04/20
                elif typ[1]=='2': p1[:,i]=pt[i]*disch1[i] ## discharge or mass discharge #EV 20/04/20
            return distz,p1,labels
       ## XY plot
        elif typ[0]=='X': 
            indCol = self.data['cols'].index(esp[0])
            labels = zlist['name']
            mes0=self.data['data'][:,indCol]; # OA 20/6 added three lines here to reorder and get labels
            npts = len(mes0)
            mes1,pt1,idx,pt0 = mes0*0,zeros((1,npts)),0,pt[0]
            for i in range(len(labels)): 
               if labels[i] in self.data['rows']:
                   mes1[idx] = mes0[self.data['rows'].index(labels[i])]
                   pt1[0,idx]=pt0[0,self.data['rows'].index(labels[i])]
                   idx+=1
            p1=zeros((len(mes1),1));p1[:,0]=pt1;
            return mes1,p1,'correl'
            
    def diffvect(self,v):
        v= v[1:]-v[:-1]
        return concatenate([v,v[-1:]])
            
    def diffvect(self,v):
        v= v[1:]-v[:-1]
        return concatenate([v,v[-1:]])
        
    def moyvect(self,v):
        v1= (v[2:]+v[1:-1]+v[:-2])/3
        return concatenate([v[:1],v1,v[-1:]])
        
class dicZone:
    """it is a dictionnary of zones, ordered by lines, each model has one dic"""
    def __init__(self,parent,modName):
        self.dic = {}
        groups = parent.dickword[modName].groups
        lines = parent.dickword[modName].lines;#print modName,groups,lines
        self.dicLines,self.dicLinesComm = {},{}
        for g in list(groups.keys()):
            self.dicLines[g] = []
            self.dicLinesComm[g] = [] ;#print g,groups[g]
            for ll in groups[g]:
                if lines[ll]['type'][0][:3] == 'arr':
                    self.dicLines[g].append(ll)
                    self.dicLinesComm[g].append(lines[ll]['comm'])#EV 25/10/2018 removed limit size of comm [:25])
            if len(self.dicLines[g])==0 : 
                self.dicLines.pop(g,None)
                self.dicLinesComm.pop(g,None)
        #print('core zones',modName,self.dicLines)
                
    def setDic(self,dic) :
        self.dic = dic           
    def getLinesDic(self):
        return self.dicLines
    def getLinesCommDic(self):
        return self.dicLinesComm    
        
    def getNbZones(self,line):
        if line in self.dic:
            return len(self.dic[line]['name'])
        else : return 0
        
    def setValueFromName(self,line,zname,val):
        for i,n in enumerate(self.dic[line]['name']):
            if n == zname : self.dic[line]['value'][i] = val
            
    def getIndexFromName(self,line,zname):
        for i,n in enumerate(self.dic[line]['name']):
            if n == zname : return i
            
    def getLineAndIndexFromName(self,zname):
        for line in list(self.dic.keys()):
            for i,n in enumerate(self.dic[line]['name']):
                if n == zname : return line,i
                               
    def getValue(self,line,parameter,nb):
        return self.dic[line][parameter][nb]
        
    def setValue(self,line,parameter,nb,value):
        self.dic[line][parameter][nb] = value

    def getMediaList(self,line,nb):
        md =  self.dic[line]['media'][nb]
        if type(md)==type([5,6]): return md
        else : return [md]
                   
    def addZone(self,line):
        if line not in list(self.dic.keys()):
            self.dic[line]={'number':[''],'name':[''],'coords':[''],'media':[''],'value':[''],'type':['']}
        else :
            for k in list(self.dic[line].keys()):
                self.dic[line][k].append('')
        self.dic[line]['number'][-1] = self.getNbZones(line)
            
    def delZone(self,line,iz):
        for k in list(self.dic[line].keys()):
            self.dic[line][k].pop(iz)
        nbz = self.getNbZones(line)
        self.dic[line]['number'] = list(range(nbz))
        if nbz==0: #EV 14/08/19
            del self.dic[line]
        
    def getTableOfZones(self,line):
        table = list(zip(self.dic[line]['name'],self.dic[line]['media'],self.dic[line]['value']))
        return table