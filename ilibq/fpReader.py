# -*- coding: utf-8 -*-
"""
Created on Sun Sep 21 08:25:52 2014

@author: olive
"""
import os
from array import array as arr2


class fpReader:
    
    def __init__(self, fDir, fName):
        self.fDir,self.fName = fDir,fName
    
    def readHeadFile(self,core,iper=0):
        return self.readFile(core,'fhd',iper)
        
    def readUCN(self,core,typ,iper,iesp,species):
        '''typ is to differentiate between transport and chem (like mt3d pht3d)
        '''
        return self.readFile(core,'ucn',iper)
        
    def readFile(self,core,typ,iper):
        f1=open(self.fDir+os.sep+self.fName+'.'+typ,'rb')
        ncells = core.addin.fipy.ncells
        f1.seek(iper*ncells*4)
        data=arr2('f');data.fromfile(f1,ncells)
        f1.close()
        return data
