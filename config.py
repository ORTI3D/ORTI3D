# -*- coding: utf-8 -*-
"""
Created on Sun Dec 15 20:53:49 2013

@author: olive
"""
sclist = ['amax','amin','arange','arctan','argmin','argsort','around','array','c_','ceil','clip','compress',
    'concatenate','cos','cumsum','dot','equal','exp','floor','linspace','loadtxt','log','log10','logspace','maximum','mean','median',
    'meshgrid','minimum','mod','nonzero','ones','put','putmask','r_','rand','ravel','reshape','savetxt','shape',
    'sign','sin','sort','sqrt','sum','take','transpose','unique','where','zeros','zeros_like','prod','split',
    'delete']
for n in sclist: exec('from scipy import '+n)

import matplotlib.tri as mptri
from numpy.linalg import solve
import warnings,os
warnings.filterwarnings("ignore")

def nice(x):
    """makes a flot in a nice format"""
    if x!=0:om=int(round(log10(abs(x)))) # EV 27/08/19
    else:om=1 
    dec = max(5-om,0)
    if abs(om)<4: c = ' %+10.'+str(dec)+'f'
    else : c = '%.4e'
    return (c%x).strip()

class Config():
    def __init__(self,core):
        if core.gui != None :
#            if core.gui.gtyp=='wx':
#                import wxDialogs,wxShow,wxValueDialog
#                self.dialogs,self.show,self.valDlg = wxDialogs,wxShow,wxValueDialog
#                self.gtyp = 'wx'
#            elif core.gui.gtyp in ['qt','qgis']:
            import qtDialogs,qtValueDialog
            self.dialogs,self.valDlg = qtDialogs,qtValueDialog
            self.gtyp = core.gui.gtyp
            lfi=os.listdir(core.baseDir) # oa 29/3/17
            if 'ilibq' in lfi: self.typInstall = 'python' # oa 29/3/17
            else : self.typInstall = 'exe' # oa 29/3/17
            self.curVar={}
            

