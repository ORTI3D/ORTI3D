# -*- coding: utf-8 -*-
"""
Created on Wed Sep  8 07:10:49 2021

@author: olivi
"""

class  OpF:
    def __init__(self):
        self.grpList=['FPRM','DIS','HEAD','KHY','WEL','RCH','GHB','RIV','DRN','UNS','FSLV']
        self.groups={
        'FPRM':['fprm.1'],
        'DIS':['dis.'+str(i) for i in range(1,6)],
        'HEAD':['head.1','head.2'],
        'KHY':['khy.'+str(i) for i in range(1,4)],
        'WEL':['wel.1'],
        'RCH':['rch.1'],
        'GHB':['ghb.1'],
        'RIV':['riv.1'],
        'DRN':['drn.1'],
        'UNS':['uns.1','uns.2','uns.3'], # unsaturated
        'FSLV':['fslv.1'], # solver
        }
        self.lines={
        'fprm.1':{'comm':'Flow parameters','cond':'','kw':['OFP1'],
                'detail':['parm1'],'type':['int'],'default':[0]},
        'dis.1':{'comm':'Type of mesh','cond':'','kw':['MshType'],
                'detail':[['','Rectangle','Nested','Triangle','Voronoi']],
                'type':['choice'],'default':[0]},
        'dis.2':{'comm':'Mesh','cond':'','kw':['O_Mesh'],'detail':[],
                  'type':['arrfloat'],'default':[0.]},
        'dis.3':{'comm':'Top','cond':'','kw':['TOP(NCELL)'],'detail':[],
                        'type':['arrfloat'],'default':[10.],'units':['L']},
        'dis.4':{'comm':'Bottom','cond':'','kw':['BOTM(NLAY,NCELL)'],'detail':[],
                        'type':['arrfloat'],'default':[0.],'units':['L']},
        'dis.5':{'comm':'Periods characteristics','cond':'','kw':['PERLENu','NSTPu','TSMULTu', 'SsTru'], # OA 3/11/18 
                'detail':['Period length','internal steps','multiplier',['type','Steady state','transient']],
                'type':['float','int','float','choice'],
                'default':[1.,1,1.,0],'units':['T','','','']},
        'head.1':{'comm':'Initial heads','cond':'','kw':['OIH'],
                 'detail':[],'type':['arrfloat'],'default':[10.]},
        'head.2':{'comm':'Fixed heads','cond':'','kw':['OFH'],
                 'detail':[],'type':['arrfloat'],'default':[10.]},
        'khy.1':{'comm':'Hydraulic cond. parms','cond':'','kw':['OKTYP'],
                 'detail':['Kv type','ratio','value'],
                'type':['choice'],'default':[0]},
        'khy.2':{'comm':'Hor hydraulic cond.','cond':'','kw':['OKH'],
                 'detail':[],'type':['arrfloat'],'default':[10.]},
        'khy.3':{'comm':'Vert hydraulic cond.','cond':'','kw':['OKV'],
                 'detail':[],'type':['arrfloat'],'default':[1.]},
        'wel.1':{'comm':'Wells','cond':'','kw':['OWELL'],
                 'detail':[],'type':['arrfloat'],'default':[0.]},
        'rch.1':{'comm':'Recharge','cond':'','kw':['ORCH'],
                 'detail':[],'type':['arrfloat'],'default':[0.]},
        'ghb.1':{'comm':'Gen Head Bound','cond':'','kw':['OGHB'],
                 'detail':[],'type':['arrfloat'],'default':[0.]},
        'riv.1':{'comm':'River','cond':'','kw':['ORIV'],
                 'detail':[],'type':['arrfloat'],'default':[0.]},
        'drn.1':{'comm':'Drain','cond':'','kw':['ODRN'],
                 'detail':[],'type':['arrfloat'],'default':[0.]},
        'uns.1':{'comm':'Unsat parameters','cond':'','kw':['OUNPR'],
                 'detail':[],'type':['arrfloat'],'default':[0.]},
        'uns.2':{'comm':'alpha_vg','cond':'','kw':['OAVG'],
                 'detail':[],'type':['arrfloat'],'default':[1.]},
        'uns.3':{'comm':'n_vg','cond':'','kw':['ONVG'],
                 'detail':[],'type':['arrfloat'],'default':[1.]},
        'fslv.1':{'comm':'solver for h','cond':'','kw':['SOFH'],
                 'detail':[['solver','PBiCGStab'],['preconditioner','DIC','DILU'],
                           'tolerance','relTol'],
                 'type':['choice'],'default':[0,0,1e-12,0]},
        'fslv.2':{'comm':'Picard','cond':'','kw':['SOFPI'],
                 'detail':['tolerance','minIter','maxIter','nIterStability'],
                 'type':['float','int','int','int'],'default':[0.01,3,10,5]},
                  
        }

#Transport
class  OpT:
    def __init__(self):
        self.grpList=['TPRM','PORO','DSP','CONC','TSLV']
        self.groups={
        'TPRM':['tprm.1'],
        'PORO':['poro.1'],
        'DSP':['dsp.1'],
        'CONC':['conc.1','conc.2','conc.3'],
        'TSLV':['tslv.1'],
        }
        self.lines={
        'tprm.1':{'comm':'Transp parameters','cond':'','kw':['OTP1'],
                'detail':['parm1'],'type':['int'],'default':[0]},
        'poro.1':{'comm':'Porosity','cond':'','kw':['EPS'],
                'detail':[],'type':['arrfloat'],'default':[0.3]},
        'dsp.1':{'comm':'Dispersion','cond':'','kw':['ALPHL','ALPHT'],
                'detail':[],'type':['float','float'],'default':[1.,0.1]},
        'conc.1':{'comm':'active zone','cond':'','kw':['OCACT'],
                'detail':[],'type':['arrfloat'],'default':[0.]},
        'conc.2':{'comm':'Initial conc','cond':'','kw':['OICONC'],
                'detail':[],'type':['arrfloat'],'default':[0.]},
        'conc.3':{'comm':'Fixed conc','cond':'','kw':['OFCONC'],
                'detail':[],'type':['arrfloat'],'default':[0.]},
        'tslv.1':{'comm':'solver for Cw','cond':'','kw':['SOTC'],
                 'detail':[['solver','PBiCG'],['preconditioner','DIC','DILU'],
                           'tolerance','relTol'],
                 'type':['choice'],'default':[0,1,1e-12,0]},
        }
    
#Chemistry
class  OpC:
    def __init__(self):
        pass