# -*- coding: utf-8 -*-
"""
Created on Wed Sep  8 07:10:49 2021

@author: olivi
"""

class  OpF:
    def __init__(self):
        self.grpList=['FPRM','DIS','HorP','KHY','WEL','RCH','GHB','RIV','DRN','UNS','FSLV']
        self.groups={
        'FPRM':['fprm.1'],
        'DIS':['dis.'+str(i) for i in range(1,9)],
        'HorP':['head.1','head.2','press.1','press.2'],
        'KHY':['khy.'+str(i) for i in range(1,4)],
        'WEL':['wel'],
        'RCH':['rch'],
        'GHB':['ghb'],
        'RIV':['riv'],
        'DRN':['drn'],
        'UNS':['uns.'+str(i) for i in range(1,6)], # unsaturated
        'FSLV':['fslv.'+str(i) for i in range(1,4)] # solver
        }
        self.lines={
        'fprm.1':{'comm':'Flow parameters','cond':'','kw':['OFPMS','OFPMEQ'],
                'detail':[['type of flow','saturated','unsaturated','2phase'],
                    ['pressure equilibrium','no','yes']],
                'type':['choice','choice'],
                'default':[0,0]},
        'dis.1':{'comm':'Type of mesh','cond':'','kw':['MshType'],
                'detail':[['','Rectangle','Nested','Triangle','Voronoi']],
                'type':['choice'],'default':[0]},
        'dis.2':{'comm':'Mesh','cond':'','kw':['O_Mesh'],'detail':[],
                  'type':['arrfloat'],'default':[0.]},
        'dis.3':{'comm':'Model properties','cond':'','kw':['ONLAY','ONROW','ONCOL','ONPER','OTUNI','OLUNI'],
               'detail':['Nb of layers','Nb of rows','Nb of columns','Nb of periods',['time units','-','sec','min','hours','days','years'],
                ['length units','-','ft','m','cm']],'type':['int','int','int','int','choice','choice'],
                'default':[1,10,10,1,4,2]},
        'dis.4':{'comm':'Col width','cond':'','kw':['ODELR(ONCOL)'],'detail':['Col width'],
                        'type':['vecfloat'],'default':[10],'units':['L']},
        'dis.5':{'comm':'Row height','cond':'','kw':['ODELC(ONROW)'],'detail':['Row height'],
                        'type':['vecfloat'],'default':[10],'units':['L']},
        'dis.6':{'comm':'Top','cond':'','kw':['OTOP(ONCELL)'],'detail':[],
                        'type':['arrfloat'],'default':[10.],'units':['L']},
        'dis.7':{'comm':'Bottom','cond':'','kw':['OBOTM(ONLAY,ONCELL)'],'detail':[],
                        'type':['arrfloat'],'default':[0.],'units':['L']},
        'dis.8':{'comm':'Periods characteristics','cond':'','kw':['PERLENu','NSTPu','TSMULTu', 'SsTru'], # OA 3/11/18 
                'detail':['Period length','internal steps','multiplier',['type','Steady state','transient']],
                'type':['float','int','float','choice'],
                'default':[1.,1,1.,0],'units':['T','','','']},
        'head.1':{'comm':'Initial heads','cond':'OFPMS<2','kw':['OIH'],
                 'detail':[],'type':['arrfloat'],'default':[10.]},
        'head.2':{'comm':'Fixed heads','cond':'OFPMS<2','kw':['OFH'],
                 'detail':[],'type':['arrfloat'],'default':[10.]},
        'press.1':{'comm':'Initial pressure','cond':'OFPMS==2','kw':['OIP'],
                 'detail':[],'type':['arrfloat'],'default':[0.]},
        'press.2':{'comm':'Fixed pressure','cond':'OFPMS==2','kw':['OFPR'],
                 'detail':[],'type':['arrfloat'],'default':[0.]},
        'khy.1':{'comm':'Hydraulic cond. parms','cond':'','kw':['OKTYP'],
                 'detail':[['','Kv type','ratio','value']],
                'type':['choice'],'default':[0]},
        'khy.2':{'comm':'Hor hydraulic cond.','cond':'','kw':['OKH'],
                 'detail':[],'type':['arrfloat'],'default':[10.]},
        'khy.3':{'comm':'Vert K ratio/value','cond':'','kw':['OKV'],
                 'detail':[],'type':['arrfloat'],'default':[1.]},
        'wel':{'comm':'Wells','cond':'','kw':['OWELL'],
                 'detail':[],'type':['arrfloat'],'default':[0.],'units':'m3/d'},
        'rch':{'comm':'Recharge','cond':'','kw':['ORCH'],
                 'detail':[],'type':['arrfloat'],'default':[0.]},
        'ghb':{'comm':'Gen Head Bound','cond':'','kw':['OGHB'],
                'names':['Elevation','Conductance'],
                 'detail':[],'type':['arrfloat'],'default':[0.]},
        'riv':{'comm':'River','cond':'','kw':['ORIV'],
                 'names':['Stage','Conductance','Elevation'],
                 'detail':[],'type':['arrfloat'],'default':[0.]},
        'drn':{'comm':'Drain','cond':'','kw':['ODRN'],
                 'names':['Elevation','Conductance'],
                 'detail':[],'type':['arrfloat'],'default':[0.]},
        'uns.1':{'comm':'Unsat parameters','cond':'','kw':['OUNPR','OSWMAX'],
                 'detail':[['capillary type','van Genuchten'],'Sw_max'],
                 'type':['choice','float'],'default':[0,0.999]},
        'uns.2':{'comm':'sw_min','cond':'','kw':['OSWMIN'],
                 'detail':[],'type':['arrfloat'],'default':[0.1]},
        'uns.3':{'comm':'alpha_vg','cond':'','kw':['OAVG'],
                 'detail':[],'type':['arrfloat'],'default':[1.]},
        'uns.4':{'comm':'n_vg','cond':'','kw':['ONVG'],
                 'detail':[],'type':['arrfloat'],'default':[1.5]},
        'uns.5':{'comm':'sw_init','cond':'','kw':['OSWINI'],
                 'detail':[],'type':['arrfloat'],'default':[1.]},
        'fslv.1':{'comm':'solver for h','cond':'','kw':['SOFS','SOFPRE','SOFTOL','SOFRTOL','SOFSIMP'],
                 'detail':[['solver','PBiCGStab'],['preconditioner','DIC','DILU'],
                           'tolerance','relTol','SteadyTol'],
                 'type':['choice','choice','float','float','float'],
                 'default':[0,0,1e-12,0,1e-6]},
        'fslv.2':{'comm':'Picard','cond':'','kw':['SOFPI'],
                 'detail':['tolerance','minIter','maxIter','nIterStability'],
                 'type':['float','int','int','int'],'default':[0.01,3,10,5]},
        'fslv.3':{'comm':'Timing','cond':'','kw':['SODT0','SOMXDT','SOMXCO'],
                'detail':['deltaT0','maxDeltaT','maxCourant'],
                'type':['float','float','float'],
                'default':[1,100,0.75]}
        }

#Transport
class  OpT:
    def __init__(self):
        self.grpList=['TPRM','PORO','DSP','DIFFU','CONC','TEMP','RCT','TSLV']
        self.groups={
        'TPRM':['tprm'],
        'PORO':['poro'],
        'DSP':['dsp'],
        'DIFFU':['diffu.1','diffu.2','diffu.3'],
        'CONC':['cactiv','cinit','cfix','cghb','cwel','crch'],
        'TEMP':['temp.1','tcps','tlbds','tactiv','tinit','tfix','tghb','twel','trch'],
        'RCT':['rct.1','rct.2a','rct.2b','rct.2c','rct.3','rct.4','rct.5','rct.6'],
        'TSLV':['tschm','tslv'],
        }
        self.lines={
        'tprm':{'comm':'Transp parameters','cond':'','kw':['OTTYP','OTSTDY'],
                'detail':[['transp type','Concentration','Temperature','both'],
                        ['transp. steady','no','yes']],
                'type':['choice','choice'],'default':[0,0]},
        'poro':{'comm':'Porosity','cond':'','kw':['EPS'],
                'detail':[],'type':['arrfloat'],'default':[0.3]},
        'dsp':{'comm':'Dispersion','cond':'','kw':['ALPHL','ALPHT'],
                'detail':[],'type':['float','float'],'default':[1.,0.1]},
        'diffu.1':{'comm':'Diffusion model','cond':'','kw':['ODFMOD'],
                'detail':[['model','millington']],
                'type':['choice'],'default':[0]},
        'diffu.2':{'comm':'Diffusion model coeffs','cond':'','kw':['ODFCOF'],
                'detail':[],
                'type':['float'],'default':[3.333]},
        'diffu.3':{'comm':'Diffusion','cond':'','kw':['ODFFW','ODFFG'],
                'detail':[],'type':['float','float'],'default':[1e-10,1e-6],
                'units':['m2/s','m2/s']},
                   
        'cactiv':{'comm':'active zone','cond':'','kw':['OCACT'],
                'detail':[],'type':['arrfloat'],'default':[0.]},
        'cinit':{'comm':'Initial conc','cond':'','kw':['OICONC'],
                'detail':[],'type':['arrfloat'],'default':[0.]},
        'cfix':{'comm':'Fixed conc','cond':'','kw':['OFCONC'],
                'detail':[],'type':['arrfloat'],'default':[0.]},
        'cghb':{'comm':'Ghb conc','cond':'','kw':['OGCONC'],
                'detail':[],'type':['arrfloat'],'default':[0.]},
        'cwel':{'comm':'well conc','cond':'','kw':['OWCONC'],
                'detail':[],'type':['arrfloat'],'default':[0.]},
        'crch':{'comm':'Recharge conc','cond':'','kw':['ORCONC'],
                'detail':[],'type':['arrfloat'],'default':[0.]},
                
        'temp.1':{'comm':'Temp parameters','cond':'','kw':['OTCPW:','OTLBDW'],
                'detail':['water heat capacity','water thermal conductivity'],
                'type':['float','float'],'default':[4182,0.6]},
        'tcps':{'comm':'solid heat capacity','cond':'','kw':['OTCPS'],
                'detail':[],'type':['arrfloat'],'default':[840],'units':'J.kg-1.C-1'},
        'tlbds':{'comm':'solid thermal cond','cond':'','kw':['OTLBDS'],
                'detail':[],'type':['arrfloat'],'default':[3.5],'units':'J.s-1.m-1.C-1'},
        'tactiv':{'comm':'active zone for temp','cond':'','kw':['OTACT'],
                'detail':[],'type':['arrfloat'],'default':[0.]},
        'tinit':{'comm':'Initial temp','cond':'','kw':['OITMP'],
                'detail':[],'type':['arrfloat'],'default':[0.]},
        'tfix':{'comm':'Fixed temp','cond':'','kw':['OFTMP'],
                'detail':[],'type':['arrfloat'],'default':[0.]},
        'tghb':{'comm':'Ghb temp','cond':'','kw':['OGTMP'],
                'detail':[],'type':['arrfloat'],'default':[0.]},
        'twel':{'comm':'well temp','cond':'','kw':['OWTMP'],
                'detail':[],'type':['arrfloat'],'default':[0.]},
        'trch':{'comm':'Recharge temp','cond':'','kw':['ORTMP'],
                'detail':[],'type':['arrfloat'],'default':[0.]},
                
        'tschm':{'comm':'scheme for C or T','cond':'','kw':['OTSCH'],
                 'detail':[['scheme','Gauss vanLeer','Gauss limitedLinear01 1','Gauss SuperBee']],
                 'type':['choice'],'default':[0]},
        #RCT
        'rct.1':{'comm':'major flags','cond':'','kw':['OISOTH','OIREAC','OIGETSC'],
                 'detail':[['sorption','no sorption','linear'],
                        ['kinetic reaction','no reaction','1st order'],
                        ['initial conc.','at equilibrium','given']],
                'type':['choice','choice','choice'],
                'default':[0,0,2,0]},
        #removed 'freundlich','langmuir','kinetic sorption','dual domain no sorp','dual domain with sorp'],#,'0th order'],
        'rct.2a':{'comm':'density','cond':'OISOTH in [1,2,3,4,6]','kw':['ORHOB(NLAY,NCELL)'],
                 'detail':['bulk density'],'type':['arrfloat'],'default':[1.8]},
        'rct.2b':{'comm':'poro immob domain','cond':'OISOTH>4','kw':['OPRSITY2(NLAY,NCELL)'],
                 'detail':['porosity'],'type':['arrfloat'],'default':[1.]},
        'rct.2c':{'comm':'concentrations','cond':'OIGETSC>0','kw':['ORCONC(NLAY,NCELL)'],
                 'detail':[''],'type':['arrfloat'],'default':[0.]},
        'rct.3':{'comm':'1st sorption parm','cond':'OISOTH>0','kw':['OSP1(NLAY,NCELL)'],
                 'detail':[''],'type':['arrfloat'],'default':[1.]},
        'rct.4':{'comm':'2nd sorption parm','cond':'OISOTH>0','kw':['OSP2(NLAY,NCELL)'],
                 'detail':[''],'type':['arrfloat'],'default':[1.]},
        'rct.5':{'comm':'1st order rate','cond':'OIREAC>0','kw':['ORC1'],
                 'detail':[''],'type':['float'],'default':[0.]},
        'rct.6':{'comm':'1st rate sorbed','cond':'OIREAC>0','kw':['ORC2'],
                 'detail':[''],'type':['float'],'default':[0.]},
        #Solver
        'tslv':{'comm':'solver for C','cond':'','kw':['OSTSO','OSTPRE','OSTTOL','OSTRTOL','OSTDCMX'],
                 'detail':[['solver','PBiCG'],['preconditioner','DIC','DILU'],
                           'tolerance','relTol','dCmax'],
                 'type':['choice','choice','float','float','float'],
                 'default':[0,1,1e-12,0,0.02]},
        }
    
#Chemistry
class  OpC:
    def __init__(self):
        self.grpList=['CHPRM','SOLU','GAS','CHSLV']
        self.groups={
        'CHPRM':['chprm'],
        'SOLU':['sactiv','sinit','sfix','swel','srch','sghb','sriv'],
        'GAS':['ginit','gfix','gwel'],
        'CHSLV':['chslv'],
        }
        self.lines={
        'chprm':{'comm':'Chem parameters','cond':'','kw':['OCHNSTP'],
                'detail':['reaction stp interval'],
                'type':['int'],'default':[10]},
        'sactiv':{'comm':'active zone','cond':'','kw':['OSACT'],
                'detail':[],'type':['arrfloat'],'default':[0.]},
        'sinit':{'comm':'Initial solutions','cond':'','kw':['OISOL'],
                'detail':[],'type':['arrfloat'],'default':[0.]},
        'sfix':{'comm':'Sources','cond':'','kw':['OSSOL'],
                'detail':[],'type':['arrfloat'],'default':[0.]},
        'swel':{'comm':'well solutions','cond':'','kw':['OWSOL'],
                'detail':[],'type':['arrfloat'],'default':[0.]},
        'sghb':{'comm':'ghb solutions','cond':'','kw':['OWSOL'],
                'detail':[],'type':['arrfloat'],'default':[0.]},
        'sriv':{'comm':'river solutions','cond':'','kw':['OWSOL'],
                'detail':[],'type':['arrfloat'],'default':[0.]},
        'srch':{'comm':'Recharge solu','cond':'','kw':['ORSOL'],
                'detail':[],'type':['arrfloat'],'default':[0.]},
        'ginit':{'comm':'Initial gas','cond':'','kw':['OIGAS'],
                'detail':[],'type':['arrfloat'],'default':[0.]},
        'gfix':{'comm':'fixed gas','cond':'','kw':['OSGAS'],
                'detail':[],'type':['arrfloat'],'default':[0.]},
        'gwel':{'comm':'well gas','cond':'','kw':['OWGAS'],
                'detail':[],'type':['arrfloat'],'default':[0.]},
        'chslv':{'comm':'solver for Cwi','cond':'','kw':['OCHSO'],
                 'detail':[['solver','PBiCG'],['preconditioner','DIC','DILU'],
                           'tolerance','relTol'],
                 'type':['choice','choice','float','float'],'default':[0,1,1e-12,0]},
        }
