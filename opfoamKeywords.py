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
        'HorP':['head.1','head.2','press.1','press.2','press.3'],
        'KHY':['khy.'+str(i) for i in range(1,6)],
        'WEL':['wel'],
        'RCH':['rch'],
        'GHB':['ghb'],
        'RIV':['riv'],
        'DRN':['drn'],
        'UNS':['uns.'+str(i) for i in range(1,6)], # unsaturated
        'FSLV':['fslv.'+str(i) for i in range(1,4)] # solver
        }
        self.lines={
        'fprm.1':{'comm':'Flow parameters','cond':'','kw':['OFPMS','OFPMEQ','OFXSECT','OFRHOW','OFMUW','OFRHOG','OFMUG'],
                'detail':[['type of flow','saturated','unsaturated','2phases'],
                    ['pressure equilibrium','no','yes'],'Xsection',
                    'rhow (kg/m3)','muw (kg/m/s)','rhog (kg/m3)','mug (kg/m/s)'],
                'type':['choice','choice','int','float','float','float','float'],
                'default':[0,0,0,1e3,1e-3,3.5,1e-5]},
        'dis.1':{'comm':'Type of mesh','cond':'','kw':['MshType'],
                'detail':[['','Rectangle','Nested','Triangle','Voronoi']],
                'type':['choice'],'default':[0]},
        'dis.2':{'comm':'Mesh','cond':'','kw':['O_Mesh'],'detail':[],
                  'type':['arrfloat'],'default':[0.]},
        'dis.3':{'comm':'Model properties','cond':'','kw':['ONLAY','ONROW','ONCOL','ONPER','OTUNI','OLUNI'],
               'detail':['Nb of layers','Nb of rows','Nb of columns','Nb of periods',
                ['time units','-','sec','min','hours','days','years'],
                ['length units','-','cm','m','km','ft']],'type':['int','int','int','int','choice','choice'],
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
                 'detail':[],'type':['arrfloat'],'default':[10.],'units':['L']},
        'head.2':{'comm':'Fixed heads','cond':'OFPMS<2','kw':['OFH'],
                 'detail':[],'type':['arrfloat'],'default':[10.],'units':['L']},
        'press.1':{'comm':'Initial pressure','cond':'OFPMS==2','kw':['OIP'],
                 'detail':[],'type':['arrfloat'],'default':[0.]},
        'press.2':{'comm':'Fixed pressure','cond':'OFPMS==2','kw':['OFPR'],
                 'detail':[],'type':['arrfloat'],'default':[0.]},
        'press.3':{'comm':'gravity norm','cond':'OFPMS==2','kw':['OFPR'],
                 'detail':[],'type':['float'],'default':[9.81]},
        'khy.1':{'comm':'Hydraulic cond. parms','cond':'','kw':['OKTYP'],
                 'detail':[['','Kh/Kv ratio','value']],
                'type':['choice'],'default':[0]},
        'khy.2':{'comm':'Hor hydraulic cond.','cond':'OFXSECT==0 or ONCOL>1','kw':['OKH'],
                 'detail':[],'type':['arrfloat'],'default':[10],
                 'units':['L/T']},
        'khy.3':{'comm':'Vert K or ratio/value','cond':'','kw':['OKV'],
                 'detail':[],'type':['arrfloat'],'units':['L/T'],'default':[1.]},
        'khy.4':{'comm':'Storage confined.','cond':'','kw':['OSTOC'],
                 'detail':[],'type':['arrfloat'],'default':[1e-4],
                 'units':['-']},
        'khy.5':{'comm':'Storage unconfined','cond':'','kw':['OSTOU'],
                 'detail':[],'type':['arrfloat'],'default':[0.25],
                 'units':['-']},
        'wel':{'comm':'Wells','cond':'','kw':['OWELL'],
                 'detail':[],'type':['arrfloat'],'default':[0.],'units':['L3/T']},
        'rch':{'comm':'Recharge','cond':'','kw':['ORCH'],
                 'detail':[],'type':['arrfloat'],'default':[0.],'units':['L/T']},
        'ghb':{'comm':'Gen Head Bound','cond':'','kw':['OGHB'],
                'names':['Elevation','Conductance'],
                 'detail':[],'type':['arrfloat'],'default':[0.]},
        'riv':{'comm':'River','cond':'','kw':['ORIV'],
                 'names':['Stage','Conductance','Elevation'],
                 'detail':[],'type':['arrfloat'],'default':[0.]},
        'drn':{'comm':'Drain','cond':'','kw':['ODRN'],
                 'names':['Elevation','Conductance'],
                 'detail':[],'type':['arrfloat'],'default':[0.]},
        'uns.1':{'comm':'Unsat parameters','cond':'','kw':['OUNPR','OSWMAX','OFCAPL'],
                 'detail':[['capillary type','van Genuchten'],'Sw_max',
                           ['Capill. active','no','yes']],
                 'type':['choice','float','choice'],'default':[0,0.999,1]},
        'uns.2':{'comm':'sw_min','cond':'','kw':['OSWMIN'],
                 'detail':[],'type':['arrfloat'],'default':[0.1]},
        'uns.3':{'comm':'alpha_vg','cond':'','kw':['OAVG'],
                 'detail':[],'type':['arrfloat'],'default':[1.],
                 'units':['L-1']},
        'uns.4':{'comm':'n_vg','cond':'','kw':['ONVG'],
                 'detail':[],'type':['arrfloat'],'default':[1.5]},
        'uns.5':{'comm':'sw_init','cond':'','kw':['OSWINI'],
                 'detail':[],'type':['arrfloat'],'default':[1.]},
        'fslv.1':{'comm':'solver for h','cond':'','kw':['SOFS','SOFPRE','SOFTOL','SOFRTOL','SOFSIMP'],
                 'detail':[['solver','PBiCGStab'],['preconditioner','DIC','DILU'],
                           'tolerance','relTol','SteadyTol'],
                 'type':['choice','choice','float','float','float'],
                 'default':[0,0,1e-12,0,1e-6]},
        'fslv.2':{'comm':'Picard','cond':'','kw':['SOFPI1','SOFPI2','SOFPI3'],
                 'detail':['tolerance','minIter','maxIter','nIterStability'],
                 'type':['float','int','int','int'],'default':[0.01,1,20,5]},
        'fslv.3':{'comm':'Timing','cond':'','kw':['SODT0','SOMXDT','SOMXCO'],
                'detail':['deltaT0','maxDeltaT','maxCourant'],
                'type':['float','float','float'],
                'default':[1e-3,1,0.75]}
        }

#Transport
class  OpT:
    def __init__(self):
        self.grpList=['TPRM','PORO','DSP','DIFFU','CONC','TEMP','RCT','TSLV']
        self.groups={
        'TPRM':['tprm'],
        'PORO':['poro'],
        'DSP':['dsp'],
        'DIFFU':['diffu.1','diffu.2'],
        'CONC':['cactiv','cinit','cfix','cghb','cwel','crch'],
        'TEMP':['temp.1','tcps','tlbds','tactiv','tinit','tfix','tghb','twel','trch'],
        'RCT':['rct.1','rct.2a','rct.2b','rct.2c','rct.3','rct.4','rct.5'],
        'TSLV':['tschm','tslv'],
        }
        self.lines={
        'tprm':{'comm':'Transp parameters','cond':'',
                'kw':['OTTYP','OTSTDY','OTDET1','OTDET2'],
                'detail':[['transp type','Concentration','Temperature','both'],
                        ['transp. steady','no','yes'],
                        ['De vary with T','no','linear'],['De f(T) fact.']],
                'type':['choice','choice','choice','Textlong'],'default':[0,0,0,['']]},
        'poro':{'comm':'Porosity','cond':'','kw':['EPS'],
                'detail':[],'type':['arrfloat'],'default':[0.3]},
        'dsp':{'comm':'Dispersion','cond':'','kw':['ALPHL','ALPHT'],
                'detail':[],'type':['float','float'],'default':[1.,0.1]},
        'diffu.1':{'comm':'Diffusion model','cond':'','kw':['ODFMOD'],
                'detail':[['model','none','millington10','millington7']],
                'type':['choice'],'default':[1]},
        'diffu.2':{'comm':'D0 values','cond':'','kw':['ODFFW','ODFFG'],
                'detail':['Dw m2/s','Dg m2/s'],'type':['float','float'],'default':[1e-10,1e-6],
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
                
        'temp.1':{'comm':'Temp parameters','cond':'','kw':['OTCPW','OTLBDW'],
                'detail':['water thermal conductivity','water heat capacity'],
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
                'default':[0,0,2]},
        #removed 'freundlich','langmuir','kinetic sorption','dual domain no sorp','dual domain with sorp'],#,'0th order'],
        'rct.2a':{'comm':'density','cond':'OISOTH in [1,2,3,4,6]','kw':['ORHOB(NLAY,NCELL)'],
                 'detail':['bulk density'],'type':['arrfloat'],'default':[1.8]},
        'rct.2b':{'comm':'poro immob domain','cond':'OISOTH>4','kw':['OPRSITY2(NLAY,NCELL)'],
                 'detail':['porosity'],'type':['arrfloat'],'default':[1.]},
        'rct.2c':{'comm':'concentrations','cond':'OIGETSC>0','kw':['ORCONC(NLAY,NCELL)'],
                 'detail':[''],'type':['arrfloat'],'default':[0.]},
        'rct.3':{'comm':'foc','cond':'OISOTH>0','kw':['OSFOC(NLAY,NCELL)'],
                 'detail':[''],'type':['arrfloat'],'default':[0.]},
        #'rct.4a':{'comm':'AW sorption a','cond':'OISOAW>0','kw':['OSAWA(NLAY,NCELL)'],
        #         'detail':[''],'type':['arrfloat'],'default':[0.]},
        #'rct.4b':{'comm':'AW sorption b','cond':'OISOAW>0','kw':['OSAWB(NLAY,NCELL)'],
        #         'detail':[''],'type':['arrfloat'],'default':[0.]},
        'rct.4':{'comm':'1st order rate','cond':'OIREAC>0','kw':['ORC1'],
                 'detail':[''],'type':['float'],'default':[0.]},
        'rct.5':{'comm':'1st rate sorbed','cond':'OIREAC>0','kw':['ORC2'],
                 'detail':[''],'type':['float'],'default':[0.]},
        #Solver
        'tslv':{'comm':'solver for C','cond':'','kw':['OSTSO','OSTPRE','OSTTOL','OSTRTOL','OSTDCMX'],
                 'detail':[['solver','PBiCG'],['preconditioner','DIC','DILU'],
                           'tolerance','relTol','dCmax'],
                 'type':['choice','choice','float','float','float'],
                 'default':[0,1,1e-12,0,0.02]},
        }
        
    def addKeyword(self,lname,kdef):
        for n in kdef:
            self.lines[lname][n].append(kdef[n])
        
#Chemistry
class  OpC:
    def __init__(self):
        self.grpList=['CHPRM','SOLU','GAS','CHSLV']
        self.groups={
        'CHPRM':['chprm'],
        'SOLU':['sactiv','sinit','sfix','swel','srch','sghb','sriv'],
        'GAS':['ginit','gfix','gwel','grch'],
        'CHSLV':['chslv'],
        }
        self.lines={
        'chprm':{'comm':'Chem parameters','cond':'','kw':['OCHNSTP','OCKOPT'],
                'detail':['reaction stp interval','kinet options'],
                'type':['int','Textlong'],'default':[10,'-tol 1e-8']},
        'sactiv':{'comm':'active zone','cond':'','kw':['OSACT'],
                'detail':[],'type':['arrfloat'],'default':[0.]},
        'sinit':{'comm':'Initial solutions','cond':'','kw':['OISOL'],
                'detail':[],'type':['arrfloat'],'default':[0.]},
        'sfix':{'comm':'Sources','cond':'','kw':['OSSOL'],
                'detail':[],'type':['arrfloat'],'default':[0.]},
        'swel':{'comm':'well solu nb','cond':'','kw':['OWSOL'],
                'detail':[],'type':['arrfloat'],'default':[0.]},
        'sghb':{'comm':'ghb solu nb','cond':'','kw':['OWSOL'],
                'detail':[],'type':['arrfloat'],'default':[0.]},
        'sriv':{'comm':'river solu nb','cond':'','kw':['OWSOL'],
                'detail':[],'type':['arrfloat'],'default':[0.]},
        'srch':{'comm':'Recharge solu nb','cond':'','kw':['ORSOL'],
                'detail':[],'type':['arrfloat'],'default':[0.]},
        'ginit':{'comm':'Initial gas','cond':'','kw':['OIGAS'],
                'detail':[],'type':['arrfloat'],'default':[0.]},
        'gfix':{'comm':'fixed gas','cond':'','kw':['OSGAS'],
                'detail':[],'type':['arrfloat'],'default':[0.]},
        'gwel':{'comm':'well gas','cond':'','kw':['OWGAS'],
                'detail':[],'type':['arrfloat'],'default':[0.]},
        'grch':{'comm':'Recharge gas','cond':'','kw':['ORGAS'],
                'detail':[],'type':['arrfloat'],'default':[0.]},
        'chslv':{'comm':'solver for Cwi','cond':'','kw':['OCHSO'],
                 'detail':[['solver','PBiCG'],['preconditioner','DIC','DILU'],
                           'tolerance','relTol'],
                 'type':['choice','choice','float','float'],'default':[0,1,1e-12,0]},
        }
