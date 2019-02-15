# -*- coding: utf-8 -*-
"""
Created on Wed Aug 27 10:54:19 2014

@author: olive
"""
class ogF:
    def __init__(self):
        self.grpList=['DOMN','FLOW','UNSAT'] #,'DENS','XIPHASE']
        self.groups={
        'DOMN':['domn.'+str(a) for a in range(1,9)],
        'FLOW':['flow.'+str(a) for a in range(1,9)],
        'UNSAT':['unsat.'+str(a) for a in range(1,8)],
        }
        self.lines={
        'domn.1': {'comm':'dimensions','cond':'','kw':['O_GRID','NLAY','NCOL','NROW','EPS'],
                'detail':[['type of grid','rectangular','triangle'],
                          'nb of layers','nb of columns','nb of rows','epsilon distance'],
                'type':['choice','int','int','int','float']},
        'domn.2':{'comm':'Col width','cond':'','kw':['O_DELR(NCOL)'],'detail':['Col width'],
                        'type':['vecfloat'],'default':[10],'units':['L']},
        'domn.3':{'comm':'Row height','cond':'','kw':['O_DELC(NROW)'],'detail':['Row height'],
                        'type':['vecfloat'],'default':[10],'units':['L']},
        'domn.4': {'comm':'mesh','cond':'','kw':['O_MESH'],
                'detail':[''],'type':['arrfloat']},
        'domn.5': {'comm':'top','cond':'','kw':['O_TOP'],
                'detail':[''],'type':['arrfloat'],'default':[10.],'units':['L']},
        'domn.6': {'comm':'bottom','cond':'','kw':['O_BOTM'],
                'detail':[''],'type':['arrfloat'],'default':[0.],'units':['L']},
        'domn.7': {'comm':'units','cond':'','kw':['TUNIT','LUNIT','MUNIT'],
                'detail':[['time units','-','years','days','hours','seconds'],['Length unit','-','M','ft'],['Mass unit','-','kg','lb']],
                'type':['choice','choice','choice'],'default':[2,1,1]},
        'domn.8': {'comm':'saturation type','cond':'','kw':['O_UNSAT'],
                'detail':[['saturation','-','saturated','unsaturated']],
                'type':['choice',],'default':[1]},
        #flow : 
        'flow.1': {'comm':'initial head','cond':'','kw':['O_HINIT'],
                'detail':[''],'type':['arrfloat']},
        'flow.2': {'comm':'fixed head (1st)','cond':'','kw':['O_HFIXED'],
                'detail':[''],'type':['arrfloat']},
        'flow.3': {'comm':'water flux (2nd)','cond':'','kw':['O_FLUX'],
                'detail':[''],'type':['arrfloat']},
        'flow.4': {'comm':'head gradient (3rd)','cond':'','kw':['O_HGRAD'],
                'detail':[''],'type':['arrfloat']},
        'flow.5': {'comm':'Medium number','cond':'','kw':['O_MEDN'],
                'detail':[''],'type':['arrfloat'],
                'names':['K. Hydraul','Storage','Porosity'],
                'default':[1e-4,1e-6,0.25]},
        'flow.6': {'comm':'Ratio Horiz Vertic K','cond':'','kw':['O_KHKV'],
                'detail':[''],'type':['float'],'default':[1.]},
        'flow.7': {'comm':'Flux (wells)','cond':'','kw':['O_WELL'],
                'detail':[''],'type':['arrfloat']},
        'flow.8': {'comm':'Recharge','cond':'','kw':['O_RECH'],
                'detail':[''],'type':['arrfloat']},
        # unsaturated medium (not mixed with flow because of different units)
        'unsat.1': {'comm':'initial pressure','cond':'','kw':['O_PINIT'],
                'detail':[''],'type':['arrfloat']},
        'unsat.2': {'comm':'fixed pressure (1st)','cond':'','kw':['O_PFIXED'],
                'detail':[''],'type':['arrfloat']},
        'unsat.3': {'comm':'water flux (2nd)','cond':'','kw':['O_FLUX'],
                'detail':[''],'type':['arrfloat']},
        'unsat.4': {'comm':'pressure gradient (3rd)','cond':'','kw':['O_PGRAD'],
                'detail':[''],'type':['arrfloat']},
        'unsat.5': {'comm':'Medium properties','cond':'','kw':['O_MEDN'],
                'detail':[''],'type':['arrfloat'],
                'names':['Permeability','Porosity','Smin','Smax',  'alpha_VG','m_VG','k_exponent'],
                'default':[1e-10,0.25,0.05,1.0,5.,0.37,0.5]},
        'unsat.6': {'comm':'Ratio Horiz Vertic K','cond':'','kw':['O_KHKV1'],
                'detail':[''],'type':['float'],'default':[1.]},
        'unsat.7': {'comm':'Wells','cond':'','kw':['O_WELL1'],
                'detail':[''],'type':['arrfloat']},
        }
