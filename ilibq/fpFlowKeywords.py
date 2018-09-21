# -*- coding: utf-8 -*-
"""
Created on Wed Aug 27 10:54:19 2014

@author: olive
what is needed 
flow : K and Kh, BC(head), wells (can be distributed), averaging?
"""
grpList=['DOMN','FLOW'] #,'DENS','UNSAT','XIPHASE']
groups={
'DOMN':['domn.'+str(a) for a in range(1,8)],
'FLOW':['flow.'+str(a) for a in range(1,9)],
}
lines={
'domn.1': {'comm':'dimensions','cond':'','kw':['F_GRID','NLAY','NCOL','NROW'],
        'detail':[['type of grid','rectangular','triangle'],
                  'nb of layers','nb of columns','nb of rows'],
        'type':['choice','int','int','int']},
'domn.2':{'comm':'Col width','cond':'','kw':['DELR(NCOL)'],'detail':['Col width'],
                'type':['vecfloat'],'default':[10],'units':['L']},
'domn.3':{'comm':'Row height','cond':'','kw':['DELC(NROW)'],'detail':['Row height'],
                'type':['vecfloat'],'default':[10],'units':['L']},
'domn.4': {'comm':'mesh','cond':'','kw':['F_MESH'],
        'detail':[''],
        'type':['arrfloat']},
'domn.5': {'comm':'top','cond':'','kw':['F_TOP'],
        'detail':[''],'type':['arrfloat'],'default':[10.],'units':['L']},
'domn.6': {'comm':'bottom','cond':'','kw':['F_BOTM'],
        'detail':[''],'type':['arrfloat'],'default':[0.],'units':['L']},
'domn.7': {'comm':'units','cond':'','kw':['TUNIT','LUNIT','MUNIT'],
        'detail':[['time units','-','days','hour'],['Length unit','-','M','ft'],['Mass unit','-','kg','lb']],
        'type':['choice','choice','choice']},
'flow.1': {'comm':'initial head','cond':'','kw':['F_HINIT'],
        'detail':[''],'type':['arrfloat']},
'flow.2': {'comm':'fixed head','cond':'','kw':['F_HFIXED'],
        'detail':[''],'type':['arrfloat']},
'flow.3': {'comm':'Horiz hydraul conductivity','cond':'','kw':['F_KHOR'],
        'detail':[''],'type':['arrfloat'],'default':[1.]},
'flow.4': {'comm':'Ratio Horiz Vertic K','cond':'','kw':['F_KHKV'],
        'detail':[''],'type':['arrfloat'],'default':[1.]},
'flow.5': {'comm':'specific storage','cond':'','kw':['F_SS'],
        'detail':[''],'type':['arrfloat'],'default':[1e-6]},
'flow.6': {'comm':'porosity','cond':'','kw':['F_PORO'],
        'detail':[''],'type':['arrfloat'],'default':[.25]},
'flow.7': {'comm':'Flux (wells)','cond':'','kw':['F_WELL'],
        'detail':[''],'type':['arrfloat']},
'flow.8': {'comm':'Recharge','cond':'','kw':['F_RECH'],
        'detail':[''],'type':['arrfloat']},
}
