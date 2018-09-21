# -*- coding: utf-8 -*-
"""
Created on Wed Aug 27 10:54:19 2014

@author: olive
what is needed 
transport : porosity, disp (aL array, other ratio like in mt3d)
"""
grpList=['TRANS','RCT']
groups={
'TRANS':['trans.'+str(a) for a in range(1,10)],
'RCT':['rct.'+str(a) for a in range(1,3)],
}
lines={
'trans.1':{'comm':'Transport parameters','cond':'','kw':['F_P1'],
           'detail':['parm1','val1','val2'],
            'type':['choice'],'default':[0]},
'trans.2':{'comm':'Longit dispersivity','cond':'','kw':['F_AL'],'detail':['aL'],
            'type':['float'],'default':[1.],'units':['L']},
'trans.3':{'comm':'Transverse dispersivity','cond':'','kw':['F_AT'],'detail':['aT'],
            'type':['float'],'default':[0.1],'units':['L']},
'trans.4':{'comm':'Vertical dispersivity','cond':'','kw':['F_AZ'],'detail':['aZ'],
            'type':['float'],'default':[0.1],'units':['L']},
'trans.5':{'comm':'Diffusion coeff','cond':'','kw':['F_DIFF'],'detail':['Diffusion'],
            'type':['float'],'default':[0.],'units':['L^2/T']},
'trans.6':{'comm':'Initial concentration','cond':'','kw':['F_CINIT'],'detail':[''],
            'type':['arrfloat'],'default':[0.]},
'trans.7':{'comm':'Fixed concentration (1st)','cond':'','kw':['F_CFIXED'],'detail':[''],
            'type':['arrfloat'],'default':[0.]},
'trans.8':{'comm':'Fixed mass flux (2nd)','cond':'','kw':['F_CMASS'],'detail':[''],
            'type':['arrfloat'],'default':[0.]},
'trans.9':{'comm':'Fixed gradient (3rd)','cond':'','kw':['F_CGRAD'],'detail':[''],
            'type':['arrfloat'],'default':[0.]},
'rct.1': {'comm':'sorption coeff','cond':'','kw':['F_KD'],
        'detail':['Kd'],'type':['float']},
'rct.2': {'comm':'degradation rate','cond':'','kw':['F_LAMBDA'],
        'detail':['deg rate'],'type':['float']},
}
