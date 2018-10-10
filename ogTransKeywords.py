# -*- coding: utf-8 -*-
"""
Created on Wed Aug 27 10:54:19 2014

@author: olive
"""
class ogT:
    def __init__(self):
        self.grpList=['Trans'] #,'DENS','UNSAT','XIPHASE']
        self.groups={
        'Trans':['trans.'+str(a) for a in range(1,10)],
        }
        self.lines={
        'trans.1': {'comm':'initial concentration','cond':'','kw':['O_CINIT'],
                'detail':[''],'type':['arrfloat']},
        'trans.2': {'comm':'fixed conc.','cond':'','kw':['O_CFIXED'],
                'detail':[''],'type':['arrfloat']},
        'trans.3': {'comm':'fixed flux','cond':'','kw':['O_CFLUX'],
                'detail':[''],'type':['arrfloat']},
        'trans.4': {'comm':'conc. gradient','cond':'','kw':['O_CGRAD'],
                'detail':[''],'type':['arrfloat']},
        'trans.5': {'comm':'Longitud. dispersivity','cond':'','kw':['O_AL'],
                'detail':[''],'type':['float'],'default':[1.]},
        'trans.6': {'comm':'lat. dispers.','cond':'','kw':['O_AT'],
                'detail':[''],'type':['float'],'default':[0.1]},
        'trans.7': {'comm':'vert. dispers.','cond':'','kw':['O_AZ'],
                'detail':[''],'type':['float'],'default':[0.1]},
        'trans.8': {'comm':'diffusion coeff','cond':'','kw':['O_DIFFU'],
                'detail':[''],'type':['float'],'default':[1e-9]},
        'trans.9': {'comm':'Recharge','cond':'','kw':['O_TRECH'],
                'detail':[''],'type':['arrfloat']},
        }
