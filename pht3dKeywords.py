# -*- coding: utf-8 -*-
"""
Created on Tue Dec 24 09:45:48 2013

@author: olive
"""
class Ph:
    def __init__(self):
        self.grpList=['PH']
        self.groups = {
        'PH':['ph.'+str(i) for i in range(1,9)]
        }
        self.lines ={
        'ph.1':{'comm':'major flags','cond':'',
                'kw':['OS','TMP_LOC','RED_MODE','TEMP','ASBIN','EPSAQU','EPSPH','PACK_SZ'],
                'detail':[['splitting','-','iterative','calc flow steps','calc transp steps'],
                    ['write in','temp','local dir'],['fixing pe pH','all free','fix pe','fix pe and pH'],
                    'temperature',['output','only bin','bin and ascii'],'minimum conc','minimum pH','packet size'],
                'type':['choice','choice','choice','float','choice','float','float','int'],
                'default':[2,1,0,25,0,1e-10,1e-3,4000]},
        'ph.2':{'comm':'Charge imbalance','cond':'','kw':['CB_OFFSET'],
                'detail':'','type':['float'],'default':[0]},
        'ph.3':{'comm':'Initial chemistry','cond':'','kw':['PH_BTN(NLAY,NCOL,NROW)'],'detail':[],'type':['arrint']},
        'ph.4':{'comm':'Source / sink','cond':'','kw':['PH_SSM(NLAY,NCOL,NROW)'],'detail':[],'type':['arrint']},
        'ph.5':{'comm':'Recharge','cond':'','kw':['PH_RECH(NLAY,NCOL,NROW)'],'detail':[],'type':['arrint']},
        'ph.6':{'comm':'specific options','cond':'',
                'kw':['SU_OPT','NB_SOLU','NB_PHASE','NB_EXCH','NB_SURF'],
                'detail':['surface option','nb solutions','nb phases','nb exchange','nb surfaces'],
                'type':['string','int','int','int','int'],'default':['no_edl',4,4,4,4]},
        'ph.7':{'comm':'Evapotransp','cond':'','kw':['PH_EVT(NLAY,NCOL,NROW)'],'detail':[],'type':['arrint']},
        'ph.8':{'comm':'Gwet','cond':'','kw':['PH_GWET(NLAY,NCOL,NROW)'],'detail':[],'type':['arrint']},
        }