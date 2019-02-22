# -*- coding: utf-8 -*-
"""
Created on Sun Oct 19 10:01:50 2014

@author: olive
""" 
class m3C:
    def __init__(self):
        self.grpList=['geoch','inic','bcc','conc']
        self.groups={
        'geoch':['geoch.1','geoch.2','geoch.3'],
        'inic':['inic.1'],
        'bcc':['bcc.1','bcc.2','bcc.3','bcc.4'],
        'conc':['conc.1','conc.2','conc.3','conc.4']
        }
        s='--------------------------------------------------------------\n'
        self.longNames={'geoch': '!-'+s+'!Data Block 2: Geochemical System '+s+s+'geochemical system', #GC 21/02/2019 (+OA)
            'inic' : '!-'+s+'!Data Block 16: Initial Condition - Reactive Transport '+s+s+' initial condition - reactive transport', #GC 21/02/2019
            'bcc'  : '!-'+s+'!Data Block 15: Boundary Condition - Reactive Transport '+s+s+' boundary conditions - reactive transport', #GC 21/02/2019 
            'conc' : '!-'+s+'!Data Block 5: Control Parameters - Local Geochemistry '+s+s+' control parameters - local geochemistry', #GC 21/02/2019
                        }
        self.lines={
        #geochemistry # components are obtained from chem database
        'geoch.1':{'comm':'use new database format','cond':'','kw':['Geoch1'],
                  'detail':[],'type':['string'],'default':['\n']},
        'geoch.2':{'comm':'database directory','cond':'','kw':['Geoch2'],
                  'detail':[],'type':['string'],'default':['\'\'']},
        'geoch.3':{'comm':'components','cond':'','kw':['Geoch3'],
                  'detail':[],'type':['string'],'default':['1 \n\'tracer\'']},
        #'initial condition – reactive transport' #here we store solution numbers
        'inic.1':{'comm':'concentration input','cond':'','kw':['Cinit'],
                  'detail':[],'type':['arrint'],'default':[0]},
        #'boundary condition – reactive transport' #here we store solution numbers
        'bcc.1':{'comm':'first(BC fixed)','cond':'','kw':['BCfix'],
                  'detail':[],'type':['arrint'],'default':[0]},
        'bcc.2':{'comm':'second(free exit)','cond':'','kw':['BCfree'],
                  'detail':[],'type':['arrint'],'default':[0]},
        'bcc.3':{'comm':'third(mass flux)','cond':'','kw':['BCmass'],
                  'detail':[],'type':['arrint'],'default':[0]},
        'bcc.4':{'comm':'mixed(gas...)','cond':'','kw':['BCmix'],
                  'detail':[],'type':['arrint'],'default':[0]},
        #control parameters chemistry
        'conc.1':{'comm':'newton iteration settings','cond':'','kw':['conc1a','conc1b'],
                  'detail':[],'type':['float','float'],'default':[1e-6,1e-6]},
        'conc.2':{'comm':'output time untis','cond':'','kw':['conc2'],
                  'detail':[['years','days','minutes','seconds']],
                'type':['choice'],'default':[1]},
        'conc.3':{'comm':'maximum ionic strength','cond':'','kw':['conc3'],
                  'detail':[],'type':['float'],'default':[2]},
        'conc.4':{'comm':'minimum activity for h2o','cond':'','kw':['conc4'],
                  'detail':[],'type':['float'],'default':[0.5]},
        }
