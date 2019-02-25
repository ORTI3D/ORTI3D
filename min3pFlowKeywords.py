# -*- coding: utf-8 -*-
"""
Created on Sun Oct 19 10:01:50 2014

@author: olive
""" 
class m3F:
    def __init__(self):
        self.grpList=['glo','spat','time','out','poro','flow','inif','bcf','conf']
        self.groups={
        'glo':['glo.1'],
        'spat':['spat.'+str(x) for x in range(1,9)],
        'time':['time.1'],
        'out':['out.1'], #,'out.2','out.3'],
        'poro':['poro.1'],
        'flow':['flow.1','flow.2'],
        'inif':['inif.1'],
        'bcf':['bcf.1','bcf.2','bcf.3','bcf.4','bcf.5'],
        'conf':['conf.1','conf.2','conf.2a','conf.3','conf.4','conf.5','conf.6','conf.7','conf.8']
        }
        s='!--------------------------------------------------------------\n'
        self.longNames={'glo': s+'!Data Block 1: Global Control Parameters \n'+s+s+ '\'global control parameters', #GC 21/02/2019 + OA
            'spat' : s+'!Data Block 3: Spatial Discretization \n'+s+s+'\'spatial discretization', #GC 21/02/2019
            'time': s+'!Data Block 4: Time Step Control - Global System \n'+s+s+'\'time step control - global system', #GC 21/02/2019
            'out' : s+'!Data Block 8: Output Control \n'+s+s+'\'output control', #GC 21/02/2019
            'poro': s+'!Data Block 9: Physical Parameters - Porous Medium \n'+s+s+'\'physical parameters - porous medium', #GC 21/02/2019
            'flow': s+'!Data Block 10: Physical Parameters - Variably Saturated Flow \n '+s+s+'\'physical parameters - variably saturated flow', #GC 21/02/2019
            'inif': s+'!Data Block 12: Initial Condition - Variably Saturated Flow \n '+s+s+'\'initial condition - variably saturated flow', #GC 21/02/2019
            'bcf' : s+'!Data Block 13: Boundary Condition - Variably Saturated Flow \n '+s+s+'\'boundary conditions - variably saturated flow', #GC 21/02/2019 
            'conf': s+'!Data Block 6: Control Parameters - Variably Saturated Flow \n -'+s+s+'\'control parameters - variably saturated flow', #GC 21/02/2019
        }
        self.lines={
        #'global control parameters' 
        'glo.1':{'comm':'Global','cond':'',
                 'kw':['isFlow','Fsteady','fully_saturated','isREaction','xiDiffusion','dens_dependt','energy_balance','ice_sheet','use_evapo'],
                'detail':[['Flow','true','false'],['Steady flow','true','false'],['Saturated Flow','true','false'],
                    ['Reactions','true','false'],
                    ['multicomponent diffusion','no','yes'],['density dependent flow','no','yes'],
                    ['energy balance','no','yes'],['compute ice sheet loading','no','yes'],
                    ['use evaporation block','no','yes']],
                'type':['choice','choice','choice','choice','choice','choice','choice','choice','choice']},
        #'spatial discretization' 
        'spat.1':{'comm':'x discretization','cond':'','kw':['NINTX','NX','Xmin','Xmax'],
                        'detail':['nb of intervals','nb of cells','X minimum','X maximum'],
                        'type':['int','int','float','float'],'default':[1,10,0.,1.],'units':['L']},
        'spat.2':{'comm':'y discretization','cond':'','kw':['NINTY','NY','Ymin','Ymax'],
                        'detail':['nb of intervals','nb of cells','Y minimum','Y maximum'],
                        'type':['int','int','float','float'],'default':[1,1,0.,1.],'units':['L']},
        'spat.3':{'comm':'z discretization','cond':'','kw':['NINTZ','NZ','Zmin','Zmax'],
                        'detail':['nb of intervals','nb of cells','Z minimum','Z maximum'],
                        'type':['int','int','float','float'],'default':[1,1,0.,1.],'units':['L']},
        'spat.4':{'comm':'unstructured','cond':'','kw':['P_Uns'],
                        'detail':[['','no','yes']],'type':['choice'],'default':[0]},
        'spat.5':{'comm':'mesh building','cond':'P_Uns!=0','kw':['P_Umake','P_Unumb','P_Umeth','P_Uobtus'],
                  'detail':[['','make from parms','make from bc','make from file',
                        'read unstructured grid from file'],
                        ['renumber from file','no','yes'],
                        ['control volume method','voronoi diagram','median dual'],
                        ['allow obtuse cells','no','yes']],
                    'print':[1,0,1,0],
                    'type':['choice','choice','choice','choice'],'default':[3,0,0,1]},
        'spat.6': {'comm':'mesh','cond':'P_Uns!=0','kw':['P_MESH'],'detail':[''],
                   'type':['arrfloat']},
        'spat.7': {'comm':'top','cond':'','kw':['P_TOP'],'detail':[''],
                   'type':['arrfloat']},
        'spat.8': {'comm':'bottom','cond':'','kw':['P_BOTM'],'detail':[''],
                   'type':['arrfloat']},
        #'time step control' 
        'time.1':{'comm':'time control','cond':'','kw':['Tunit','Tstart','Tfinal','Tmaxstep','Tminstep'],
                  'detail':[['time unit','years','days','hours','minutes'],
                        'time at start','final time','max time step','min time step'],
                  'type':['choice','float','float','float','float'],
                  'default':[1,0.,10.,.1,1e-7]},
        #'output control
        'out.1':{'comm':'output of spatial data','cond':'','kw':['Outs'],
                  'detail':[''],'type':['string'],'default':['']},
        #'physical parameters – porous medium'  
        'poro.1':{'comm':'porosity','cond':'','kw':['Poro'],'detail':[],'type':['arrfloat'],'default':[0.25]},
        #'physical parameters – variably saturated flow'  
        'flow.1':{'comm':'hydraulic conductivity','cond':'Uperm==0','kw':['Kxyz'],
                  'detail':[''],'type':['arrint'],
                  'default':[1e-4,1e-4,1e-4,0,.1,25.,3.,.5,0.],
                  'names':['K x (ms-1)','K y (ms-1)','K z (ms-1)','Ss','Swr','a_vg','n_vg','l_vg','Pe'],
                  'longNames':['hydraulic conductivity in x-direction',
                    'hydraulic conductivity in y-direction',
                    'hydraulic conductivity in z-direction',
                    'specific storage coefficient',
                    'soil hydraulic function parameters','','','',''
                    ]},
        'flow.2':{'comm':'permeability zones','cond':'Uperm==1','kw':['kxx'],
                  'detail':[],'type':['arrint'],'default':[1e-10,1e-10,1e-10],
                'names':['k x (m2)','K y (m2)','K z (m2)'],      
                'longNames':['permeability in x-direction','permeability in y-direction','permeability in z-direction']},
        #'initial condition – variably saturated flow' 
        'inif.1':{'comm':'initial condition','cond':'','kw':['Hinit'],
                  'detail':[],'type':['arrfloat'],'default':[10.]},
        #'boundary condition – variably saturated flow' 
        'bcf.1':{'comm':'first(BC head)','cond':'','kw':['BChead'],
                  'detail':[],'type':['arrfloat'],'default':[10.]},
        'bcf.2':{'comm':'second(BC flux) (m/s)','cond':'','kw':['BCflux'],
                  'detail':[],'type':['arrfloat'],'default':[0.]},
        'bcf.3':{'comm':'third(seepage)','cond':'','kw':['BCseep'],
                  'detail':[],'type':['arrfloat'],'default':[0.]},
        'bcf.4':{'comm':'point(flux)','cond':'','kw':['BChdpoint'],
                  'detail':[],'type':['arrfloat'],'default':[0.]},
        'bcf.5':{'comm':'atmospheric','cond':'','kw':['BCatmo'],
                  'detail':[],'type':['arrfloat'],'default':[0.]},
        #'control parameters – variably saturated flow'
        'conf.1':{'comm':'mass balance','cond':'','kw':['conf1'],
                  'detail':[['use','no','yes']],'type':['choice'],'default':[0]},
        'conf.2':{'comm':'input units for boundary and initial conditions','cond':'',
                  'kw':['conf2'],'detail':[['type','hydraulic head','pressure head']],
                'type':['choice'],'default':[0]},
        'conf.2a':{'comm':'input units for media permeability','cond':'',
                  'kw':['Uperm'],'detail':[['type','hydraulic conductivity','permeability']],
                'type':['choice'],'default':[0]},
        'conf.3':{'comm':'centered weighting','cond':'','kw':['conf3'],
                  'detail':[''],'type':['string'],'default':['']},
        'conf.4':{'comm':'compute underrelaxation factor','cond':'','kw':['conf4'],
                  'detail':['maximum update'],'type':['float'],'default':[10.]},
        'conf.5':{'comm':'newton iteration settings','cond':'',
                  'kw':['conf5a','conf5b','conf5c','conf5d'],
                  'detail':['num increment','max nb of iteration','converg tolerance',
                            'saturation change'],
                'type':['float','int','float','float'],
                'default':[1e-6,100,1e-6,0.1]},
        'conf.6':{'comm':'solver settings','cond':'',
                  'kw':['conf6a','conf6b','conf6c','conf6d','conf6e'],
                  'detail':['fact level','nb of iteration','info level',
                        'residual tolerance','update tolerance'],
                'type':['int','int','float','float','float'],
                'default':[0,100,1,1e-7,1e-7]},
        'conf.7':{'comm':'other parameters','cond':'','kw':['conf7a','conf7b'],
                  'detail':[['natural ordering','no','yes'],['iterative solver','no','yes']],
                  'type':['choice','choice'],'default':[0,0]},
        'conf.8':{'comm':'variable density parameters','cond':'dens_dependt!=0',
                  'kw':['conf8a','conf8b','conf8c','conf8d','conf8e','conf8f'],
                  'detail':['reference density','drho/dC','max Picard it.','target Picard it.',
                            'Picard tolerance','target curant nb'],
                'type':['float','float','int','int','float','float'],
                'default':[1e3,0.,4,2,0.1,0.1]},
        }
 
