# -*- coding: utf-8 -*-
"""
Created on Sun Oct 19 10:01:50 2014

@author: olive
""" 
class m3T:
    def __init__(self):
        self.grpList=['trans','trac','init','bct','cont','engp','inie','bce','cone','conv']
        self.groups={
        'trans':['trans.1','trans.2'],
        'trac':['trac.1','trac.2','trac.3'],
        'init':['init.1'],
        'bct':['bct.1','bct.2','bct.3','bct.4'],
        'cont':['cont.'+str(x) for x in range(1,12)],
        'engp':['engp.1','engp.2','engp.3'],
        'inie':['inie.1','inie.2'],
        'bce':['bce.1','bce.2','bce.3','bce.4','bce.5'],
        'cone':['cone.'+str(x) for x in range(1,10)],
        'conv':['conv.'+str(x) for x in range(14)],
        }
        self.longNames={'trans':'physical parameters - reactive transport',
                   'trac':'geochemical system',
                   'init':'initial condition - reactive transport',
                   'bct':'boundary conditions - reactive transport', 
                   'cont':'control parameters - reactive transport',
                   'engp':'physical parameters - energy balance',
                   'inie':'initial condition - energy balance',
                   'bce':'boundary conditions - energy balance', 
                   'cone':'control parameters - energy balance',
                   'conv':'control parameters - evaporation',
        }
        self.lines={
        #'physical parameters – reactive transport'  
        'trans.1':{'comm':'diffusion coefficients','cond':'','kw':['Diff_w','Diff_a','Diff_choice'],
                  'detail':['Diff in water (m2/s)','Diff in air (m2/s)',['Diffusion type','Original','Binary','Dusty gas']],
                   'type':['float','float','choice'],
                   'default':[1e-10,1e-5,0]},
        'trans.2':{'comm':'dispersivities','cond':'','kw':['Dsp'],
                  'detail':[''],'type':['arrint'],'default':[1.,0.1,0.1],
                  'names':['DspL','DspT','DspTV'],
                  'longNames':['longitudinal dispersivity',
                       'transverse horizontal dispersivity',
                       'transverse vertical dispersivity']},
        
        #tracer geochemistry
        'trac.1':{'comm':'use new database format','cond':'','kw':['Trac1'],
                  'detail':[''],'type':['string'],'default':['\n']},
        'trac.2':{'comm':'database directory','cond':'','kw':['Trac2'],
                  'detail':[''],'type':['string'],'default':['\'\'']},
        'trac.3':{'comm':'components','cond':'','kw':['Trac3'],
                  'detail':[''],'type':['string'],'default':['1 \n\'tracer\'']},
        #'initial condition – reactive transport' 
        'init.1':{'comm':'concentration input','cond':'','kw':['Cinit'],
                  'detail':[''],'type':['arrfloat'],'default':[1e-9],
                'suffx': '\'free\''},
        #'boundary condition – reactive transport' 
        'bct.1':{'comm':'first(BC fixed)','cond':'','kw':['BCfix'],
                  'detail':[''],'type':['arrfloat'],'default':[0.],
                'prefx':'\'concentration input\'\n','suffx': '\'free\''},
        'bct.2':{'comm':'second(free exit)','cond':'','kw':['BCfree'],
                  'detail':[''],'type':['arrfloat'],'default':[0.],
                'prefx':'\'concentration input\'\n','suffx': '\'free\''},
        'bct.3':{'comm':'third(mass flux)','cond':'','kw':['BCmass'],
                  'detail':[''],'type':['arrfloat'],'default':[0.],
                'prefx':'\'concentration input\'\n','suffx': '\'free\''},
        'bct.4':{'comm':'mixed(gas...)','cond':'','kw':['BCmix'],
                  'detail':[''],'type':['arrfloat'],'default':[0.],
                'prefx':'\'concentration input\'\n','suffx': '\'free\''},
        #control parameters reactive transport
        'cont.1':{'comm':'mass balance','cond':'','kw':['cont1'],
                  'detail':[''],'type':['string'],'default':['']},
        'cont.2':{'comm':'spatial weighting','cond':'','kw':['cont2'],
                  'detail':[['type','upstream','centered','van leer']],
                'type':['choice'],'default':[0]},
        'cont.3':{'comm':'activity update settings','cond':'','kw':['cont3'],
                  'detail':[['type','no update','time lagged','double update']],
                'type':['choice'],'default':[0]},
        'cont.4':{'comm':'tortuosity correction','cond':'','kw':['cont4'],
                  'detail':[['type','millington','no correction']],
                'type':['choice'],'default':[0]},
        'cont.5':{'comm':'degassing','cond':'','kw':['cont5'],
                  'detail':['degassing rate'],'type':['float'],'default':[0.]},
        'cont.6':{'comm':'update porosity','cond':'','kw':['cont6'],
                  'detail':[['use','no','yes']],'type':['choice'],'default':[0]},
        'cont.7':{'comm':'update permeability','cond':'','kw':['cont7'],
                  'detail':[['use','no','yes']],'type':['choice'],'default':[0]},
        'cont.8':{'comm':'user-specified underrelaxation factor','cond':'','kw':['cont8'],
                  'detail':['factor'],'type':['float'],'default':[1.]},
        'cont.9':{'comm':'newton iteration settings','cond':'',
                  'kw':['cont9a','cont9b','cont9c','cont9d','cont9e','cont9f'],
                  'detail':['num increment','nb of iteration','max nb of iteration',
                        'conc updates','max conc updates','conc tolerance'],
                'type':['float','int','int','float','float','float'],
                'default':[1e-6,12,15,0.5,1.,1e-6]},
        'cont.10':{'comm':'solver settings','cond':'',
                  'kw':['cont10a','cont10b','cont10c','cont10d','cont10e'],
                  'detail':['fact level','nb of iteration','info level',
                        'residual tolerance','update tolerance'],
                'type':['float','int','int','float','float','float'],
                'default':[0,100,1,1e-7,1e-7]},
        'cont.11':{'comm':'natural ordering','cond':'','kw':['cont11'],
                  'detail':[''],'type':['string'],'default':['']},
        #'physical parameters – energy balance'  
        'engp.1':{'comm':'specific heat of water','cond':'','kw':['Wheat'],'detail':[''],
                  'type':['float'],'default':[4182]},
        'engp.2':{'comm':'specific heat of air','cond':'','kw':['Aheat'],'detail':[''],
                  'type':['float'],'default':[1005]},
        'engp.3':{'comm':'Thermal properties','cond':'','kw':['Thermp'],'detail':[''],
                  'type':['arrint'],
                  'names':['spHeat sol','wat Conduct x','wat Conduct y','wat Conduct z',
                           'sol Conduct x','sol Conduct y','sol Conduct z',
                           'long dispers','transv disp','vert disp','b. dens'],
                 'default':[840,0.58,0.58,0.58,3.5,3.5,3.5,1,0.1,0.1,2000],
                 'longNames':['specific heat of solid',
                        'water thermal conductivity in x-direction','water thermal conductivity in y-direction','water thermal conductivity in z-direction',
                        'solid thermal conductivity in x-direction','solid thermal conductivity in y-direction','solid thermal conductivity in z-direction',
                        'longitudinal dispersivity','transverse vertical dispersivity','transverse horizontal dispersivity',
                        'solid bulk density']},
        #'initial condition – energy balance' 
        'inie.1':{'comm':'initial condition(T in C)','cond':'','kw':['Tinit'],
                  'detail':[''],'type':['arrfloat'],'default':[20]},
        'inie.2':{'comm':'geothermic gradient','cond':'','kw':['GGrad','GGtop'],
                  'detail':['gradient','temp'],'type':['float','float'],'default':[0.03,17]},
        #'boundary condition – energy balance' 
        'bce.1':{'comm':'first(T fixed)','cond':'','kw':['BCEfirst'],
                  'detail':[''],'type':['arrfloat'],'default':[0.]},
        'bce.2':{'comm':'second(sp. heat flux)','cond':'','kw':['BCEscnd'],
                  'detail':[''],'type':['arrfloat'],'default':[0.]},
        'bce.3':{'comm':'point(pt heat flux)','cond':'','kw':['BCEpoint'],
                  'detail':[''],'type':['arrfloat'],'default':[0.]},
        'bce.4':{'comm':'free(free exit)','cond':'','kw':['BCEfree'],
                  'detail':[''],'type':['arrfloat'],'default':[0.]},
        'bce.5':{'comm':'gradient','cond':'','kw':['BCEgrad'],
                  'detail':[''],'type':['arrfloat'],'default':[0.]},
        #control parameters energy balance
        'cone.1':{'comm':'Main parameters','cond':'','kw':['cone1a','cone1b','cone1c'],
                'detail':[['energy balance','no','yes'],['compute evaporation','no','yes'],
               ['non-linear density','no','yes']],
                'type':['choice','choice','choice'],
                'default':[0,0,0]},
        'cone.2':{'comm':'update viscosity','cond':'','kw':['cone2'],
                  'detail':[['viscocity model','no','sutra','diersch']],
                'type':['choice'],'default':[0]},
        'cone.3':{'comm':'spatial weighting','cond':'','kw':['cone3'],
                  'detail':[['choose','upstream','centered','van leer']],
                'type':['choice'],'default':[0]}, 
        'cone.4':{'comm':'reference tds','cond':'','kw':['cone4'],
                  'detail':[''],'type':['float'],'default':[0]},
        'cone.5':{'comm':'reference temperature for density','cond':'','kw':['cone5'],
                  'detail':[''],'type':['float'],'default':[20]},
        'cone.6':{'comm':'energy balance parameters','cond':'','kw':['cone6'],
                  'detail':[''],'type':['float'],'default':[-0.375]},
        'cone.7':{'comm':'thermal conductivity model','cond':'','kw':['cone7a','cone7b'],
                  'detail':[['type','model 1','model 2','model 3','model 4','model 5','model 6'],'value'],
                'type':['choice','float'],'default':[0,0.]},
        'cone.8':{'comm':'newton iteration settings','cond':'',
                  'kw':['cone10a','cone10b','cone10c','cone10d','cone10e','cone10f'],
                  'detail':['num increment','nb of iteration','max nb of iteration',
                        'conc updates','max conc updates','conc tolerance'],
                'type':['float','int','int','float','float','float'],
                'default':[1e-6,12,15,0.5,1.,1e-6]},
        'cone.9':{'comm':'solver settings','cond':'',
                  'kw':['cone11a','cone11b','cone11c','cone11d','cone11e'],
                  'detail':['fact level','nb of iteration','info level',
                        'residual tolerance','update tolerance'],
                'type':['float','int','int','float','float','float'],
                'default':[0,100,1,1e-7,1e-7]}, 
        # evapotranspiration control parameters
        'conv.0':{'comm':'write transient evaporation info','cond':'','kw':['conv0'],
                  'detail':[['use','no','yes']],'type':['choice'],'default':[0]},
        'conv.1':{'comm':'vapour density model','cond':'','kw':['conv1a','conv1b'],
                  'detail':[['model','saito et al. (2006)','default'],['capillarity correction','true','false']],
                 'type':['choice','choice'],'default':[0,0]},
        'conv.2':{'comm':'temperature gain factor for soil','cond':'','kw':['conv2'],
                  'detail':['value'],'type':['float'],'default':[7.]},
        'conv.3':{'comm':'update vapor density derivatives','cond':'','kw':['conv3'],
                  'detail':[['use','no','yes']],'type':['choice'],'default':[0]},
        'conv.4':{'comm':'reference vapor diffusivity','cond':'','kw':['conv4'],
                  'detail':['value'],'type':['float'],'default':[2.12e-5]},
        'conv.5':{'comm':'compute enhanced factor in thermal vapor fluxes','cond':'',
                  'kw':['conv5a','conv5b','conv5c'],
                  'detail':[['use','no','yes'],'clay fraction','empirical factor'],
                  'type':['choice','float','float'],'default':[0,0.5,8.]},
        'conv.6':{'comm':'soil surface resistance to vapor flow','cond':'','kw':['conv6'],
                  'detail':[['model','model 1','model 2','model 3','model 4']],
                  'type':['choice'],'default':[0]},
        'conv.7':{'comm':'split divergence of vapor density','cond':'','kw':['conv7'],
                  'detail':[['use','no','yes']],'type':['choice'],'default':[0]},
        'conv.8':{'comm':'tortuosity model to vapor flow','cond':'','kw':['conv8'],
                  'detail':[['type','millington','not consider']],'type':['choice'],'default':[0]},
        
        'conv.9':{'comm':'relative humidity parameters','cond':'',
                  'kw':['conv9a','conv9b','conv9c','conv9d','conv9e','conv9f'],
                  'detail':[['read from file','false','true'],'average Hr value (yr)',
                        'year amplitude','day amplitude','time max along year','time max along day'],
                'type':['choice','float','float','float','float','float','float'],
                'default':[0,0.85,0,0,182.5,0.5]},
        'conv.10':{'comm':'temperature parameters','cond':'',
                  'kw':['conv10a','conv10b','conv10c','conv10d','conv10e','conv10f'],
                  'detail':[['read from file','false','true'],'average Temp value (yr)',
                        'year amplitude','day amplitude','time max along year','time max along day'],
                'type':['choice','float','float','float','float','float','float'],
                'default':[0,20,0,0,182.5,0.5]},
        'conv.11':{'comm':'solar radiation parameters','cond':'',
                  'kw':['conv11a','conv11b','conv11c','conv11d','conv11e','conv11f','conv11g','conv11h','conv11i','conv11j'],
                  'detail':[['read radiation from file','false','true'],['read cloud I from file','false','true'],
                            'fact. to multiply radiation','latitude','time of expt start',
                        'time corresp to noon','time autumn starts','cloud index',
                        'albedo dry soil','albedo wet soil'],
                'type':['choice','choice','float','float','float','float','float','float','float','float'],
                'default':[0,0,0.,45.,0,0.5,273.5,0.5,0.2,0.1]},
        'conv.12':{'comm':'rain parameters','cond':'','kw':['conv12a','conv12b'],
                  'detail':[['read from file','false','true'],'leakage coefficient'],
                  'type':['choice','float'],'default':[0,-100.0]},
        'conv.13':{'comm':'evaporation parameters','cond':'',
                  'kw':['conv13a','conv13b','conv13c','conv13d','conv13e','conv13f','conv13g','conv13h','conv13i','conv13j',
                        'conv13k','conv13l','conv13m','conv13n'],
                  'detail':[['compute from aerodynamic','false','true'],['read from file','false','true'],
                            'imposed value','multiply factor','roughness length','screen height','stability factor',
                        'density of air',['Read wind from file','false','true'],
                        'wind speed average','wind year amplitude','day amplitude','yearly max','daily max'],
                'type':['choice','choice','float','float','float','float','float','float','choice','float','float','float','float','float'],
                'default':[0,0,0.,1.,1e-3,1e-3,1.,1.112,0,0.,0,0,0,0]},
        }