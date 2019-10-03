class  Mt:
    def __init__(self):
        self.grpList=['BTN','ADV','DSP','GCG','SSMs','RCT','VDF','VSC','UZT']
        self.groups={
        'BTN':['btn.'+str(a) for a in range(1,24)],
        'ADV':['adv.1','adv.2','adv.3','adv.4','adv.5'],
        'DSP':['dsp.1','dsp.2','dsp.3','dsp.4','dsp.5'],
        'GCG':['gcg.1','gcg.2'],
        'SSMs':['ssms.1'],
        'RCT':['rct.1','rct.2a','rct.2b','rct.2c','rct.3','rct.4','rct.5','rct.6'],
        'VDF':['vdf.1','vdf.2','vdf.3','vdf.4','vdf.4a','vdf.4b','vdf.4c','vdf.5','vdf.6','vdf.7'],
        'VSC':['vsc.1','vsc.2','vsc.3a','vsc.3b','vsc.3c','vsc.3d'],
        'UZT':['uzt.1','uzt.2','uzt.3','uzt.4','uzt.5','uzt.6','uzt.7','uzt.8','uzt.9','uzt.10'],
        }
        self.longNames={'BTN':'physical parameters - reactive transport',
                   'ADV':'advection','DSP':'Dispersion','GCG':'Solver',
                   'SSMs':'Specific source and sinks','RCT':'Reaction','VDF':'Variable density',
                   'VSC':'Variable viscosisty','UZT':'Unsaturated zone'
        }
        self.lines={
            #BTN
        'btn.1': {'comm':'Title','cond':'','kw':['tiBtn1'],'detail':[],'type':['title'],'default':['# title 1']},
        'btn.2': {'comm':'Title2','cond':'','kw':['tiBtn2'],'detail':[],'type':['title'],'default':['# title 2']},
        'btn.3': {'comm':'dimensions','cond':'',
                'kw':['NLAY','NROW','NCOL','NPER','NCOMP','MCOMP','IGSFLG','GCOMPN','KCOMPN'],
                'detail':['Nb of layers','Nb of rows','Nb of columns','Nb of periods','Nb of components',
                          'Nb mobile comp.',['gas used','no','yes'],'nb of gas compounds','nb of immob kinetic sp.'],
                'type':['int','int','int','int','int','int','choice','int','int']},
        'btn.4': {'comm':'units','cond':'','kw':['TUNIT','LUNIT','MUNIT'],
                'detail':[['time units','-','days','hour'],['Length unit','-','M','ft'],['Mass unit','-','kg','lb']],
                'type':['choice','choice','choice']},
        'btn.5': {'comm':'flags','cond':'','kw':['TRNOP'],'detail':['flags (adv, disp, ssm , chem gcg)'],
                'type':['string']},
        'btn.6':{'comm':'type of layer','cond':'','kw':['TLAYCON(NLAY)'],'detail':[],
                'type':['layint']},
        'btn.7':{'comm':'Col width','cond':'','kw':['DELR(NCOL)'],'detail':[],'type':['vecfloat']},
        'btn.8':{'comm':'Row height','cond':'','kw':['DELC(NROW)'],'detail':[],'type':['vecfloat']},
        'btn.9':{'comm':'Top','cond':'','kw':['HTOP(NCOL,NROW)'],'detail':[],'type':['arrfloat']},
        'btn.10':{'comm':'Thickness','cond':'','kw':['DZ(NLAY,NCOL,NROW)'],'detail':[],'type':['arrfloat']},
        'btn.11':{'comm':'Porosity','cond':'','kw':['PRSTY(NLAY,NCOL,NROW)'],'detail':[],
                'type':['arrfloat'],'default':[0.25]},
        'btn.12':{'comm':'Boundary conditions','cond':'','kw':['ICBUND(NLAY,NCOL,NROW)'],'detail':[],
                'type':['arrint'],'default':[1]},
        'btn.13':{'comm':'Concentrations','cond':'','kw':['SCONC(NLAY,NCOL,NROW)'],'detail':[],'type':['arrfloat']},
        'btn.14':{'comm':'Inactive cells and min saturation','cond':'','kw':['CINACT','THKMIN'],
                'detail':['Value for conc inactive cells','Minimum saturation in a cell'],
                'type':['float','float'],'default':[0,.01]},
        'btn.15':{'comm':'Printing flags','cond':'','kw':['IFMTCN','IFMTNP','IFMTRF','IFMTDP','SAVUCN'],
                'detail':[['print conc.','no','yes'],['print nb particles','no','yes'],
                          ['print retard coeff','no','yes'],['print dispersion','no','yes'],'Save to UCN'],
                'type':['choice','choice','choice','choice','string'],'default':[0,0,0,0,'T']},
        'btn.16':{'comm':'flag printing times','cond':'','kw':['NPRS'],
                  'detail':[['print','no','yes']],'type':['choice'],'default':[0]},
        'btn.17':{'comm':'Print times','cond':'NPRS>0','kw':['TIMPRS(NPRS)'],'detail':[''],'type':['string']},
        'btn.18':{'comm':'Observation points','cond':'','kw':['NOBS','NPROBS'],
                'detail':['Nb of observation points','Freq of printing'],'type':['int','int']},
        'btn.19':{'comm':'Obs points cells','cond':'NOBS>0','kw':['KOBS','IOBS','JOBS'],
                  'detail':[''],'type':['int','int','int']},
        'btn.20':{'comm':'mass printing','cond':'','kw':['CHKMAS','NPRMAS'],'detail':['',''],
                  'type':['string','int'],'default':['T',1]}, # OA 3/10/19 changed to true and 1 to remov. for file
        'btn.21':{'comm':'Periods','cond':'','kw':['PERLEN(NPER)','NSTP(NPER)','TSMULT(NPER)','SSFLAG'],
                'detail':['Period length','nb of time steps','multiplying factor','Flag steady'],
                'type':['float','int','float','string'],'default':[1.,1,1.05,'']},
        'btn.22':{'comm':'Periods','cond':'',
                  'kw':['DT0(NPER)','MXSTRN(NPER)','TTSMULT(NPER)','TTSMAX(NPER)'],
                 'detail':['start time','max T step in F steps','2nd multipl','2nd max tstep'],
                 'type':['float','int','float','float'],'default':[0.,10000,1.0,1.0]},
        'btn.23':{'comm':'Recharge','cond':'','kw':['MTRECH(NLAY,NCOL,NROW)'],'detail':[],'type':['arrfloat']},
        #'btn.24':{'comm':'Mass flux','cond':'','kw':['MTMFlux(NLAY,NCOL,NROW)'],'detail':[],'type':['arrfloat']},
        #ADVECTION
        'adv.1': {'comm':'Major data','cond':'','kw':['MIXELM','PERCEL','MXPART','NADVFD','IALTFM'],
                'detail':[['Method','TVD','Finite difference','MOC','MMOC','HMOC'],'Courant nb',
                          'Maximum nb particles','weighing scheme','formulation transient term'],
                'type':['choice','float','int','int','int'],'default':[0,.75,50000,1,0]},
        'adv.2': {'comm':'Particles','cond':'MIXELM>1','kw':['ITRACK','WD'],
                  'detail':['Type of tracking','weighing factor'],'type':['int','float'],'default':[1,.75]}, #EV 07/02/19
        'adv.3': {'comm':'Particles details','cond':'MIXELM in [2,4]',
                  'kw':['DCEPS','NPLANE','NPL','NPH','NPMIN','NPMAX'],
                'detail':['minimum gradient','initial placement','initial number 1','initial number 2','min nb allowed','max nb allowed'],
                'type':['float','int','int','int','int','int'],'default':[1e-5,1,0,16,2,32]},
        'adv.4': {'comm':'Particles details','cond':'MIXELM in [3,4]',
                  'kw':['INTERP','NLSINK','NPSINK'],
                'detail':['Interpol scheme','Flag inital sink part','Nb inital sink part'],
                'type':['int','int','int'],'default':[1,1,16]},
        'adv.5':{'comm':'Gradient to choose MOC scheme','cond':'MIXELM==4',
                 'kw':['DCHMOC'],'detail':[],'type':['float'],'default':[1.]},
        #dsp
        'dsp.1':{'comm':'multidiffusion','cond':'','kw':['IMDIFF'],'detail':[],'type':['textlong'],
                 'default':['#']},
        'dsp.2':{'comm':'Longitudinal disp','cond':'','kw':['AL(NLAY,NCOL,NROW)'],'detail':[],
                    'type':['arrfloat'],'units':['L'],'default':[1]},
        'dsp.3':{'comm':'Ratio Transv/Long disp','cond':'','kw':['TRPT(NCOL,NROW)'],'detail':[],
                    'type':['arrfloat'],'default':[.1]},
        'dsp.4':{'comm':'Ratio Vert/Long disp','cond':'','kw':['TRPV(NCOL,NROW)'],'detail':[],
                    'type':['arrfloat'],'default':[.05]},
        'dsp.5':{'comm':'Diffusion coeff','cond':'','kw':['DMCOEF(NCOL,NROW)'],'detail':[],
               'type':['arrfloat'],'units':['L2/T']},
        #GCG Solver
        'gcg.1':{'comm':'Main flags','cond':'','kw':['TMXITER','TITER1','ISOLVE','NCRS'],
                 'detail':['nb outer iterations','inner iterations',['type','-','Jacobi','SSOR','MIC'],['disp. tensor','lump','full']],
                'type':['int','int','choice','choice'],
                'default':[1,40,2,0]},
        'gcg.2':{'comm':'convergence','cond':'','kw':['ACCL','CCLOSE','IPRGCG'],
                 'detail':['acceleration','convergence','printing'],
                'type':['float','float','int'],
                'default':[1.,1e-11,0]},
        # SSM special
        'ssms.1':{'comm':'Specific Source/sink','cond':'','kw':['SSMSPEC(NCOL,NROW)'],'detail':[],
                    'type':['arrint'],'default':[-1]},
        #RCT
        'rct.1':{'comm':'major flags','cond':'','kw':['ISOTHM','IREACT','IRCTOP','IGETSC'],
                 'detail':[['sorption','no sorption','linear','freundlich','langmuir','kinetic sorption',
                            'dual domain no sorp','dual domain with sorp'],
                            ['kinetic reaction','no reaction','1st order','0th order'],'array format',
                        ['initial conc.','at equilibrium','given']],
                'type':['choice','choice','int','choice'],
                'default':[0,0,2,0]},
        'rct.2a':{'comm':'density','cond':'ISOTHM in [1,2,3,4,6]','kw':['RHOB(NLAY,NCOL,NROW)'],
                 'detail':['bulk density'],'type':['arrfloat'],'default':[1.8]},
        'rct.2b':{'comm':'immobile domain','cond':'ISOTHM>4','kw':['PRSITY2(NLAY,NCOL,NROW)'],
                 'detail':['porosity'],'type':['arrfloat'],'default':[1.]},
        'rct.2c':{'comm':'concentrations','cond':'IGETSC>0','kw':['SRCONC(NLAY,NCOL,NROW)'],
                 'detail':[''],'type':['arrfloat'],'default':[0.]},
        'rct.3':{'comm':'1st sorption parm','cond':'ISOTHM>0','kw':['SP1(NLAY,NCOL,NROW)'],
                 'detail':[''],'type':['arrfloat'],'default':[1.]},
        'rct.4':{'comm':'2nd sorption parm','cond':'ISOTHM>0','kw':['SP2(NLAY,NCOL,NROW)'],
                 'detail':[''],'type':['arrfloat'],'default':[1.]},
        'rct.5':{'comm':'1st order rate','cond':'IREACT>0','kw':['RC1(NLAY,NCOL,NROW)'],
                 'detail':[''],'type':['arrfloat'],'default':[1.]},
        'rct.6':{'comm':'1st rate sorbed','cond':'IREACT>0','kw':['RC2(NLAY,NCOL,NROW)'],
                 'detail':[''],'type':['arrfloat'],'default':[1.]},
        #VDF for seawat version 4
        ##        1. MTDNCONC MFNADVFD NSWTCPL IWTABLE*
        'vdf.1' : {'comm':'major flags','cond':'','kw':['MT3DRHOFLG','MFNADVFD','NSWTCPL','IWTABLE'],
                 'detail':['specified (0) or species nb','internode density calc.',
                           'coupling : iterative (0/1) or nsteps',
                        ['water table correct.','not applied','applied']],
                'type':['int','int','int','choice'],
                'default':[1,2,0,0]},
        'vdf.2' : {'comm':'density limits','cond':'','kw':['DENSEMIN','DENSEMAX'],
                 'detail':['minimum','maximum'],'type':['float','float'],'default':[0.,0.]},
        'vdf.3' : {'comm':'density tolerance','cond':'NSWTCPL>1','kw':['DNSCRIT'],
                 'detail':[''],'type':['float'],'default':[.1]},
        'vdf.4' : {'comm':'density curves','cond':'','kw':['DENSEREF','DRHODC'],
                 'detail':['reference density','slope of dens line'],
                'type':['float','float'],'default':[1000.,0.714]},
        'vdf.4a': {'comm':'density curves/species','cond':'MT3DRHOFLG<0',
                   'kw':['DENSEREF','DRHODPRHD','PRHDREF'],
                 'detail':['reference density','slope of dens line','density vs head'],
                'type':['float','float','float'],'default':[1000.,0.714,4.4e-3]}, 
        'vdf.4b' : {'comm':'nb of species','cond':'MT3DRHOFLG<0','kw':['NSRHOEOS'],
                 'detail':[''],'type':['int'],'default':[1]},
        'vdf.4c' : {'comm':'per species','cond':'MT3DRHOFLG<0',
                    'kw':['MTRHOSPEC(NSRHOEOS)','DRHODC(NSRHOEOS)','CRHOREF(NSRHOEOS)'],
                 'detail':['species nb','slope of dens line','ref conc'],
                'type':['int','float','float'],'default':[0,0.714,0.]},
        'vdf.5' : {'comm':'1st time step','cond':'','kw':['FIRSTDT'],
                 'detail':[''],'type':['float'],'default':[1e-2]},
        'vdf.6' : {'comm':'type of input','cond':'MT3DRHOFLG==0','kw':['INDENSE'],
                 'detail':[['type','from previous step','Denseref','specif as dens','specif as conc']],
                'type':['choice'],'default':[0]},
        'vdf.7' : {'comm':'dens value','cond':'MT3DRHOFLG==0','kw':['DENSE'],
                 'detail':[''],'type':['float'],'default':[1000]},
        #VSC viscosity package for sewat
        'vsc.1' : {'comm':'species for viscosity','cond':'','kw':['MT3DMUFLG'],
                 'detail':[''],'type':['int'],'default':[-1]},
        'vsc.2' : {'comm':'viscosity limits','cond':'','kw':['VISCMIN','VISCMAX'],
                 'detail':['minimum','maximum'],'type':['float','float'],'default':[0.,0.]},
        'vsc.3a' : {'comm':'reference viscosity','cond':'','kw':['VISCREF'],
                 'detail':['ref'],'type':['float'],'default':[8.904e-4]},
        'vsc.3b' : {'comm':'used species','cond':'','kw':['NSMUEOS','MUTEMPOPT'],
                 'detail':['nb of species','calculation mode'],'type':['int','int'],'default':[1,1]},
        'vsc.3c' : {'comm':'ref species','cond':'','kw':['MTMUSPEC','DMUDC','CMUREF'],
                 'detail':['species nb','visc slope','ref conc.'],
                'type':['int','float','float'],'default':[1,1.923e-6,0.]},
        'vsc.3d' : {'comm':'coefficients','cond':'','kw':['MTMUTEMPSPEC','A1','A2','A3','A4'],
                 'detail':['temper. nb','coeff A1','coeff A2','coeff A3','coeff A4'],
                'type':['int','float','float','float','float'],
                'default':[2,239.4e-7,10.,248.37,133.15]},
        #UZT for Unsaturated PHT3D Ming Wu
        'uzt.1' : {'comm':'main flags','cond':'',
                   'kw':['MXUZCON','ICBCUZ','IETFLG','IDTORT','IFNAPL','NAPLNN'],
                 'detail':['?','?',['Evapotranspiration','no','yes'],['Tortuosity','no','yes'],
                           ['Napl layer moving','no','yes'],'Nb of NAPL species'],
                'type':['int','int','choice','choice','choice','int'],
                'default':[0,0,0,.5,0,0]},
        'uzt.2' : {'comm':'Concentration BC','cond':'','kw':['IUZFBND(NCOL,NROW)'],
                 'detail':[''],'type':['arrint'],'default':[0]},
        'uzt.3' : {'comm':'Initial water content','cond':'','kw':['SWC(NCOL,NROW)'],
                 'detail':[''],'type':['arrfloat'],'default':[0.1]},
        'uzt.4' : {'comm':'saturated thickness','cond':'','kw':['SDH(NCOL, NROW)'],
                 'detail':[''],'type':['arrfloat'],'default':[0.]},
        'uzt.5' : {'comm':'flag infiltration conc','cond':'','kw':['INCUZINF'],
                 'detail':[['present','no','yes']],'type':['choice'],'default':[0]},
        'uzt.6' : {'comm':'infiltration concentrations','cond':'INCUZINF>0',
                   'kw':['CUZINF(NCOL,NROW)'],
                 'detail':[''],'type':['float'],'default':[0.]},
        'uzt.7' : {'comm':'flag evapotransp. conc','cond':'','kw':['INCUZEVT'],
                 'detail':[['present','no','yes']],'type':['choice'],'default':[0]},
        'uzt.8' : {'comm':'evapotransp. concentrations','cond':'INCUZEVT>0',
                   'kw':['CUZEVT(NCOL,NROW)'],
                 'detail':[''],'type':['float'],'default':[0.]},
        'uzt.9' : {'comm':'flag gwet conc','cond':'','kw':['INCGWET'],
                 'detail':[['present','no','yes']],'type':['choice'],'default':[0]},
        'uzt.10' : {'comm':'gwet concentrations','cond':'INCGWET>0',
                   'kw':['CGWET(NCOL,NROW)'],
                 'detail':[''],'type':['float'],'default':[0.]},
                }
