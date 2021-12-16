class  Mtu:
    def __init__(self):
        self.grpList=['BCT','PCB','CRCH','CSS'] # EV 3/12/21
        bclist=['1a','1b'];
        bclist.extend([str(a) for a in range(2,21)])
        #bclist.extend([str(a) for a in range(14,21)])
        self.groups={
        'BCT':['bct.'+a for a in bclist],
        'PCB':['pcb.1','pcb.2'], # conc at Boundaries
        'CRCH':['crch.1'], # conc at recharge zones
        'CSS':['css.1'], # conc at source and sinks
        #'CGHB':['cghb.1'], # conc at g hed boundaries # EV 3/12/21
        #'DPT':['dpt.1'], # dual porosity
        #'DDF':['ddf.1'], # variable density
        }
        self.longNames={'BCT':'major transport parameters',
        'PCB': 'fixed concentrations',
        'CRCH':' Conc. in recharge',
        #'CRIV':' Conc. in rivers', # EV 3/12/21
        #'CCHD':' Conc. at transient head BC',
        'CSS': 'Conc at sources and sink'}
        
        self.lines={
            #BCT
        'bct.1a': {'comm':'Main flags','cond':'',
        'kw': ['ITRNSP','IBCTCB','MCOMP','ICBNDFLG','ITVD','IADSRB','ICT','CINACT','CICLOSE',\
        'IDISP','IXDISP','DIFFNC','IZOD','IFOD','IFMBC','IHEAT','NIMCOMP','IDISPCLN','NSEQITR',\
        'ITRNSP'],
        'detail' : [
        ['simulation','no transport','steady or each flow step','not implemented'], # -1 not set
        'save budget','nb of mobile species',
        ['domain','diff from flow','equal to flow'],
        'adv scheme, TVD nb',
        ['sorption','no sorption','linear','freundlich'],
        ['sorbed transport','no','total'],
        'Conc at inactive cells',
        'Solver Conc tolerance',
        ['Dispersion','no','isotropic','anisotropic'],
        ['cross dispersion','no','yes'],
        'Diffusion coeff',
        ['zero order decay','no','water','soil','soil+water'],
        ['1st order decay','no','water','soil','soil+water'],
        ['mass bal errors','no','computed'],
        'index for heat','nb of immobile species','CLN-GW equation','Nb sequential iterations',
        'Transport follows flow'
        ],
  
        'type':['choice','int','int','choice','int',
        'choice','choice','float','float','choice',
        'choice','float','choice','choice','choice',
        'int','int','int','int','int','int','int'
        ],
        'default':[1,0,1,1,0,0,0,-1.0,1e-9,1,0,1e-10,0,0,0,0,0,0,1,1]},  # OA 22/8/19 modif 2 to 1 for idisp

        'bct.1b': {'comm':'Print flags','cond':'IFMBC>0',
        'kw': ['MBEGWUNF','MBEGWUNT','MBECLNUNF', 'MBECLNUNT'],       
        'detail': [
        'unit nb for flow','unit nb transpt','unit nb CLN','unit nb imbalance'],
        'type':['int','int','int','int'],
        'default':[72,73,74,75]},
        # Boundary conditions
        'bct.2': {'comm':'Boundary conditions','cond':'ICBNDFLG==0',
        'kw': ['ICBUND'],'detail': ['transpt BC'],
        'type':['arrint'],'default':[1]},
        # porosity
        'bct.3': {'comm':'Porosity','cond':'',
        'kw': ['PRSITY'],'detail': ['Porosity'],
        'type':['arrfloat'],'default':[0.25]},
        # bulk density
        'bct.4': {'comm':'Bulk density','cond':'IADSRB>0',
        'kw': ['BULKD'],'detail': ['Bulk density'],
        'type':['arrfloat'],'default':[1.8]},
        # angle dispersivity
        'bct.5': {'comm':'dispersivity  angle','cond':'IDISP>0',
        'kw': ['ANGLEX'],'detail': ['Disp. angle'],
        'type':['float'],'default':[0.]},
        # long dispersivity
        'bct.6': {'comm':'long. dispersivity','cond':'IDISP==1',
        'kw': ['DL'],'detail': ['long. Disp.'],
        'type':['arrfloat'],'default':[1.]},
        # transverse dispersivity
        'bct.7': {'comm':'transv. dispersivity','cond':'IDISP==1',
        'kw': ['DT'],'detail': ['transv Disp.'],
        'type':['arrfloat'],'default':[0.1]},
        # x-direction longitudinal dispersivity for an anisotropic medium #EV 8/12/21
        'bct.8': {'comm':'x-direction long Disp.','cond':'IDISP==2',
        'kw': ['DLX'],'detail': ['x-direction long Disp.'],
        'type':['arrfloat'],'default':[1.]},
        # y-direction longitudinal dispersivity for an anisotropic medium
        'bct.9': {'comm':'y-direction long Disp.','cond':'IDISP==2',
        'kw': ['DLY'],'detail': ['y-direction long Disp. '],
        'type':['arrfloat'],'default':[1.]},
        # z-direction longitudinal dispersivity for an anisotropic medium
        'bct.10': {'comm':'z-direction long Disp.','cond':'IDISP==2',
        'kw': ['DLZ'],'detail': ['z-direction long Disp. '],
        'type':['arrfloat'],'default':[1.]},
        # xy-direction transverse dispersivity for an anisotropic medium
        'bct.11': {'comm':'xy-direction trans Disp.','cond':'IDISP==2',
        'kw': ['DTXY'],'detail': ['xy-direction trans Disp. '],
        'type':['arrfloat'],'default':[0.1]},
        # yz-direction transverse dispersivity for an anisotropic medium
        'bct.12': {'comm':'yz-direction trans Disp.','cond':'IDISP==2',
        'kw': ['DTYZ'],'detail': ['yz-direction trans Disp. '],
        'type':['arrfloat'],'default':[0.1]},
        # xz-direction transverse dispersivity for an anisotropic medium
        'bct.13': {'comm':'xz-direction trans Disp.','cond':'IDISP==2',
        'kw': ['DTXZ'],'detail': ['xz-direction trans Disp. '],
        'type':['arrfloat'],'default':[0.1]},
        # linear sorption
        'bct.14': {'comm':'Sorption', 'cond':'IADSRB!=0',
        'kw': ['ADSORB'],'detail': ['Kd'],
        'type':['arrfloat'],'default':[0.]},
        # freundlich
        'bct.15': {'comm':'Sorption 2nd', 'cond':'IADSRB==2',
        'kw': ['FLICH'],'detail': ['Freundlich'],
        'type':['arrfloat'],'default':[0.]},
        # zero order decay water
        'bct.16': {'comm':'zero order decay (W)', 'cond':'IZOD in [1,3]',
        'kw': ['ZODRW'],'detail': ['Oth coeff'],
        'type':['arrfloat'],'default':[0.]},
        # zero order decay soil
        'bct.17': {'comm':'zero order decay (S)', 'cond':'(IADSRB!=0) and (IZOD in [2,3])',
        'kw': ['ZODRS'],'detail': ['Oth coeff'],
        'type':['arrfloat'],'default':[0.]},
        # first order decay water
        'bct.18': {'comm':'1st order decay (W)', 'cond':'IFOD in [1,3]',
        'kw': ['FODRW'],'detail': ['1sto coeff'],
        'type':['arrfloat'],'default':[0.]},
        # first order decay soil
        'bct.19': {'comm':'1st order decay (S)', 'cond':'(IADSRB!=0) and (IFOD in [2,3])',
        'kw': ['FODRS'],'detail': ['1sto coeff'],
        'type':['arrfloat'],'default':[0.]},
        # initial concentrations (including T and immobile species
        'bct.20': {'comm':'initial concentrations', 'cond':'',
        'kw': ['CONC'],'detail': ['conc.'],
        'type':['arrfloat'],'default':[0.]},
        
        ## PCB
        # zones to specify concentrations
        'pcb.1': {'comm':'write budget', 'cond':'',
        'kw': ['IPCBCB'],'detail': ['print'],
        'type':['int'],'default':[0]},
        'pcb.2': {'comm':'conc. at boundaries', 'cond':'',
        'kw': ['CONC_BC'],'detail': ['conc.'],
        'type':['arrfloat'],'default':[0.]},

        ## CRCH
        'crch.1': {'comm':'Recharge concentrations', 'cond':'',
        'kw': ['RCHCONC'],'detail': [''],
        'type':['arrfloat'],'default':[0.]},
        ## CSS
        'css.1': {'comm':'source sink concentrations', 'cond':'',
        'kw': ['WELLCONC'],'detail': [''],
        'type':['arrfloat'],'default':[0.]},
        ## CRIV # EV 3/12/21
        #'criv.1': {'comm':'River concentrations', 'cond':'',
        #'kw': ['RIVCONC'],'detail': [''],
        #'type':['arrfloat'],'default':[0.]},
        ## CGHB
        #'cghb.1': {'comm':'Well concentrations', 'cond':'',
        #'kw': ['GHBCONC'],'detail': [''],
        #'type':['arrfloat'],'default':[0.]},
        ## DPT
        #'dpt.1': {'comm':'dummy', 'cond':'',
        #'kw': ['dummy1'],'detail': [''],
        #'type':['int'],'default':[0]},

        ## DDF
        #'ddf.1': {'comm':'dummy', 'cond':'',
        #'kw': ['dummy2'],'detail': [''],
        #'type':['int'],'default':[0]}, # EV 3/12/21
        }
