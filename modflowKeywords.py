class Mf:
    def __init__(self):
        self.grpList=['DIS','DISU','BAS6','BCF6','LPF','WEL','RCH','EVT','DRN','RIV','GHB','SIP','PCG','SOR','DE4',\
        'UPW','NWT','UZF','MNWT','SMS','GMG','HFB6'] # CHD written in auto mode
        self.groups={
        'DIS':['dis.'+str(i) for i in range(1,9)],
        'DISU':['disu.'+str(i) for i in range(1,10)],
        'BAS6':['bas.1','bas.2','bas.3','bas.4','bas.5'],
        'BCF6':['bcf.1','bcf.2','bcf.3','bcf.4','bcf.5','bcf.6','bcf.7','bcf.8','bcf.9'],
        'LPF':['lpf.'+str(i) for i in range(1,14)],
        'WEL':['wel.1'],
        'RCH':['rch.1','rch.2'],
        'EVT':['evt.1','evt.2','evt.3','evt.4'],
        'DRN':['drn.1'],
        'RIV':['riv.1'],
        'CHD':['chd.1'],
        'GHB':['ghb.1'],
        'SIP':['sip.1','sip.2'],
        'PCG':['pcg.1','pcg.2'],
        'SOR':['sor.1','sor.2'],
        'DE4':['de4.1','de4.2'],
        'MNWT':['mnwt.0','mnwt.1','mnwt.2a','mnwt.2b','mnwt.2c','mnwt.2d1','mnwt.2d2','mnwt.2f','mnwt.2g',
                'mnwt.2h','mnwt.3','mnwt.4a','mnwt.4b'],
        'UPW':['upw.'+str(i) for i in range(1,13)],
        'NWT':['nwt.1','nwt.2','nwt.3'],
        'UZF':['uzf.'+str(i) for i in range(1,17)],
        'SMS':['sms.1a','sms.1b','sms.2','sms.3'],
        'GMG':['gmg.1','gmg.2','gmg.3'],
        'HFB6':['hfb.3'], #'hfb.1','hfb.2',
        }
        # types
        # title : adds a # at the beginning
        # int or float : classical integer or float
        # vecint, vecfloat, arrint, arrfloat : vector or array of integers or floats
        # vecstring : a vector of strings
        # choice : choice list
        # layint layfloat : a vector of integers for layers that a re-written in lnes of 40 values
        self.lines={
        #DIS
        'dis.1':{'comm':'Domain','cond':'','kw':['MFDOMN'],'detail':[],'type':['arrfloat'],'default':['#']},
        'dis.2':{'comm':'Model properties','cond':'','kw':['NLAY','NROW','NCOL','NPER','ITMUNI','LENUNI'],
               'detail':['Nb of layers','Nb of rows','Nb of columns','Nb of periods',['time units','-','sec','min','hours','days','years'],
                ['length units','-','ft','m','cm']],'type':['int','int','int','int','choice','choice'],
                'default':[1,10,10,1,4,2]},
        'dis.3':{'comm':'Quasi_3D confining bed','cond':'','kw':['LAYCBD(NLAY)'],'detail':[], #['Quasi_3D confining bed','No','Yes']],
                'type':['layint'],'default':[0]},
        'dis.4':{'comm':'Col width','cond':'','kw':['DELR(NCOL)'],'detail':['Col width'],
                        'type':['vecfloat'],'default':[10],'units':['L']},
        'dis.5':{'comm':'Row height','cond':'','kw':['DELC(NROW)'],'detail':['Row height'],
                        'type':['vecfloat'],'default':[10],'units':['L']},
        'dis.6':{'comm':'Top','cond':'','kw':['TOP(NROW,NCOL)'],'detail':[],
                        'type':['arrfloat'],'default':[10.],'units':['L']},
        'dis.7':{'comm':'Bottom','cond':'','kw':['BOTM(NLAY,NROW,NCOL)'],'detail':[],
                        'type':['arrfloat'],'default':[0.],'units':['L']},
        'dis.8':{'comm':'Periods characteristics','cond':'','kw':['PERLEN','NSTP','TSMULT','SsTr'],
                'detail':['Period length','internal steps','multiplier',['type','Steady state','transient']],
                'type':['float','int','float','choice'],
                'default':[1.,1,1.,0],'units':['T','','','']},
        #DISU for modflow UNS
        'disu.1':{'comm':'Title','cond':'','kw':['title'],'detail':[],'type':['title'],'default':['#']},
        'disu.2':{'comm':'Type of mesh','cond':'','kw':['MshType'],
                'detail':[['','Rectangle','Nested','Triangle','Voronoi']],
                'type':['choice'],'default':[0]},
        'disu.3':{'comm':'Mesh','cond':'','kw':['U_Mesh'],'detail':[],
                  'type':['arrfloat'],'default':[0.]},
        'disu.4':{'comm':'Model properties','cond':'','kw':['NCELL','UNLAY','NJAG','IVSD','UNPER','ITMUNI','LENUNI'],
               'detail':['Nb of cells','Nb of layers','nb of connections','Vert. subdiscretiz.','Nb of periods',['time units','-','sec','min','hours','days','years'],
                ['length units','-','ft','m','cm']],'type':['int','int','int','int''int','int','choice','choice'],
                'default':[1,1,0,-1,1,4,2]},  # OA 19/8/19
        'disu.5':{'comm':'type of layer','cond':'','kw':['LAYCBD(NLAY)'],'detail':[['confining bed below','No','Yes']],
                'type':['layint'],'default':[0]},
        'disu.6':{'comm':'nb of nodes per layer','cond':'','kw':['NODELAY(NLAY)'],'detail':['confining bed below'],
                'type':['arrint'],'default':[0]},
        'disu.7':{'comm':'Top','cond':'','kw':['TOP(NROW,NCOL)'],'detail':[],
                        'type':['arrfloat'],'default':[10.],'units':['L']},
        'disu.8':{'comm':'Bottom','cond':'','kw':['BOTM(NLAY,NROW,NCOL)'],'detail':[],
                        'type':['arrfloat'],'default':[0.],'units':['L']},
        'disu.9':{'comm':'Periods characteristics','cond':'','kw':['PERLENu','NSTPu','TSMULTu', 'SsTru'], # OA 3/11/18 
                'detail':['Period length','internal steps','multiplier',['type','Steady state','transient']],
                'type':['float','int','float','choice'],
                'default':[1.,1,1.,0],'units':['T','','','']},
        #BA6
        'bas.1':{'comm':'Bas title 1','cond':'','kw':['title'],'detail':[],'type':['title']},
        'bas.2':{'comm':'Bas options','cond':'','kw':['bas_opt'],'detail':[],
                'type':[''],'default':['FREE']},  # OA 10/2/20 changed to free
        'bas.3':{'comm':'Boundary conditions','cond':'','kw':['IBOUND(NLAY,NROW,NCOL)'],'detail':[],
                        'type':['arrint'],'default':[1]},
        'bas.4':{'comm':'value of head for no flow','cond':'','kw':['HNOFLO'],'detail':[],
                 'type':['float']},
        'bas.5':{'comm':'Initial and fixed head','cond':'','kw':['STRT(NLAY,NROW,NCOL)'],'detail':[],
                'type':['arrfloat'],'units':['L']},
        
        #BCF
        'bcf.1':{'comm':'General flags','cond':'','kw':['IBCFCB','HDRY','IWDFLG','WETFCT','IWETIT','IHDWET'],
               'detail':[['write budget','no','yes'],'head for dry cells',['wetting','inactive','active'],
                         'factor to wet','iteration interval for wetting',['wetting equation','eq1','eq2']],
                'type':['int','float','int','float','int','int']},
        'bcf.2':{'comm':'Type of layer (average/confinement)','cond':'','kw':['LAYCON(NLAY)'],'detail':[],'type':['layint']},
        'bcf.3':{'comm':'anisotropy factor','cond':'','kw':['TRPY(NLAY)'],'detail':[],'type':['vecint']},
        'bcf.4':{'comm':'First storage coeff','cond':'NPER>1','kw':['Sf1(NCOL,NROW)'],'detail':[],'type':['arrfloat']},
        'bcf.5':{'comm':'Transmissivity','cond':'LAYCON in [0,2]','kw':['Tran(NCOL,NROW)'],'detail':[],'type':['arrfloat']},
        'bcf.6':{'comm':'Hydraulic conductivity','cond':'LAYCON in [1,3]','kw':['HY(NCOL,NROW)'],'detail':[],'type':['arrfloat']},
        'bcf.7':{'comm':'Vertical conductance','cond':'NLAY>1','kw':['Vcont(NCOL,NROW)'],'detail':[],'type':['arrfloat']},
        'bcf.8':{'comm':'2nd storage coeff','cond':'NPER>1 and LAYCON>=2','kw':['Sf2(NCOL,NROW)'],'detail':[],'type':['arrfloat']},
        'bcf.9':{'comm':'','cond':'IWDFLG!=0 and LAYCON in [1,3]','kw':['WETDRY(NCOL,NROW)'],'detail':[],'type':['arrfloat']},
        #LPF
        'lpf.0':{'comm':'Lpf title 1','cond':'','kw':['title'],'detail':[],'type':['title']},
        'lpf.1':{'comm':'General flags','cond':'','kw':['ILPFCB','HDRY','NPLPF','NWOPT1','NWOPT2','NWOPT3'],
               'detail':['Write budget','head for dry cells','Nb Lpf params','Newt opt1','Newt opt2','Newt opt3'],#EV 25/02/20 + OA 7/5/21
               'type':['int','float','int','int','string','string'],'default':[31,0.01,0,0,' ',' ']},
        'lpf.2':{'comm':'Type of layer: confined(0) or convertible(1)','cond':'','kw':['LAYTYP(NLAY)'],'detail':[],'type':['layint']}, #[['confinement','confined','convertible']],
        #'lpf.2':{'comm':'Type of layer: confined(0) or convertible(1)','cond':'','kw':['LAYTYP(NLAY)'],'detail':['Confined or convertible'],'type':['arrint'],'default':[0]},
        'lpf.3':{'comm':'Type of averaging','cond':'','kw':['LAYAVG(NLAY)'],'detail':[['type of average','harmonic','logarithmic','harm/log']],'type':['layint']},
        'lpf.4':{'comm':'Horizontal anisotropy','cond':'','kw':['CHANI(NLAY)'],'detail':[],
                    'type':['layfloat'],'default':[' 1']},
        'lpf.5':{'comm':'flag for vertical cond: VKA(0) or ratio(1)','cond':'','kw':['LAYVKA(NLAY)'],
                'detail':[['what is Vka','vertical K','ratio Kh/Kv']],'type':['layint']},
        'lpf.6':{'comm':'wetting active or inactive','cond':'','kw':['LAYWET(NLAY)'],
                'detail':[['wetting','inactive','active']],'type':['layint']},
        'lpf.7':{'comm':'calculation of wetting','cond':'LAYWET>0','kw':['WETFCT','IWETIT','IHDWET'],'detail':[],
               'type':['float','int','int'],'default':[0.1,1,1]},
        'lpf.8':{'comm':'Hydraulic conductivity','cond':'','kw':['HK(NLAY,NROW,NCOL)'],
                'detail':[],'type':['arrfloat'],'units':['L/T'],'default':[1]},
        'lpf.9':{'comm':'Vertical cond. or ratio Kh/Kv','cond':'','kw':['VKA(NLAY,NROW,NCOL)'],
                'detail':[],'type':['arrfloat'],'units':['L/T'],'default':[1]},
        'lpf.10':{'comm':'First storage coeff','cond':'NPER>1','kw':['Sf1(NLAY,NROW,NCOL)'],
                    'detail':[],'type':['arrfloat'],'default':[1e-6]},
        'lpf.11':{'comm':'2nd storage coeff','cond':'NPER>1','kw':['Sf2(NLAY,NROW,NCOL)'], # EV 25/10/18 removed cond NPER>1 and LAYTYP>0 EV 03/05/19 added cond NPER>1 
                    'detail':[],'type':['arrfloat'],'default':[.25]},
        'lpf.12':{'comm':'Vertical cond confined bed','cond':'LAYCBD>0','kw':['VKCB(NLAY,NROW,NCOL)'],
                    'detail':[],'type':['arrfloat']},
        'lpf.13':{'comm':'Threshold wetting and neighbours','cond':'LAYWET!=0 and LAYTYP!=0',
                  'kw':['WETDRY(NLAY,NROW,NCOL)'],'detail':[],'type':['arrfloat'],
                    'default':[0.1]},
        #WEL
        'wel.1': {'comm':'zones for Wells','cond':'','kw':['WELLS(NLAY,NROW,NCOL)'],
                'detail':[],'type':['arrfloat'],'units':['L3/T']},
        #'wel.1':{'comm':'Flags for Wells','cond':'','kw':['MXACTW','IWELCB'],'detail':['total nb of wells','save budget'],'type':['','']},
        #'wel.2':{'comm':'Wells','cond':'','kw':['WLayer','WRow','WColumn','WQ'],
        #       'detail':['well layer','well row','well column','wel discharge'],'type':['','','','']},
        #RECH
        'rch.1':{'comm':'Flags for recharge','cond':'','kw':['NRCHOP','IRCHCB'],
                 'detail':[['position of recharge','-','top cell','not implemented','highest active'],'save to budget'],
                'type':['choice','int'],'default':[3,0]}, #EV 22/07/2019
        'rch.2':{'comm':'Recharge value','cond':'','kw':['RECH(NPER,NROW,NCOL)'],'detail':[],
                    'type':['arrfloat'],'units':['L/T']},
        #EVT
        'evt.1':{'comm':'Flags for evapotransp','cond':'','kw':['NEVTOP','IEVTCB'],
                 'detail':[['position of evapo','-','top cell','specified'],'save to budget'],
                'type':['',''],'default':[1,0]},
        'evt.2':{'comm':'Evapo. ref. surface','cond':'','kw':['SURF(NPER,NROW,NCOL)'],'detail':[],
                    'type':['arrfloat'],'units':['L']},
        'evt.3':{'comm':'Evapotransp value','cond':'','kw':['EVTR(NPER,NROW,NCOL)'],'detail':[],
                    'type':['arrfloat'],'units':['L/T']},
        'evt.4':{'comm':'Evapo depth','cond':'','kw':['EXPD(NPER,NROW,NCOL)'],'detail':[],
                    'type':['arrfloat'],'units':['L'],'default':[0.5]},
        #Solvers
        'sip.1':{'comm':'SIP globals','cond':'','kw':['MXITER','NPARM'],'detail':['Mx nb of iteration','Nb of iteration variables'],
                 'type':['',''],'default':[500,5]},
        'sip.2':{'comm':'Steps','cond':'','kw':['SACCL','SHCLOSE','IPCALC','WSEED','IPRSIP'],
               'detail':['Acceleration','tolerance for head',['seed','given','calculated'],'seed value','print interval'],
               'type':['float','float','choice','float','int'],
               'default':[1.,1e-9,0,1e-3,1]},
        'pcg.1':{'comm':'PCG globals','cond':'','kw':['MXITER','ITER1','NPCOND'],
               'detail':['Mx nb of iteration','inner iterations',['solution','-','Cholesky','Polynomial']],
                'type':['int','int','choice'],'default':[10,30000,1]},
        'pcg.2':{'comm':'Steps','cond':'','kw':['PHCLOSE','RCLOSE','RELAX','NBPOL','IPRPCG','MUTPCG','DAMP'],
               'detail':['Head tolerance','residual volume convergence','Relaxation',['upper bound','-','calculated','=2'],
                         'print interval',['print','all tables','nb of iterations','no print','if fail'],'Damping factor'],
               'type':['float','float','float','choice','int','choice','float'],
                'default':[1e-4,1e-3,1.,2,0,2,1.]},
        'sor.1':{'comm':'SOR globals','cond':'','kw':['MXITER'],'detail':['Mx nb of iteration'],
                 'type':['int'],'default':[1000]},
        'sor.2':{'comm':'SOR globals','cond':'','kw':['OACCL','OHCLOSE','IPRSOR'],
               'detail':['Accelaeration','head tolerance','print interval'],
               'type':['float','float','int'],'default':[1.,1e-9,0]},
        'de4.1':{'comm':'DE4 globals','cond':'','kw':['ITMX','MXUP','MXLOW','MXBW'],
               'detail':['Mx nb of iteration','Mx eq upper','Mx eq lower','Mx eq Bwidth'],
               'type':['int','int','int','int'],'default':[1,0,0,0]},
        'de4.2':{'comm':'Steps','cond':'','kw':['IFREQ','MUTD4','DACCL','DHCLOSE','IPRD4'],
                'detail':['Flag for Freq','Flag for print','acceleration','head tolerance','time step to print'],
               'type':['int','int','float','float','int'],'default':[1,2,1.,1e-9,5]},
        'gmg.1':{'comm':'GMG globals','cond':'','kw':['G_RCLOSE','G_IITER','G_HCLOSE','G_MXITER'],
               'detail':['relative tol','inner iterations','head tol','Mx nb of iteration'],
                'type':['float','int','float','int'],'default':[1e-5,100,1e-5,50]},
        'gmg.2':{'comm':'Damping','cond':'','kw':['G_DAMP','G_IADAMP','G_IOUTGMG'],
               'detail':['Damping factor','adaptative damping','output'],
                'type':['float','int','int'],'default':[.5,1,0]},
        'gmg.3':{'comm':'Steps','cond':'','kw':['G_ISM','G_ISC','G_DUP','G_DLOW','CHGLIMIT'],
               'detail':[['smoother','ILU','SGS'],'semicoarsening','upper damp','lower damp','max outer head change'],
               'type':['choice','int','float','float','float'],
                'default':[0,0,1.,0.75,0.001,0.01]},
        'gmg.4':{'comm':'relaxation','cond':'G_ISC==4','kw':['G_RELAX'],
               'detail':['relaxation factor'],'type':['float'],'default':[1.]},
        #SMS sovler for modflow USG
        #1a. OPTIONS SIMPLE MODERATE COMPLEX
        'sms.1a':{'comm':'SMS options','cond':'','kw':['SMS_opt'],
               'detail':['options','No opt','SIMPLE','MODERATE','COMPLEX'],
               'type':['choice'],'default':[0]},
        'sms.1b':{'comm':'SMS globals','cond':'',
                'kw':['sHCLOSE','sHICLOSE','sMXITER','sITER1','sIPRSMS','sNONLINMTH','sLINMTH'],
               'detail':['H convergence (inner)','Hconvergence(outer)','Max. iterations (outer)','Max. iterations (inner)',
                    'Print converg.','Nonlinear method',['Non lin solver','-','xMD','PCGU']],
               'type':['float','float','int','int','int','int','choice'],
               'default':[1e-6,1e-15,10,10000,1,0,1]},
        #If NONLINMTH != 0 and OPTIONS is not specified then read item 2
        'sms.2':{'comm':'SMS globals','cond':'sNONLINMTH!=0 and SMS_opt==0',
                'kw':['sTHETA','sAKAPPA','sGAMMA','sAMOMENTUM','sNUMTRACK','sBTOL','sBREDUC','sRESLIM'],
                'detail':['Reduc fact DbD','Increment DbD','Memory term','Fraction of memory','Max bactracking',
                          'Residual reduction tol','Reuction step size','limit for resid backtrack'],
               'type':['float','float','float','float','int','float','float','int'],
               'default':[0.7,0.1,0.2,0.001,10,1e4,0.2,100]},
        #If LINMTH = 1 and OPTIONS is not specified then read item 3 for the XMD solver
        'sms.3':{'comm':'XMD solver','cond':'sLINMTH==1 and SMS_opt==0',
                 'kw':['sIACL','sNORDER','sLEVEL','sNORTH','sIREDSYS','sRRCTOL','sIDROPTOL','sEPSRN'],
                'detail':[['Accel method','Conjug grad','Orthomin','BicGstab'],['ordering','original','reverse','min degree'],
                          'level of fill','Nb accel for orthomin',['reduc system','no','red_black'],
                        'Residual tolerance',['drop tolerance','no','yes'],'drop tolerance value'],
               'type':['choice','choice','int','int','choice','float','choice','float'],
               'default':[2,0,0,7,0,0.,0,1e-3]},
        #If LINMETH = 2 and OPTIONS is not specified then read item 4 for the PCGU solver
        #4. [CLIN] IPC ISCL IORD RCLOSEPCGU [RELAXPCGU]
        # Multinode well
        # general
        'mnwt.0':{'comm':'Title','cond':'','kw':['title'],'detail':[],'type':['title'],'default':['#']},
        'mnwt.1':{'comm':'General flags','cond':'','kw':['MNWMAX','IWL2CB','MNWPRNT','OPTION'],
                  'detail':['Max nb of MN wells','units or print','Level of printing','option'],
                  'type':['int','int','int','int'],'default':[0,0,0,0]},
        # for each well
        'mnwt.2a':{'comm':'name and nodes','cond':'','kw':['WELLID','NNODES'],
                  'detail':['well name','nb of nodes (> or < 0)'],
               'type':['arrint','int'],'default':[0,0]},
        'mnwt.2b':{'comm':'model and options','cond':'','kw':['LOSSTYPE','PUMPLOC','Qlimit','PPFLAG','PUMPCAP'],
                  'detail':[['loss model','NONE','THIEM','SKIN','GENERAL','SPECIFYcwc'],
                            'specif pump position','limit for Q?','Adjust for partial pen.','Adjust Q?'],
               'type':['choice','int','int','int','int'],'default':[0,0,0,0,0]},
        'mnwt.2c':{'comm':'well characteristics','cond':'LOSSTYPE>0',
                   'kw':['Rw','Rskin','Kskin','Bparm','Cparm','Pparm','CWC'],
                  'detail':['well radius','skin radius','skin K','Parm B loss','Parm C loss','Power P loss','CellWell conductance'],
               'type':['float','float','float','float','float','float','float'],
                'default':[0.1,0.15,1e-4,1.,1.,2.,0.1]},
        ## for each node
        'mnwt.2d1':{'comm':'well node characteristics 1','cond':'NNODES>0',
                   'kw':['Lay','Row','Col','Rw','Rskin','Kskin','Bparm','Cparm','Pparm','CWC','PP'],
                  'detail':['Layer','Row','Column','well radius','skin radius','skin K','Parm B loss','Parm C loss','Power P loss','CellWell conductance','fract penetration'],
               'type':['int','int','int','float','float','float','float','float','float','float'],
                'default':[1,1,1,0.1,0.15,1e-4,1.,1.,2.,0.1,0.0]},
        'mnwt.2d2':{'comm':'well node characteristics 2','cond':'NNODES<0',
                   'kw':['Ztop','Zbot','Row','Col','Rw','Rskin','Kskin','Bparm','Cparm','Pparm','CWC','PP'],
                  'detail':['top','bottm','Row','Column','well radius','skin radius','skin K','Parm B loss','Parm C loss','Power P loss','CellWell conductance','fract penetration'],
               'type':['float','float','int','int','float','float','float','float','float','float','float'],
                'default':[10.,0.,1,1,0.1,0.15,1e-4,1.,1.,2.,0.1,0.0]},
        ## for complex limiting conditions
        'mnwt.2e':{'comm':'pump position','cond':'PUMPLOC!=0','kw':['PUMPLAY','PUMPROW','PUMPCOL','Zpump'],
                  'detail':['layer of pump','pump row','pump col','pump elevation'],
               'type':['int','int','int','float'],'default':[0,0,0,0.]},
        'mnwt.2f':{'comm':'limiting Q','cond':'Qlimit>0','kw':['Hlim','QCUT','Qfrcmn','Qfrcmx'],
                  'detail':['Limiting head','Specify as rate or fraction','min. rate','max. rate'],
               'type':['float','int','float','float'],'default':[10.,0,0.,1.]},
        'mnwt.2g':{'comm':'adjust for lift','cond':'PUMPCAP>0','kw':['Hlift','LIFTq0','LIFTqmax','HWtol'],
                  'detail':['max head for lift','min lift','Head diff. in calc.'],
               'type':['float','float','float','float'],'default':[10.,0.,1.,0.01]},
        'mnwt.2h':{'comm':'lift curve','cond':'PUMPCAP>0','kw':['LIFTn','Qn'],
                  'detail':['lift n','Q n'],'type':['float','float'],'default':[0.,0.]},
        ## for each period
        'mnwt.3':{'comm':'nb of wells for period','cond':'','kw':['ITMP'],
                  'detail':['i period'],'type':['int'],'default':[0]},
        'mnwt.4a':{'comm':'well rate','cond':'ITMP>0','kw':['WELLID','Qdes','CapMult','Cprime'],
                  'detail':['name of the well','Desired Q','curve Q lift','Concentration'],
                  'type':['string','float','float','float'],'default':['',0.,1.,0.]},
        'mnwt.4b':{'comm':'well conditions','cond':'Qlimit<0','kw':['Hlim','QCUT','Qfrcmn','Qfrcmx'],
                  'detail':['Limiting head','Specify as rate or fraction','min. rate','max. rate'],
                  'type':['float','int','float','float'],'default':[10.,0,0.,1.]},
        
        #Drain package
        'drn.1':{'comm':'Drains: elevation, conductance','cond':'',
                'kw':['DRAIN(NLAY,NROW,NCOL)','DR_COND(NLAY,NROW,NCOL)'],'detail':[],
                'names':['Elevation','Conductance'],
                'type':['arrfloat','arrfloat'],'default':[0,0]},
        #River package
        'riv.1':{'comm':'River: stage, conduct., Rbot','cond':'','kw':['RIV_Stage(NLAY,NROW,NCOL)','RIV_Cond(NLAY,NROW,NCOL)','RIV_Botm(NLAY,NROW,NCOL)'],'detail':[],
                 'names':['Stage','Conductance','Elevation'], # added OA 25/4/19
               'type':['arrfloat','arrfloat','arrfloat'],'default':[0,0,0]}, # modif OA 25/4/19
        #Variable head package
        'chd.1':{'comm':'Heads: first, last','cond':'','kw':['CHD(NLAY,NROW,NCOL)'],'detail':[],
               'type':['arrfloat'],'default':[0]},
        #Ghb package
        'ghb.1':{'comm':'General Head: head, conductance','cond':'',
                'kw':['GHB_HD(NLAY,NROW,NCOL)','GHB_COND(NLAY,NROW,NCOL)'],'detail':[],
                'names':['Elevation','Conductance'],
               'type':['arrfloat','arrfloat'],'default':[0,0]},
        #UPW (close to Lpf for modflow nwt
        'upw.1':{'comm':'General flags','cond':'','kw':['IUPWCB','HDRY','NPUPW','IPHDRY'],
               'detail':['write budget','head for dry cells','Nb upw params',['set to dry','no','yes']],
               'type':['int','float','int','choice'],
                'default':[0,-999,0,1]},
        'upw.2':{'comm':'Type of layer (confinement)','cond':'','kw':['UlAYTYP(NLAY)'],'detail':[['confinement','confined','convertible']],'type':['layint']},
        'upw.3':{'comm':'Type of averaging','cond':'','kw':['UlAYAVG(NLAY)'],'detail':[['type of average','harmonic','logarithmic','arith/log']],'type':['layint']},
        'upw.4':{'comm':'anisotropy flag','cond':'','kw':['UcHANI(NLAY)'],'detail':[],'type':['layfloat']},
        'upw.5':{'comm':'flag for vertical cond','cond':'','kw':['UlAYVKA(NLAY)'],'detail':[['what is Vka','vertical K','ratio Kv/Kh']],'type':['layint']},
        'upw.6':{'comm':'wetting active or inactive','cond':'','kw':['UlAYWET(NLAY)'],'detail':[['wetting','inactive','active']],'type':['layint']},
        'upw.7':{'comm':'Hydraulic conductivity','cond':'','kw':['UhKV(NLAY,NROW,NCOL)'],'detail':[],'type':['arrfloat']},
        'upw.8':{'comm':'anisotropy coeff','cond':'UcHANI<=0','kw':['UhANI(NLAY,NROW,NCOL)'],'detail':[],'type':['arrfloat']}, #UCHANI<0
        'upw.9':{'comm':'Vertical conductivity','cond':'','kw':['UvKA(NLAY,NROW,NCOL)'],'detail':[],'type':['arrfloat'],'default':[10]},
        'upw.10':{'comm':'First storage coeff (Ss)','cond':'NPER>1','kw':['Uss(NLAY,NROW,NCOL)'],
                    'detail':[],'type':['arrfloat'],'default':[1e-7]},
        'upw.11':{'comm':'2nd storage coeff (Sy)','cond':'NPER>1 and UlAYTYP>0','kw':['Usy(NLAY,NROW,NCOL)'],
                    'detail':[],'type':['arrfloat'],'default':[.4]},
        'upw.12':{'comm':'Vertical conductivity confined bed','cond':'LAYCBD>0',
                  'kw':['UvKCB(NLAY,NROW,NCOL)'],'detail':[],'type':['arrfloat']},
        #NWT solver keywords
        'nwt.1':{'comm':'General flags','cond':'',
                 'kw':['HEADTOL','FLUXTOL','MAXITEROUT','THICKFACT','LINMETH','IPRNWT','IBOTAV','N_OPTS',
                      'DBDTHETA','DBDKAPPA','DBDGAMMA','MOMFACT','BACKFLAG','MAXBACKITER',
                      'BACKTOL','BACKREDUCE'],
                 'detail':['max head change','max flux change','max outer iterations','portion to adjust storage',
                    'mat solver 1 gmres, 2 xmd',['print solver convergence','no','yes'],
                    ['correction for dewaterd neighbours','no','yes'],'type of problem',
                    'coeff reduce head change','coeff incr. head change','weight factor','momentum coeff',
                    'flag residual control','max nb of reduction','resid prop decrease','resid factor'],
                 'type':['float','float','int','float','int','int','choice','string',
                         'float','float','float','float','int','int','float','float'],
                'default':[1e-4,500,100,1e-5,1,0,0,'SPECIFIED',
                           0.7,1e-4,0.,0.1,0,20,2.,.6]},
        'nwt.2':{'comm':'flags for GMRES','cond':'LINMETH==1 and N_OPTS==\'SPECIFIED\'',
                 'kw':['MAXITINNER','ILUMETHOD','LEVFILL','STOPTOL','MSDR'],
                 'detail':['max nb of iteration',['LU solver','-','drop tolerance','order k'],
                 'fill limit','tol for LU solver','nb iteration for GMRES'],
                 'type':['int','choice','int','float','int'],
                 'default':[500,2,1,1e-10,10]},
        'nwt.3':{'comm':'flags for XMD','cond':'LINMETH==2 and N_OPTS==\'SPECIFIED\'',
                 'kw':['IACL','NORDER','LEVEL','NORTH','IREDSYS','RRCTOLS','IDROPTOL','EPSRN',
                       'HCLOSEXMD','MXITERXMD'],
                 'detail':[['accel method','conj grad','orthomin','GBstab'],['ordering','original','RCM','min degree'],
                       'level of fill','nb of orthogonal',['apply reduced system','no','yes'],
                        'residual tolerance',['use drop reduction','no','yes'],
                        'drop tolerance','head closure criteria','mx nb iteration'],
                 'type':['choice','choice','int','int','choice','float','choice','float','float','int'],
                 'default':[2,1,1,2,0,0.,1,1e-3,1e-4,50]},
        
        #UZF package for unsaturated zone
        'uzf.1':{'comm':'General flags','cond':'',
                 'kw':['NUZTOP','IUZFOPT','IRUNFLG','IETFLG','IUZFCB1','IUZFCB2','NTRAIL2','NSETS2','NUZGAG','SURFDEP4'],
                 'detail':[['recharge cell','-','top','spec. layer','upper wet'],
                ['unsat K','no unsat','from VKS','from LPF'],['budget for gw to surface','removed','routed'],
                ['Evapotranspiration','not simulated','simuated'],
                'save budget 1','save budget 2','nb of trailing waves','nb of waves for transient',
                'nb of cells to print','undulation depth'],
                 'type':['choice','choice','choice','choice','int','int','int','int','int','float'],
                 'default':[0,0,0,0,0,0,15,20,0,.5]},
        'uzf.2':{'comm':'aerial extend of simulation','cond':'','kw':['IUZFBND(NCOL,NROW)'],
                'detail':[],'type':['arrint'],'default':[1]},
        'uzf.3':{'comm':'position of segments for streamflow','cond':'IRUNFLG>0',
                 'kw':['IRUNBND(NCOL,NROW)'],'detail':[],'type':['arrint']},
        'uzf.4':{'comm':'vertical hydr conductivity','cond':'IUZFOPT==1',  # OA 7/3/19
                 'kw':['VKS(NCOL,NROW)'],'detail':[],'type':['arrfloat'],'default':[10.]},
        'uzf.5':{'comm':'epsilon brooks-corey','cond':'IUZFOPT>0',
                 'kw':['EPS(NCOL,NROW)'],'detail':[],'type':['arrfloat'],'default':[2.5]},
        'uzf.6':{'comm':'water content at saturation','cond':'IUZFOPT>0',
                 'kw':['THTS(NCOL,NROW)'],'detail':[],'type':['arrfloat'],'default':[0.4]},
        'uzf.7':{'comm':'initial water content','cond':'IUZFOPT>0',
                 'kw':['THTI(NCOL,NROW)'],'detail':[],'type':['arrfloat'],'default':[0.2]},
        'uzf.8':{'comm':'observation points','cond':'NUZGAG>0',
                 'kw':['IUZROW','IUZCOL','IFTUNIT','IUZOPT'],
                'detail':['row number','col number','unit to print on','options'],
                'type':['int','int','int','int'],'default':[0,0,0,0]},
        'uzf.9':{'comm':'specify (0) or reuse (-1) infiltration','cond':'',
                 'kw':['NUZF1'],'detail':[],'type':['int'],'default':[0]},
        'uzf.10':{'comm':'infiltration rate','cond':'NUZF1>-1',
                 'kw':['FINF(NPER,NCOL,NROW)'],'detail':[],'type':['arrfloat'],'default':[0]},
        'uzf.11':{'comm':'specify (0) or reuse (-1) ETP','cond':'IETFLG>0',
                 'kw':['NUZF2'],'detail':[],'type':['int'],'default':[0]},
        'uzf.12':{'comm':'evapotranspiration rate','cond':'IETFLG>0 and NUZF2>-1',
                 'kw':['PET(NPER,NCOL,NROW)'],'detail':[],'type':['arrfloat'],'default':[0]},
        'uzf.13':{'comm':'specify (0) or reuse (-1) extinc depth','cond':'IETFLG>0',
                 'kw':['NUZF3'],'detail':[],'type':['int'],'default':[-1]},
        'uzf.14':{'comm':'extinction depth','cond':'IETFLG>0 and NUZF3>-1',
                 'kw':['EXTDP(NPER,NCOL,NROW)'],'detail':[],'type':['arrfloat'],'default':[0]},
        'uzf.15':{'comm':'specify or reuse extinc wat content','cond':'IETFLG>0',
                 'kw':['NUZF4'],'detail':[],'type':['int'],'default':[-1]},
        'uzf.16':{'comm':'extinction water content','cond':'IETFLG>0 and NUZF4>-1',
                 'kw':['EXTWC(NPER,NCOL,NROW)'],'detail':[],'type':['arrfloat'],'default':[0]},
                  
        # HFB for Horizontal-Flow-Barrier package
        #'hfb.1':{'comm':'Nb of HFB in the model','cond':'','kw':['NHFB'],'detail':[],
              # 'type':['int'],'default':[0]},
        #'hfb.2':{'comm':'Nb of HFB in a layer','cond':'','kw':['NBRLAY'],'detail':[],
               #'type':['int'],'default':[0]},
        'hfb.3':{'comm':'Hydraulic characteristic Kh/e','cond':'','kw':['HYDCHR(IROW1,ICOL1,IROW2,ICOL2)'],'detail':[],
               'type':['arrint'],'default':[0]}, #NBRLAY>0
        }