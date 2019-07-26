class Pst:
    def __init__(self):
        self.grpList=['SYS','CTD','SVD','PGR','OBS','REG','PPP'] 
        self.groups={
        'SYS':['sys.1','sys.2'],    
        'CTD':['ctd.1','ctd.2','ctd.3','ctd.4','ctd.5','ctd.6'],
        'SVD':['svd.1','svd.2','svd.3'],
        'PGR':['pgr.1','pgr.2','pgr.3','pgr.4','pgr.5','pgr.6'],
        'OBS':['obs.1'],
        'REG':['reg.1','reg.2','reg.3'],
        'PPP':['ppp.1','ppp.2','ppp.3'],
        }
        self.lines={
        #SYS
        'sys.1':{'comm':'system','cond':'','kw':['SYSTEM'],
                 'detail':[['operating system','windows','linux']],'type':
                 ['choice'],'default':[0]},
        'sys.2':{'comm':'system','cond':'SYSTEM > 0','kw':['ORTI3D_path'],
                 'detail':['ORTI3D path on linux'],'type':
                 ['string'],'default':['/home/ORTI3D']},
        #CTD
        'ctd.1':{'comm':'Control data l1','cond':'','kw':['RSTFLE','PESTMODE'],
                 'detail':[['write restart data','restart','norestart'],
                 ['pest mode','estimation','not_implemented','regularisation']],
                 'type':['choice','choice'],'default':[0,0]},
        'ctd.2':{'comm':'Control data l2','cond':'','kw':['NPAR'],
               'detail':['nb of parameters'],'type':['int']},
        'ctd.3':{'comm':'Control data l3','cond':'','kw':['NUMCOM'],
               'detail':['nb of model command'],'type':['int'],'default':[1]},
        'ctd.4':{'comm':'Control data l5','cond':'','kw':['RELPARMAX','FACPARMAX',
                'FACORIG'],'detail':['param relative chg limit', 'param factor chg limit'
                , 'minimum fraction'],'type':['float','float','float'],
                 'default':[5.0, 5.0, 0.001]},
        'ctd.5':{'comm':'Control data l6','cond':'','kw':['PHIREDSWH'],
               'detail':['obj fct chg'],'type':['float'],'default':[0.1]},
        'ctd.6':{'comm':'Control data l7','cond':'','kw':['NOPTMAX', 'PHIREDSTP',
                    'NPHISTP', 'NPHINORED','RELPARSTP', 'NRELPAR'],'detail':
                 ['nb opti iteration', 'obj fct chg stop','nb iteration stop',
                'nb iteration stop 2','max param chg','nb iteration stop 3'],'type':['int'
                ,'float','int','int','float','int'],'default':[30, .005, 4, 4, .005, 4]},
        #SVD
        'svd.1':{'comm':'svd l1','cond':'','kw':['SVDMODE'],'detail':[['activate svd',
                'no','yes']],'type':['choice'],'default':[1]}, # EV 07/11
        'svd.2':{'comm':'svd l2','cond':'SVDMODE > 0','kw':['MAXSING', 'EIGTHRESH'], # EV 07/11
               'detail':['nb of singular value','eigenvalue ratio threshold'],'type':
                 ['string','float'],'default':['',5e-7]},
        'svd.3':{'comm':'svd l3','cond':'SVDMODE > 0','kw':['EIGWRITE'], # EV 07/11
               'detail':['svd output file'],'type':['int'],'default':[0]},
        #PGR
        'pgr.1':{'comm':'Nb of parameter group','cond':'','kw':['NBPARGP'],'detail':
                 ['nb of parameter grp (max=5)'],'type':['int'],'default':[1]},
        'pgr.2':{'comm':'Parameter group 1','cond':'NBPARGP > 0','kw':['PARGPN_1',
                'INCTYP_1','DERINC_1','DERINCLB_1','FORCEN_1','DERINCMUL_1',
                'DERMTHD_1'],'detail':['parameter grp name',['parameter increments method',
                'relative','absolute','rel_to_max'],'parameter increment','derivative increment',
                ['derivatives calculation','switch','always_2','always_3','switch_5',
                 'always_5'],'lower bound incretment',
                ['central derivatives','parabolic','outside_pts','best_fit','minvar',
                 'maxprec']],'type':['title','choice','float','float','choice','float', #EV 25/07/19
                'choice'],'default':['',0,0.01,0,0,2.0,0]},
        'pgr.3':{'comm':'Parameter group 2','cond':'NBPARGP > 1','kw':['PARGPN_2',
                'INCTYP_2','DERINC_2','DERINCLB_2','FORCEN_2','DERINCMUL_2',
                'DERMTHD_2'],'detail':['parameter grp name',['parameter increments method',
                'relative','absolute','rel_to_max'],'parameter increment','derivative increment',
                ['derivatives calculation','switch','always_2','always_3','switch_5',
                 'always_5'],'lower bound incretment',
                ['central derivatives','parabolic','outside_pts','best_fit','minvar',
                 'maxprec']],'type':['title','choice','float','float','choice','float',
                'choice'],'default':['',0,0.01,0,0,2.0,0]},
        'pgr.4':{'comm':'Parameter group 3','cond':'NBPARGP > 2','kw':['PARGPN_3',
                'INCTYP_3','DERINC_3','DERINCLB_3','FORCEN_3','DERINCMUL_3',
                'DERMTHD_3'],'detail':['parameter grp name',['parameter increments method',
                'relative','absolute','rel_to_max'],'parameter increment','derivative increment',
                ['derivatives calculation','switch','always_2','always_3','switch_5',
                 'always_5'],'lower bound incretment',
                ['central derivatives','parabolic','outside_pts','best_fit','minvar',
                 'maxprec']],'type':['title','choice','float','float','choice','float',
                'choice'],'default':['',0,0.01,0,0,2.0,0]},
        'pgr.5':{'comm':'Parameter group 4','cond':'NBPARGP > 3','kw':['PARGPN_4',
                'INCTYP_4','DERINC_4','DERINCLB_4','FORCEN_4','DERINCMUL_4',
                'DERMTHD_4'],'detail':['parameter grp name',['parameter increments method',
                'relative','absolute','rel_to_max'],'parameter increment','derivative increment',
                ['derivatives calculation','switch','always_2','always_3','switch_5',
                 'always_5'],'lower bound incretment',
                ['central derivatives','parabolic','outside_pts','best_fit','minvar',
                 'maxprec']],'type':['title','choice','float','float','choice','float',
                'choice'],'default':['',0,0.01,0,0,2.0,0]},
        'pgr.6':{'comm':'Parameter group 5','cond':'NBPARGP > 4','kw':['PARGPN_5',
                'INCTYP_5','DERINC_5','DERINCLB_5','FORCEN_5','DERINCMUL_5',
                'DERMTHD_5'],'detail':['parameter grp name',['parameter increments method',
                'relative','absolute','rel_to_max'],'parameter increment','derivative increment',
                ['derivatives calculation','switch','always_2','always_3','switch_5',
                 'always_5'],'lower bound incretment',
                ['central derivatives','parabolic','outside_pts','best_fit','minvar',
                 'maxprec']],'type':['title','choice','float','float','choice','float',
                'choice'],'default':['',0,0.01,0,0,2.0,0]},
        #OBS
        'obs.1':{'comm':'obs tranformation','cond':'','kw':['OB_TRANS'],
                 'detail':[['obs trans','no trans','log10','root square']],'type':
                 ['choice'],'default':[0]},
        #REG
        'reg.1':{'comm':'reg l1','cond':'PESTMODE ==2','kw':['PHIMLIM','PHIMACCEPT',
                'FRACPHIM','MEMSAVE'],'detail':['target obj fct','acceptable obj fct',
                'fraction of obj fct',['conservation of memory','-','memsave','nomemsave']],
                 'type':['float','float','float','choice'],'default':[1e-12,1.1e-12,0.1,0]},
        'reg.2':{'comm':'reg l2','cond':'PESTMODE==2','kw':['WFINIT','WFMIN','WFMAX',
                'LINREG','REGCONTINUE'],'detail':['reg weight factor',
                'min reg weight factor','max reg weight factor',['reg constraint','-',
                'linreg','nonlinreg'],['continue reg','-','continue','nocontinue']],
                 'type':['float','float','float','choice','choice'],
                 'default':[1.0,1.0e-10,1.0e+10,0,0]},
        'reg.3':{'comm':'reg l3','cond':'PESTMODE==2','kw':['WFFAC','WFTOL','IREGADJ',
                'NOPTREGADJ','REGWEIGHTRAT','REGSINGTHRESH'],'detail':['weight factor',
                'convergence criterion','inter regularization','reg iteration interval',
                'ratio regularization weight', 'singular value'],'type':['float','float',
                'int','string','string','string'],'default':[ 1.3,1.0e-2,1,'','','']},#'int','float','float'
        #PPP
        'ppp.1':{'comm':'singular value decomposition','cond':'','kw':['MAX_N_SUPER',
                'SUPER_EIGTHRES','N_ITER_BASE','N_ITER_SUPER','SVD_PACK',
                'MAX_SUPER_FRZ_ITER','MAT_INV','SUPER_RELPARMAX'],'detail':[],'type':
                 ['string','string','string','string','string','string','string','string'], #['int','float','int','int','title','int','title','float'],
                 'default':['','','','','','','','']},
        'ppp.2':{'comm':'other parameter','cond':'','kw':['AUTO_NORM',
                'LAMBDAS','MAX_REG_ITER','MAX_RUN_FAIL','ITERATION_SUMMARY','DER_FORGIVE']
                 ,'detail':[],'type':['string','string','string','string','string','string'], #['int','title','int','int','int','int'],
                 'default':['','','','','','']},
        'ppp.3':{'comm':'uncertainties','cond':'','kw':['UNCERTAINTY',
                'FORECASTS','PARAMETER_COVARIANCE'],'detail':[],'type':['string','string', #int title title
                'string'],'default':['','','']},
        }
