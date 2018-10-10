# -*- coding: utf-8 -*-
"""
Created on Thu Jun 29 18:23:47 2017

@author: oatteia
This class provides a way to calculate a pe-ph or stability diagram from a 
chemical composition given at an observation point.
The drawing will be done through a plot routine 
its midification can be called by an event (movement of a point), this part is
done in the addin module
"""

class Diagram:
    def __init__(self,core):
        self.dchem = core.addin.pht3d # the chemical database
        
    def getElts(self, dic_options):
        '''get the existing elements in the chemistry dic to expose them in
        a dialog for the user to choose, returns a list of elements
        dic_options['type'] in pepH or stability'''
        lmajor = ['Ca','Na','Mg','K']
        lmetal = ['Fe','Mn','As','Cu','Cr','Pb','Zn','U']
        lout = []
        if dic_options['type']=='stability':
            for elt in lmajor :
                idx = self.dchem['Solutions']['rows'].index(elt)
                if self.dchem['Solutions']['data'][idx][0]: lout.append(elt)
        if dic_options['type']=='pepH':
            for elt in lmetal :
                idx = self.dchem['Solutions']['rows'].index(elt)
                if self.dchem['Solutions']['data'][idx][0]: lout.append(elt)
        return lout
    
    def getReactionsFromElt(self,elt):
        '''get all possible phases and complexes that correspond to the elements
        and returns them as a dict with keys for each phase or complex
        and content of the key a list of two components : the reaction
        and the logk
        the reactions are already stored in tempDbase in pht3d
        '''
        dicDB = self.dchem.tempDbase
        master = dicDB['SOLUTION_MASTER_SPECIES']
        lelt = []
        # search for redox states or major ion
        for k in master.keys():
            k1 = k.split('(')[0]
            if elt == k1: lelt.append(master[k][0])
        self.lelt = lelt
        self.dcplx = self.splitReactions(dicDB['SOLUTION_SPECIES'],'complex',lelt)
        self.dminx = self.splitReactions(dicDB['PHASES'],'mineral',lelt)
        self.onlyMasterSpecies()
        
    def splitReactions(self,dicR,typR,lelt):
        '''from a dictionnary containing reactions written in pht3D style
        retrives the species and stoechio
        in case of complex the complex is the first after = sign
        in case of phases the phase is the first species'''
        dicOut = {}
        for k in dicR.keys():
            for e in lelt:
                if e in dicR[k][0]: 
                    reac,logk = dicR[k]
                    r1,r2 = reac.split('=')
                    dct,sgn = {},-1.0
                    if typR == 'complex': 
                        name,dct['Cplx'] = r2.split()[0],1 ##print cname
                    elif typR == 'mineral': 
                        name,dct['Cplx'] = r1.split()[0],0
                    dct = self.splitSide(dct,r1,sgn)
                    dct = self.splitSide(dct,r2,-sgn)
                    dct.pop(name)
                    if self.allEltsInMaster(dct) == False : continue
                    dct['logk'] = float(logk)
                    dct['flag'] = True
                    dicOut[name] = dct
        return dicOut
        
    def splitSide(self,dct,r1,sgn):
        r1a = r1.split(' + ')
        for sp in r1a: 
            if len(sp.split())>1 : 
                nb, sp1 = sp.split()
                dct[sp1.strip()] = sgn*float(nb)
            else :
                dct[sp.strip()] = sgn
        return dct
        
    def allEltsInMaster(self,dct):
        dbsolu = self.dchem.Base['Chemistry']['Solutions']
        dmaster = self.dchem.tempDbase['SOLUTION_MASTER_SPECIES']
        nsp, count = len(dct.keys()),0
        for sp in dct.keys():
            if sp in ['logk','flag','H+','OH-','e-','H2O','Cplx']: 
                count+=1;continue
            for km in dmaster.keys():
                if len(dmaster[km])>1 : 
                    if sp == dmaster[km][0]:
                        indx = dbsolu['rows'].index(km)
                        if dbsolu['data'][indx][0] == 1: # we got it!
                            count +=1
        return count==nsp
        
    def onlyMasterSpecies(self):
        '''this routines allow to search for the reactions that do not
        conatin only master species and transform them to on which contains
        only master species'''
            
        def replaceSpecies(db_cplx,reac,spec):
            # replace the bad species by a master one using the logk
            if spec == 'Cplx': return reac
            r2, logk2 = db_cplx[spec]
            d1 = self.splitSide({},r2.split('=')[0],1.0)
            for k in d1.keys():
                if k not in reac.keys(): reac[k] = d1[k]*1
                else : reac[k] += d1[k] # adding the stoechio
            reac.pop(spec)
            reac['logk'] += float(logk2)
            return reac
        # create the list of solutions master species in lmaster
        d0 = self.dchem.tempDbase['SOLUTION_MASTER_SPECIES']
        lmaster = []
        for k in d0.keys(): 
            if type(d0[k])==type(range(3)): lmaster.append(d0[k][0])
        lmaster.extend(['flag','logk','H2O','e-'])
        # identify if all species are in master
        db_cplx = self.dchem.tempDbase['SOLUTION_SPECIES']
        for k in self.dminx.keys():
            for k1 in self.dminx[k].keys():
                if k1 not in lmaster:
                    self.dminx[k] = replaceSpecies(db_cplx,self.dminx[k],k1)
        for k in di.dcplx.keys():
            for k1 in di.dcplx[k].keys():
                if k1 not in lmaster:
                    self.dcplx[k] = replaceSpecies(db_cplx,self.dcplx[k],k1)
        
    def makeBaseArray(self):
        ''' form a dict of the seledted elements creates a talbe of all possible 
        reactions '''
        def addLine(lsp,dsp,lelt,data,nb,istart):
            # put the elements coeff at the right place for one line
            for i,sp in enumerate(lsp): 
                if dsp[sp]['flag']:
                    data[i+istart,0] = nb
                    for j,el in enumerate(lelt):
                        if dsp[sp].has_key(el):
                            data[i+istart,j+1] = dsp[sp][el]
                        else : 
                            data[i+istart,j+1] = 0
                    #data[i+istart,-1] = dsp[sp]['logk']
        # first make the list of elements in all reactions
        lelt = []
        for k in self.dcplx.keys(): 
            if self.dcplx[k]['flag']:lelt.extend(self.dcplx[k].keys())
        for k in self.dminx.keys(): 
            if self.dminx[k]['flag']:lelt.extend(self.dminx[k].keys())
        lelt = list(unique(lelt))
        lelt.remove('flag')
        nelts = len(lelt)
        # then make a table of reactions coeffs
        lcplx,lminx = self.dcplx.keys(),self.dminx.keys()
        ncplx = len(lcplx)
        lspecies = lcplx*1
        lspecies.extend(lminx)
        nsp =len(lspecies)
        data = zeros((nsp,nelts+1))
        addLine(lcplx,self.dcplx,lelt,data,1,0)
        addLine(lminx,self.dminx,lelt,data,2,ncplx)
        lelt.insert(0,'type')
        self.tblSpec = {'rows':lspecies,'cols':lelt,'data':data}
        
    def makeReacArray(self):   
        '''from the liine array, create all possible reactions, keeping track
        of the type of the reactants : at the beginning keep the previous reactions
        with the base species (type 1 and 2)
        11: 2 complex boundary, 12: complex/mineral
        boundary, 22: mineral/mineral boundary
        in phreeqc reacions cplx is ont he right side and mineral on the left one
        so for mineral/complex just add them and for smae type substract
        '''
        data = self.tblSpec['data']
        nsp,nelt = shape(data)
        self.tblSpec['cols'].insert(1,'sp2') # insert species number
        self.tblSpec['cols'].insert(1,'sp1')
        ipH,ipe = self.tblSpec['cols'].index('H+'),self.tblSpec['cols'].index('e-')        
        # create the borders of the diagram (for pe-pH onle now!!!)
        treac = zeros((nsp*(nsp+1)/2+4,nelt+2))
        apH,ape,alk = zeros(nelt+2),zeros(nelt+2),zeros(nelt+2)
        apH[ipH],ape[ipe],alk[-1] = 1,1,1
        treac[0,:] = apH  # left vertical side H+=10^0  
        treac[1,:] = ape + alk*14 # top side  e- = 10^-14
        treac[2,:] = apH -14*alk # right side
        treac[3,:] = ape -14*alk # bottom side 
        # add the reactions with base species
        treac[4:nsp+4,0] = data[:,0] #type 1 cplx, 2 mineral with base species
        treac[4:nsp+4,1] = range(nsp) # specie snumber
        treac[4:nsp+4,3:] = data[:,1:] # stoechio and logk
        i = nsp+4
        # find the major redox reaction (ex: Fe+2 = Fe+3 e-)
        for k in self.dcplx.keys():
            count = 0
            for k1 in self.dcplx[k].keys():
                if k1 in ['flag','logk','Cplx']: continue
                if self.dcplx[k] in self.lelt: continue
                count += abs(self.dcplx[k][k1])
            if count == 0: majreac = self.dcplx[k]
        for ie in range(nsp):
            for je in range(ie+1,nsp):
                ty1,ty2 = data[ie,0],data[je,0]
                typ = ty1*10+ty2 # 11 or 12 or 22
                treac[i,:3] = typ,ie,je
                ### we ned to keep only one speceis in Fe+2/Fe+3 in the reaction
                if ty1 == ty2: # same type of reactants
                    treac[i,3:] = data[ie,1:]-data[je,1:]
                else : # bdy btw minerla and complex
                    treac[i,3:] = data[ie,1:]+data[je,1:]
                i += 1
        self.tblReac = treac
        
    def calcDiagram(self):
        '''calculates the diagram from an array of lines, by starting on one
        corner and finding the closes line that concerns the species 
        (elt, mineral or complex)
        '''
        pass
        
    def recalculate(self):
        pass