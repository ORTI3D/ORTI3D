# -*- coding: utf-8 -*-
"""
Created on Mon Dec 09 22:30:16 2013

@author: olive
This module is to deal with transient things
"""
from config import *
import numpy as np

def makeTransientTable(core):
    """creates a transient table where all zones of all lines are gathered
    and create a common period timescale
    zone keys 'number','name','coords','layer','value','type'
    each zone key contains a list of zones"""
    dZone = {'Transient':{},'zlist':{},'tlist':[]};
    tlist=[]
    # make the list of times from zones
    for mod in core.modelList:
        diczin = core.diczone[mod].dic
        llines = list(diczin.keys()) # gives the lines where they are zones
        for line in llines: 
            nz = len(diczin[line]['name']);#print 'timeper',line
            for iz in range(nz):
                slist=text2list(diczin[line]['value'][iz]);
                if type(slist)==type([5,6]): 
                    if len(slist)>1:
                        for s in slist: tlist.append(float(s[0]))
    # makes the list for the times in addin time
    tf,step = core.dicaddin['Time']['final'],core.dicaddin['Time']['steps']
    #tmode = core.dicaddin['Time']['mode'] # EV 18/02/19
    if type(tf) != type([4,5]): tf,step = [tf],[step]
    t0,wtimes = 0,[]
    ndec=max(0,int(-floor(log10(float(step[0])))))
    for i in range(len(tf)):
        t1,st = float(tf[i]),float(step[i])
        lll=around(arange(t0,t1,st),ndec)
        if ndec==0:lll=lll.astype('int')
        wtimes.extend(list(lll)) # EV 18/02/19
        t0 = float(tf[i])
    lll=around(float(tf[-1]),ndec)
    if ndec==0: lll=int(lll)
    wtimes.append(lll)
    dZone['wtimes']=wtimes[1:] # remove time 0
    # combines both lists
    tlist.extend(wtimes)
    tlist.sort();tlist = unique(tlist);#print 'timeper',len(tflow),tflow,tlist
    tlist = tlist[tlist<=float(tf[-1])] # to shorten if final time is smaller than times in zones
    #print ('tflow',tflow,'tlist',tlist) 
    dZone['tlist']=tlist
    core.nper = len(tlist)
    #create the list of transient zones and their values along time
    for mod in core.usedModelList:    
        diczin = core.diczone[mod].dic
        llines = list(diczin.keys()) # gives the lines where they are zones
        for line in llines: # loop into var to get transient values for each zone
            dZone['Transient'][line]=False # says if the line has transient zones or not
            dZone['zlist'][line]=[] # this is the list of zones for the considered line
            #dZone[line]=None
            lZ = diczin[line];
            nbzones = len(lZ['value'])
            if nbzones<=0 : continue
            a = ['              ']*len(tlist)
            tbl = reshape(array(a*nbzones),(len(tlist),nbzones))
            for iz in range(nbzones):
                # get the list, 
                slist = text2list(lZ['value'][iz]);#print 'tper',line,iz,slist
                if slist[0] in ['',' ']: continue # void zone
                elif type(slist)==type('str'): tbl[:,iz]=slist # a string for two or more values
                elif len(slist) == 1: 
                    tbl[:,iz]=slist[0] # put only the first value
                else : # a list of several values : transient
                    vprec=0;tl2=[];vl2=[]
                    dZone['Transient'][line]=True
                    dZone['zlist'][line].append(iz)
                    for s in slist:
                        tl2.append(float(s[0]));vl2.append(s[1])
                    for it in range(len(tlist)):
                        t = tlist[it]
                        if t in tl2:
                            v=vl2[tl2.index(t)];tbl[it,iz]=str(v);vprec=v
                        else : 
                            tbl[it,iz]=str(vprec)
            dZone[line]=tbl
    #print('dzonz',dZone)
    return dZone;#print 'Aq ztrans',dZone

def text2list(txt):
    '''from the text of a zone returns a list of transient vlaues'''
    if type(txt) in [type(1),type(1.)]: return[txt]
    if len(txt)==0: return [txt]
    parms=''
    if txt[0]=='$': # OA 25/4/19 all the lines in if loop modif(several parameters case)
        a,b,c = txt.split('$')
        parms = b.replace('\n','') # only the 3rd part contains the value
        txt = c
        if len(txt.split('\n'))<2 : return parms # no transient information
    a = txt.split('\n');#print'timp', a
    if len(a)==1: return [a[0]] # one value
    if len(a[1].replace(' ',''))==0: return [a[0]] # blanks in 2nd line : one value
    lout = []
    #print txt,lout
    for n in a:
        b1 = n.split()
        if len(b1)>1:
            if len(parms)>0: lout.append((b1[0],b1[1]+' '+parms.split()[1]))
            else : lout.append((b1[0],b1[1]))
    return lout
