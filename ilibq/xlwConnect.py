# -*- coding: utf-8 -*-
"""
Created on Sat Jul 29 14:28:11 2017
to connect an orti3d model to excel
@author: olive
"""
import os,sys
import xlwings as xlw
os.chdir('E:\\iPht3d\\LibDev3\\ilibq')
import core
md = core.Core()

os.chdir('E:\ipht3d\Exv2\Clement\Kmousse2D')
wb = xlw.Book('xlConnect.xlsx')

shname = 'Mod1' # Mod2
def printParameters(md,wb,shname):
    sh = wb.sheets[shname]
    fname = sh.range('B1').value;print(fname)
    md.openModel(os.getcwd(),fname)
    # write the values of the params
    n = 3 # first line to write things
    for a in range(2):
        mname,line = sh.range('A'+str(n)).value, sh.range('B'+str(n)).value
        dcz = md.diczone[mname].dic[line]
        sh.range('D'+str(n)).value = dcz['name']
        sh.range('D'+str(n+1)).value = dcz['value']
        xc = [mean(list(zip(*xy))[0]) for xy in dcz['coords']]
        sh.range('D'+str(n+2)).value = xc
        yc = [mean(list(zip(*xy))[1]) for xy in dcz['coords']]
        sh.range('D'+str(n+3)).value = yc
        dmoy = []
        for i in range(len(xc)):
            x,y = list(zip(*dcz['coords'][i]))
            dmoy.append(mean(sqrt((x-xc[i])**2+(y-yc[i])**2)))
        sh.range('D'+str(n+4)).value = dmoy
        n += 5       
        
#printParameters(md,wb,shname)
def center2coo(dcz,xcnew,ycnew,dnew):
    coonew = []
    for i in range(len(xcnew)):
        x,y = list(zip(*dcz['coords'][i]))
        xc,yc = mean(x),mean(y)
        d = mean(sqrt((x-xc)**2+(y-yc)**2));#print(dnew[i]/d)
        xn,yn = xcnew[i]+(x-xc)*dnew[i]/d,ycnew[i]+(y-yc)*dnew[i]/d
        coonew.append(list(zip(xn,yn)))
    return coonew

lfor=['s71-38','s71-126','s71-127','s71-128','s71-130','s71-131','s71-132','s71-133']
pos=14
def runModel(md,wb,shname):
    sh = wb.sheets[shname]
    fname = sh.range('B1').value;print(fname)
    md.openModel(os.getcwd(),fname)
    for nl in [3,8]:
        mname,line = sh.range('A'+str(nl)).value, sh.range('B'+str(nl)).value
        dcz = md.diczone[mname].dic[line]
        l0 = sh.range('D'+str(nl+1)).expand('right').value
        dcz['value'] = [str(x) for x in l0]
        lx = sh.range('D'+str(nl+2)).expand('right').value
        ly = sh.range('D'+str(nl+3)).expand('right').value
        ld = sh.range('D'+str(nl+4)).expand('right').value
        dcz['coords'] = center2coo(dcz,lx,ly,ld)
    print(md.diczone['Modflow'].dic['lpf.8']['value'])
    #print(md.diczone['Modflow'].dic['lpf.8']['coords'][0])
    md.baseDir='E:\\iPht3d\\LibDev3'
    md.writeModel('Modflow')
    md.runModel('Modflow')
        
def getResults(md,wb,shname,znames,pos):
    sh = wb.sheets[shname]
    zname=znames[0]
    t,bkt,labl = md.onPtObs('B0',1,'Flow',zname,['Head'])
    sh.range('B'+str(pos)).value = t
    pos += 1
    sh.range('A'+str(pos)).value = zname
    sh.range('B'+str(pos)).value = ravel(bkt)
    for zname in znames[1:]:
        pos+=1
        sh.range('A'+str(pos)).value = zname;
        t,bkt,labl = md.onPtObs('B0',1,'Flow',zname,['Head']);
        sh.range('B'+str(pos)).value = ravel(bkt);
        
runModel(md,wb,shname)
getResults(md,wb,shname,lfor,pos)  