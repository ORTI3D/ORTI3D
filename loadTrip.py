# -*- coding: utf-8 -*-
"""
Created on Fri Jun  5 07:10:26 2020

@author: olivi
"""
import os
fdir = 'D:\\iPht3d\\Exv2\\Anchor\\trip_mdf'
os.chdir(fdir)

# read discretization file
names = ['top'+str(i) for i in range(10)];names.append('botm')
line='';l0=[]
f1=open('Trip.dis','r');s=f1.read();f1.close();
s=s.split('     29')
disx,disy = [],[]
for a in s[1].split('\n')[1:]: disx.extend(a.split())
disx = array(disx).astype('float');savetxt('disx.txt',disx)
for a in s[2].split('\n')[1:]: disy.extend(a.split())
disy = array(disy).astype('float');savetxt('disy.txt',disy)
sdx,sdy = ' '.join('%.3e' %x for x in disx)+'\n',' '.join('%.3e' %x for x in disy)+'\n'
prev = 1000;ny,nx = 253,549;nc=nx*ny
for i,n in enumerate(names): #does not work for the last one
    l0 = []
    for a in s[i+3].split('\n')[1:]: l0.extend(a.split())
    arr = reshape(l0[:nc],(ny,nx)).astype('float')
    arr = minimum(prev-1,arr); prev=arr*1
    a1 = '\n'.join(' '.join('%.3e' %x for x in y) for y in arr)
    f1=open(n+'.var','w');f1.write('-1\n');f1.write(sdx);f1.write(sdy);f1.write(a1);f1.close()

# read hydraulic conductivity
f1=open('Trip.lpf','r');s=f1.read();f1.close();
s=s.split('     11')
names = []
for i in range(10): names.extend(['kx'+str(i),'ky'+str(i)])
for i,n in enumerate(names):
    l0 = []
    for a in s[i+1].split('\n')[1:]: l0.extend(a.split())
    arr = reshape(l0,(ny,nx)).astype('float')
    a1 = '\n'.join(' '.join('%.3e' %x for x in y) for y in arr)
    f1=open(n+'.var','w');f1.write('-1\n');f1.write(sdx);f1.write(sdy);f1.write(a1);f1.close()

# read RCH file, recharge is constant
f1=open('Trip.rch','r');s=f1.read();f1.close();
s=s.split('     18')
l0 = []
for a in s[1].split('\n')[1:-2]: l0.extend(a.split())
arr = reshape(l0,(253,549)).astype('float')
a1 = '\n'.join(' '.join('%.3e' %x for x in y) for y in arr)
f1=open('rech.var','w');f1.write('-1\n');f1.write(sdx);f1.write(sdy);f1.write(a1);f1.close()

# read and write BC and then heads
f1=open('Trip.bas','r');s=f1.read();f1.close();
s1=s.split(' IBOUND')
names=[]
for i in range(9): names.append('bc'+str(i))
for i,n in enumerate(names):
    l0 = []
    for a in s1[i+1].split('\n')[1:-1]: l0.extend(a.split())
    arr = reshape(l0,(253,549)).astype('int')
    a1 = '\n'.join(' '.join('%.3e' %x for x in y) for y in arr)
    f1=open(n+'.var','w');f1.write('-1\n');f1.write(sdx);f1.write(sdy);f1.write(a1);f1.close()
s1=s.split(' HEADS')
names=[]
for i in range(9): names.append('head'+str(i))
for i,n in enumerate(names):
    l0 = []
    for a in s1[i+1].split('\n')[1:-1]: l0.extend(a.split())
    arr = reshape(l0,(253,549)).astype('float')
    a1 = '\n'.join(' '.join('%.3e' %x for x in y) for y in arr)
    f1=open(n+'.var','w');f1.write('-1\n');f1.write(sdx);f1.write(sdy);f1.write(a1);f1.close()
    
# read DRN,RIV file and build zone file
disx,disy = loadtxt('disx.txt'),loadtxt('disy.txt')
xsu,ysu = r_[0,cumsum(disx)],r_[0,cumsum(disy)];ymx = ysu[-1]
xg,yg = (xsu[1:]+xsu[:-1])/2,(ysu[1:]+ysu[:-1])/2
nx,ny = len(xg),len(yg)

def readTransient(fname,nvar):
    f1=open('Trip.'+fname,'r');
    for i in range(4): a=f1.readline()
    nbpt = int(a.split()[0])
    l = f1.readline()
    pa = ['','','']
    iz,iy,ix = l[:10],int(l[10:20]),int(l[20:30])
    pa = [l[n*10+30:n*10+40] for n in range(nvar)]
    coo = [(xg[int(ix)-1],ymx-yg[int(iy)-1])]
    mdlist = [int(iz)]
    z0 = ['z'+fname+'0','$'+'\n'.join(pa)+'$ '+pa[0],int(iz),coo]
    zlist,i0,pa1 = [],1,['','','']
    for i in range(1,nbpt):
        l = f1.readline()
        iz1,iy1,ix1 = l[:10],int(l[10:20]),int(l[20:30])
        pa1 = [l[n*10+30:n*10+40] for n in range(nvar)]
        if len(set(pa).difference(pa1))>0 or abs(ix1-ix)>1 or abs(iy1-iy)>1:
            zlist.append(z0)
            coo = [(xg[int(ix1)-1],ymx-yg[int(iy1)-1])]
            mdlist = [int(iz)]
            z0 = ['z'+fname+str(i0),'$'+'\n'.join(pa1)+'$ '+pa[0],int(iz),coo]
            i0+=1
        else :
            coo = (xg[int(ix1)-1],ymx-yg[int(iy1)-1])
            if int(iz1) not in mdlist : 
                mdlist.append(int(iz1))
                z0[2] = ','.join([str(a) for a in sort(mdlist)])
            z0[3].append(coo)         
        iz,iy,ix,pa = iz1*1,iy1*1,ix1*1,pa1*1
    f1.close();
    s = ''
    for z in zlist:
        s += 'Name\t'+z[0]+'\tValue\t'+z[1].replace('\n','\\n')+'\tMedia\t'+str(z[2])+'\tcoord\t'
        a = str(z[3]).replace(',','\t').replace('(','').replace(')','')
        s += a[1:-1]+'\n'
    f1=open(fname+'zone.txt','w');f1.write(s);f1.close()

readTransient('drn',2)
readTransient('riv',3)
readTransient('ghb',2)
readTransient('wel',1)

# read btn file
f1=open('TRIPmt.btn','r');s=f1.read();f1.close();
s=s.split(' ICs')
# try to build a mask of different solutions
ic = zeros((8,8,253,549))
for j in range(8):
    for i in range(8):
        l0 = []
        for a in s[j*10+i+1].split('\n')[1:]: l0.extend(a.split())
        if len(l0)>1000: ic[j,i] = reshape(l0[:138897],(253,549)).astype('float')
l_esp=['F','C(4)','P','Al','Ca','Na','Si','Mg']
for j in range(8):
    for i in range(8):
        a1 = '\n'.join(' '.join('%.3e' %x for x in y) for y in ic[j,i])
        f1=open(l_esp[j]+'_l'+str(i)+'.var','w');f1.write(sdx);f1.write(sdy);f1.write(a1);f1.close()

######################## new read transient to matrix ########################

disx,disy = loadtxt('disx.txt'),loadtxt('disy.txt')
nr,nc = len(disy),len(disx)
sdx,sdy = ' '.join('%.3e' %x for x in disx)+'\n',' '.join('%.3e' %x for x in disy)+'\n'
m=loadtxt('Drntransient.txt');nvar=2;name='drn'
m=loadtxt('Rivtransient.txt');nvar=3;name='riv'
arr = zeros((nvar,nr,nc))
for i in range(1,14):  #layers as nb for mdflow
    m1=m[m[:,0]==i,:]
    if len(m1)==0: continue
    f1=open(name+str(i-1)+'.gvar','w');f1.write('-1 '+str(nvar)+'\n');
    f1.write(sdx);f1.write(sdy);
    for iv in range(nvar):
        arr[iv,m1[:,1].astype('int')-1,m1[:,2].astype('int')-1]=m1[:,3+iv]
        a1 = '\n'.join(' '.join('%.3e' %x for x in y) for y in arr[iv])
        f1.write(a1+'\n');
    f1.close()
        
