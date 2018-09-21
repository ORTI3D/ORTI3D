import os
import core as corebase
from geometry import *
from importExport import *

core = corebase.Core()
##d='d://ipht3d//iqpht3d//ex//test0'
##core.openModel(d,'test')
##c1='25 50 \n'
##z1={'number':[0],'name':['t1'],'coords':[c1],'layer':[1],'value':[12.],
##    'type':'B.Condition'}
##core.diczones['Modflow']['1.6']=z1
##grd=zone2grid(core,'Modflow','1.6',z1)

imp=impFile(core,'Modflow')
fdir='d://ipht3d//ex//neuchatel//neuchatel1'
fname='acidplume.nam'
m=imp.impAsciiModflow(fdir,fname)
