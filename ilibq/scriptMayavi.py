from mayavi.mlab import *
from mayavi.modules.outline import Outline
from mayavi.modules.surface import Surface

it=38;iesp=0
#data = core.gui.guiShow.arr3
data = core.transReader.readUCN(core,'Mt3dms',it,0,'Tracer')
data= data[:,-1::-1,:]
xm,ym,zm = getMesh3Dcenters(core)
xm = (xm-amin(xm))/(amax(xm)-amin(xm))
ym = (ym-amin(ym))/(amax(ym)-amin(ym))
zm = (zm-amin(zm))/(amax(zm)-amin(zm))
obj = contour3d(xm,ym,zm,data, contours=4,transparent=True)
#xm,ym,zm = getMesh3Dcenters(core)
#u,v,w= core.gui.guiShow.get3Dvectors(it)
#print shape(xm),shape(u),shape(v),shape(w)
#obj=quiver3d(xm,ym,(zm-695)*10,u,v,w)