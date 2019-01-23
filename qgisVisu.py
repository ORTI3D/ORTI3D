# -*- coding: utf-8 -*-
"""
Created on Thu Feb 20 23:43:54 2014

@author: olive
This is the visualization tol box for qgis that creates the link between orti3d
and Qgis
"""
from .geometry import *
from .config import *
from .qtDialogs import *

from PyQt5.QtCore import *
from PyQt5.QtGui import *
from qgis.core import *
from osgeo import ogr

class qgisVisu:
    def __init__(self,iface,gui,core):
        self.iface, self.gui,self.core = iface,gui,core
        self.canvas = self.iface.mapCanvas()
        self.linelist = []
        self.cfg = Config(core)
        self.dialogs = self.cfg.dialogs
        self.showzdlg = 1
                      
    # for compatibility purpose
    def addZone(self,a,b,c,d): pass
    def redraw(self,line): # needed when one object was created and cancelled      
        pass #self.iface.mapCanvas().refresh() 
    def drawGrid(self,bool) : pass
    def drawImage(self,bool) : pass
    def drawMap(self,bool) : pass
    def drawVector(self,bool) : pass
    def drawContour(self,bool):pass
    def drawParticles(self,bool) : pass
    def changeAxesOri(self,plane): pass
    def showVar(self,line,imedia):
        '''to redraw a layer in qgis with zones only for the current media'''
        allLayers = self.canvas.layers()
        for layer in allLayers:
            #self.dialogs.onMessage(self.gui,layer.name()[1:].split('_')[0])
            if layer.name()[1:].split('_')[0] != line: continue
            #request = QgsFeatureRequest().setFilterExpression( u'"Media" = '+str(imedia) )
            strsel="gid IN ("
            for feat in features:
                if feat.geometry().intersects(geomPoly):
                    strsel=strsel+ "'" + str(feat.id()) + "',"
            #Closes the string with a parenthesis
            strsel=strsel[:-1] + ")"
            layer.setSubsetString(request)

    def retranslateUi(self, Var): pass
    def clearLayers(self):
        """clear all layers, used when opening a new model"""
        allLayers = self.canvas.layers()
        for layer in allLayers:
            QgsProject.instance().removeMapLayer(layer) # OA 17/10/18 map layer replaced by project
        
    def initDomain(self): 
        """initalize the grid and all cells in Qgis"""
        for modName in self.core.modelList:
            self.linelist.append(list(self.core.dickword[modName].lines.keys()))
        self.grd = self.core.addin.getFullGrid()
        grp = self.core.addin.getModelGroup()
        self.mesh=None
        if grp[:2]=='Op' and self.core.dicval[grp+'Flow']['domn.1'][0]>0: # opegeo
            self.mesh = self.core.addin.opgeo
        if grp[:2]=='Mi' and self.core.dicval[grp+'Flow']['spat.4'][0]>0: # min3p
            self.mesh = self.core.addin.min3p
        if self.core.mfUnstruct: 
            self.mesh = self.core.addin.mfU
        for n in ['Grid','Parameters']:
            self.createGridLayer(n)
        self.iface.digitizeToolBar().actions()[5].disconnect() # OA added 20/10/18
        self.iface.digitizeToolBar().actions()[5].triggered.connect(self.attributeForm) # OA added 20/10/18
        #if self.mesh != None: self.createNodeResultLayer(nodes)
            
    def createGridLayer(self,lname): #OA moved from initdomain to here 20/1/19 
        wkt = self.encodeGrid()
        layer = QgsVectorLayer('Polygon', lname, "memory")
        QgsProject.instance().addMapLayer(layer)
        layer.startEditing()
        layer.dataProvider().addAttributes( [QgsField("id", QVariant.Int) ] )
        layer.dataProvider().addAttributes( [QgsField("value", QVariant.Double) ] )
        layer.commitChanges()
        layer.startEditing()
        for i,w in enumerate(wkt.split('\n')[1:-1]):
            poly = QgsGeometry.fromWkt(w.split(';')[1])
            feat = QgsFeature()
            feat.setGeometry(poly)
            feat.setAttributes([i,0.])
            layer.addFeature(feat)
        layer.updateExtents()
        layer.commitChanges()
        symbols = layer.renderer()#EV 15/01/19
        mySymbol1 = QgsFillSymbol.createSimple({'color':'255,0,0,0','color_border':'#969696'})
        symbols.setSymbol(mySymbol1)
        layer.triggerRepaint()
        self.iface.layerTreeView().refreshLayerSymbology(layer.id())
        
    def createAndShowObject(self,dataM,dataV,opt,value,color):
        """create the object and show it, typObj in Grid,Contour, Vector, Map
        """
        group,name,value = self.gui.guiShow.getCurrentContour() # OA 21/1/19
        #if type(value) != type(bool) : name = value # OA 21/1/19 to find the name of the variable
        #self.dialogs.onMessage(self.gui,'layer name : '+str(name))
        if dataM == None : self.drawContour(False)
        else : 
            if self.gui.guiShow.swiImg == 'Contour':
                self.createContour(dataM,value,color,name) #‼ added name for Qgis layer
            else : 
                self.createImage(dataM)
        if dataV == None : self.drawVector(False)
        else : self.createVector(dataV)
                        
    def drawObject(self,typObj,bool):
        exec('self.draw'+typObj+'('+str(bool)+')')
        
    def createLayer(self,modName,line,typList):
        """create one new layer from a line setting attributes to zones
        called when opening a model through core2qgs
        and called alone when creating a new zone layer
        typList contain types (Linstring, or Points or both)"""
        self.linelist.append(line)
        dick = self.core.dickword[modName]
        comm = dick.lines[line]['comm'][:25]
        # create the layer in qgis
        #self.dialogs.onMessage(self.gui,'layer'+line)
        layerList = []
        for tp in typList:
            layerName = tp[0]+line + '_' + comm
            if len(self.canvas.layers())>0:
                crs = self.canvas.layers()[0].crs()# finds the coordinates system
                layer = QgsVectorLayer(tp+"?crs=" + crs.authid(), layerName, "memory")#, "delimitedtext")
            else:
                layer = QgsVectorLayer(tp, layerName, "memory")
            QgsProject.instance().addMapLayer(layer) # OA modified  2/10
            layer.startEditing()
            layer.dataProvider().addAttributes( [QgsField("id", QVariant.Int) ] )
            layer.dataProvider().addAttributes( [QgsField("name", QVariant.String) ] )
            layer.dataProvider().addAttributes( [QgsField("value", QVariant.Double) ] )
            layer.dataProvider().addAttributes( [QgsField("media", QVariant.String) ] )
            layer.dataProvider().addAttributes( [QgsField("type", QVariant.Map) ] )
            layer.commitChanges()
            layerList.append(layer)
            self.layer = layer
            layer.editingStarted.connect(self.editing_started)
            #layer.committedFeaturesAdded.connect(self.cfeature_added)
            #layer.layerModified.connect(self.lay_mod)
        return layerList
        
    def getLayer(self,modName,line):
        dick = self.core.dickword[modName]
        comm = dick.lines[line]['comm'][:25]
        layerList = []
        for n in ['point','line']:
            layerName = line + '_'+n+' ' + comm
            for layer in allLayers:
                if layer.name() == layerName : layerList.append(layer)
        return layerList
        
    def editing_started(self):
        #self.dialogs.onMessage(self.gui,'edit start')
        QSettings().setValue(
                '/qgis/digitizing/disable_enter_attribute_values_dialog', True)
        layer = QgsProject.instance().sender() # OA 17/10/18 map layer replaced by project
        layer.featureAdded.connect(self.feature_added)
        layer.beforeCommitChanges.connect(self.editing_stopped)
        #self.showzdlg = 1
        
    def editing_stopped(self):
        #self.dialogs.onMessage(self.gui,'before com change')
        layer = QgsProject.instance().sender() # OA 17/10/18 map layer replaced by project
        QSettings().setValue(
                '/qgis/digitizing/disable_enter_attribute_values_dialog', False)
        try : layer.featureAdded.disconnect()
        except TypeError: a=0
        #self.showzdlg = 0

    def feature_added(self): # OA completely modified 20/10/18
        """action when the feature editing is finished in Qgis, show the zone dialog
        and then set the values in Qgis"""
        layer = self.iface.activeLayer() # OA 20/10/18 map
        line = layer.name()[1:].split('_')[0]
        if layer.name()[0] =='L': typP ='asPolyline()'
        else: typP='asPoint()'
        feats = layer.getFeatures() # list of existing features
        for i,feat in enumerate(feats):
            if i==0 : break  # the new one is the first in the list
        coords = eval('feat.geometry().'+typP) # OA 17/10/18 exec to eval
        if type(coords) != type([5]): coords = [coords] # for points
        dicz = self.prepare(line)
        self.gui.varBox.base.onZoneCreate(None, coords)
        if line in dicz.dic:
            nf = len(dicz.dic[line]['name'])-1
        else : 
            nf = 0
        self.setAttributesInQgis(layer,dicz,nf,feat)
            
    def attributeForm(self): # OA added 20/10/18
        '''replace the attributefrom of Qgis by our zone dialog'''
        layer = self.iface.activeLayer() # OA 20/10/18 map
        line = layer.name()[1:].split('_')[0]
        dicz = self.prepare(line)
        feat = layer.selectedFeatures()[0]
        fname = feat.attributes()[0] # this is the name
        zlist = dicz.dic[line]
        nf = zlist['name'].index(fname)
        # make the dialog
        dialg = zoneDialog(self,self.core,self.gui.currentModel,line,zlist,nf)
        result = dialg.saveCurrent()
        self.core.makeTtable()
        self.setAttributesInQgis(layer,dicz,nf,feat)
        
    def prepare(self,line): # OA added 20/10/18
        categ = self.gui.currentCategory
        self.gui.currentLine = line
        self.gui.currentModel = self.gui.varBox.base.modelFromLine(categ,line)
        dicz = self.core.diczone[self.gui.currentModel]
        return dicz
    
    def setAttributesInQgis(self,layer,dicz,nf,feat): # OA added 20/10/18
        #set attributes in QGis, nf is the feature/zone number
        line = layer.name()[1:].split('_')[0]
        name = str(dicz.getValue(line,'name',nf)) # has been added in diczone in the last position
        try: value = float(dicz.getValue(line,'value',nf))
        except ValueError: value=0
        media = dicz.getValue(line,'media',nf)
        layer.beginEditCommand("modify")
        feat.setAttributes([name,value,media,0]) # OA 2/10 nf removed
        coords = dicz.getValue(line,'coords',nf)
        lcoord = [QgsPoint(float(c[0]),float(c[1])) for c in coords]
        feat.setGeometry(QgsGeometry.fromPolyline(lcoord));
        layer.updateFeature(feat)
        layer.commitChanges()   
        layer.endEditCommand()
            
    def createNodeResultLayer(self,nodes):
        """specific to MIn3p when the results are given at nodes and not for elements"""
        allLayers = self.canvas.layers();#print allLayers
        layerName = 'NodesResults'
        # search if the layer already exists, if yes, returns it
        for layer in allLayers:
            if layer.name() == layerName : return [layer]
        # creating the layer
        if len(allLayers)>0:
            crs = allLayers[0].crs()# finds the coordinates system
            layer = QgsVectorLayer("Point?crs=" + crs.authid(), layerName, "memory") # "delimitedtext")
        else:
            layer = QgsVectorLayer("Point", layerName, "memory") # "delimitedtext")
        QgsProject.instance().addMapLayer(layer) # OA 17/10/18 map layer replaced by project
        layer.startEditing()
        provider = layer.dataProvider()
        provider.addAttributes( [QgsField("id", QVariant.Int) ] )
        provider.addAttributes( [QgsField("value", QVariant.Double) ] )
        nnodes,nc = shape(nodes)
        for ir in range(nnodes):
            feat = QgsFeature();
            feat.setGeometry(QgsGeometry.fromPointXY(QgsPointXY(nodes[ir,0],nodes[ir,1])));
            feat.setAttributes([ir,0.])
            provider.addFeatures([feat]);
        layer.commitChanges()
        
        
    def zonesQgs2core(self):
        '''creates zones from existing polys in Qgis layers (done when saving qgis)'''
        allLayers = self.canvas.layers()
        for layer in allLayers:
            if layer.name()[:4] == 'Grid' : continue
            line,typP = layer.name().split()[0],'asPolyline()'
            if len(line.split('_'))>1: 
                line,typ0 = line.split('_',1) #EV 15/01/19
                if typ0=='point': typP='asPoint()'
            if line not in self.linelist : continue
            for modName in self.core.modelList:
                if line in list(self.core.dickword[modName].lines.keys()):
                    dicz = self.core.diczone[modName]
                    break
            feats = layer.getFeatures()
            dicz.dic[line]={'number':[],'name':[],'coords':[],'media':[],'value':[],'type':[]}                
            #QgsMessageLog.logMessage('linee '+str(line), 'MyPlugin')
            for i,f in enumerate(feats):
                dicz.addZone(line)
                dicz.setValue(line,'value',i,f['value'])
                dicz.setValue(line,'name',i,f['name'])
                dicz.setValue(line,'media',i,f['media'])
                coords = eval('f.geometry().'+typP) # OA 6/11/18
                if type(coords) != type([5]): coords = [coords] # for points
                dicz.setValue(line,'coords',i,coords)
                
    def removeOrtiLayers(self): # OA added 16/12/18 to remove layers when opening
        for modName in self.core.modelList:
            dick = self.core.dickword[modName]
            dicz = self.core.diczone[modName]
            lines = list(dicz.dic.keys())
            for line in lines:
                comm = dick.lines[line]['comm'][:25]
                nz = dicz.getNbZones(line)
                # get layer or layers
                typList = self.findTypList(modName,line)
                for tp in typList:
                    layerName = tp[0]+line + '_' + comm
                    ln0 = QgsProject.instance().mapLayersByName(layerName)
                    if len(ln0)>0: QgsProject.instance().removeMapLayer(ln0[0]) 
        for n in ['Grid','Parameters']: #EV added 15/01/2019
            ln0 = QgsProject.instance().mapLayersByName(n)
            if len(ln0)>0: QgsProject.instance().removeMapLayer(ln0[0]) 
        
    def zonesCore2qgs(self):
        """uses zones stored in core to fill the qgis layers with features
        done during model opening, uses createLayer and fill the attributes here"""
        for modName in self.core.modelList:
            dicz = self.core.diczone[modName]
            lines = list(dicz.dic.keys())
            for line in lines:
                nz = dicz.getNbZones(line)
                # get layer or layers
                typList = self.findTypList(modName,line)
                layerList = self.createLayer(modName,line,typList) # typlist is the list of type of layer
                providerList = [layer.dataProvider() for layer in layerList]
                # add zones in the layer(s)
                for i in range(nz):
                    coords = dicz.getValue(line,'coords',i)
                    name = str(dicz.getValue(line,'name',i))
                    try: value = float(dicz.getValue(line,'value',i))
                    except ValueError: value=0
                    media = dicz.getValue(line,'media',i)
                    feat = QgsFeature()
                    if len(coords)==1:
                        pt = coords[0]
                        feat.setGeometry(QgsGeometry.fromPointXY(QgsPointXY(float(pt[0]),float(pt[1]))));
                        provider = providerList[typList.index('Point')]
                    else :
                        lcoords = []
                        for pt in coords:
                            lcoords.append(QgsPoint(float(pt[0]),float(pt[1])))
                        feat.setGeometry(QgsGeometry.fromPolyline(lcoords));
                        provider = providerList[typList.index('LineString')]
                    feat.setAttributes([name,value,media,0])
                    provider.addFeatures([feat])
                    
    def findTypList(self,modName,line):
        '''returns a list of object types as found in coords, possible:
        ['Points'],['Linestring'],['Points','Linestring']
        '''
        dicz = self.core.diczone[modName]
        ltyp = []
        if line not in list(dicz.dic.keys()): return ['Points','Linestring'] # OA 17/10/18
        for coo in dicz.dic[line]['coords']:
            if len(coo)==1 and 'Point' not in ltyp: ltyp.append('Point')
            if len(coo)>1 and 'LineString' not in ltyp: ltyp.append('LineString')
            if len(ltyp)==2: break
        return ltyp
            
    def createContour(self,data,value,color,layName='Results'):
        """creates a coutour or image object using the Grid layer 
        of Qgis"""
        layers = QgsProject.instance().mapLayers() # OA 22/1/19 start
        layers1  = [x.split('_')[0] for x in layers]
        #QgsMessageLog.logMessage(layName+' '+str(layers1), 'MyPlugin')
        if layName not in layers1: self.createGridLayer(layName) #
        layers = QgsProject.instance().mapLayers()
        layers1 = list(layers.keys()) # 
        layers2  = [x.split('_')[0] for x in layers]
        layer = layers[layers1[layers2.index(layName)]] # OA 22/1/19 end
        pr = layer.dataProvider();
        layer.startEditing()
        colorRamp = QgsGradientColorRamp.create({'color1':'255,0,0,255', 'color2':'0,0,255,255','stops':'0.25;255,255,0,255:0.50;0,255,0,255:0.75;0,255,255,255'})
        X,Y,Z = data
        d0,i = {},1
        # if self.mesh != None: # view results on node
        #     for i in range(len(Z)):
        #         d0[i]={0:i,1:float(Z[i])}                
        #     for layer in allLayers:
        #         if layer.name()=='NodeResults': break
        #     symbol = QgsMarkerSymbol()
        #     mode = QgsGraduatedSymbolRenderer.Jenks
        #     renderer = QgsGraduatedSymbolRenderer.createRenderer( layer,'value', 25, mode, symbol, colorRamp )
        #     renderer.setMode( QgsGraduatedSymbolRenderer.Custom )
        # else :
        ny,nx = shape(Z)
        for ix in range(nx):
            for iy in range(ny):
                d0[i]={0:i,1:float(Z[iy,ix])}
                i+=1
        symbol = QgsFillSymbol()
        mode = QgsGraduatedSymbolRenderer.Jenks
        renderer = QgsGraduatedSymbolRenderer.createRenderer(layer,'value',25,mode,symbol, colorRamp)
        renderer.setMode( QgsGraduatedSymbolRenderer.Custom )
        pr.changeAttributeValues(d0)
        layer.setRenderer(renderer)        
        layer.commitChanges()
        layer.triggerRepaint()
        self.iface.mapCanvas().refresh() 
        
    def encodeWKT(self,modName,line):
        """a format that is faster to read for the grid polygons"""
        dicz = self.core.diczone[modName].dic[line]
        ip, il,idx_pt,idx_li = 0, 0,[],[]
        sp,sl ='id;WKT\n','id;WKT\n'
        lcoords = dicz['coords']
        for i, coo in enumerate(lcoords):
            strcoo = ','.join([str(a).replace('(','').replace(')','').replace(',','') for a in coo])
            if len(coo)==1: 
                ip = 1
                idx_pt.append(i)
                sp+= str(i+1)+';POINT('+ strcoo + ')\n'
            else :
                il = 1
                idx_li.append(i)
                sl+= str(i+1)+';LINESTRING('+ strcoo + ')\n'
        if ip==0 : sp = None
        if il==0 : sl = None
        return [sp,sl,idx_pt,idx_li]
        
    def encodeGrid(self):
        s ='id;WKT\n'
        g = self.grd
        if self.mesh != None :
            xcoo = self.mesh.elx
            ycoo = self.mesh.ely
            nnodes = len(xcoo)
            for ir in range(nnodes):
                coo = list(zip(xcoo[ir],ycoo[ir]))
                strcoo = ','.join([str(a).replace('(','').replace(')','').replace(',','') for a in coo])            
                s += str(ir+1)+';POLYGON(('+ strcoo+ '))\n'
        else :
            count,a = 0,g['x0']
            for i in range(g['nx']):
                b,c = g['y0'],g['dx'][i]
                for j in range(g['ny']):
                    d =  g['dy'][j]
                    #self.dialogs.onMessage(self.gui,str(a)+' '+str(b)+' '+str(c)+' '+str(d))
                    coo = [(a,b),(a+c,b),(a+c,b+d),(a,b+d),(a,b)]
                    strcoo = ','.join([str(x).replace('(','').replace(')','').replace(',','') for x in coo])            
                    s += str(count+1)+';POLYGON(('+strcoo+'))\n'
                    b += d
                    count += 1
                a += c
        return s

'''        
    def encodeKml(self,modName,line):
        dicz = self.core.diczone[modName].dic[line]
        ip, il = 0, 0
        sp ='<?xml version="1.0" encoding="utf-8" ?>'
        sp+= '<kml xmlns="http://www.opengis.net/kml/2.2">'
        sp+= '<Document id="root_doc">'
        sp+= '<Folder><name>testgeom</name>'
        sl = sp*1
        #sl1 += '<Style><LineStyle><color>ff0000ff</color></LineStyle></Style>'
        sp1 = '<Placemark><PointString><altitudeMode>clampToGround</altitudeMode><coordinates>'
        sp2 = '</coordinates></PointString></Placemark>\n'
        sl1 = '<Placemark><LineString><altitudeMode>clampToGround</altitudeMode><coordinates>'
        sl2 = '</coordinates></LineString></Placemark>\n'
        lcoords = dicz['coords']
        for coo in lcoords:
            strcoo = ' '.join([str(a) for a in coo])
            if len(coo)==1: 
                ip = 1
                sp+= sp1 + strcoo.replace('(','').replace(')','')+ sp2
            else :
                il = 1
                sl+= sl1 + strcoo.replace('(','').replace(')','')+ sl2
        if ip==0 : sp = None
        if il==0 : sl = None
        return sp,sl
'''