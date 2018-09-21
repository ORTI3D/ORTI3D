import wx
#import wx.grid as wgrid
import os
import core as corebase
from valueDialog import *
from addin import *
from menus import *

class modifyModelGui(wx.Frame):
    """ main frame"""
    def __init__(self,modName):
        wx.Frame.__init__(self,None,1,title='iQpht3d',size=(300,500))
        #self.Maximize(True)
        self.fDir,self.fName='',''
        core = corebase.Core(gui=self)
        self.core,self.modName = core,modName

        # create addin
        self.addin = core.addin
        self.addin.setGui(self)
        # panels on the left
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        self.valueGui = valueDialog(self,core,modName)
        mainSizer.Add(self.valueGui,1,wx.EXPAND)
        self.SetSizer(mainSizer)

        #Menus
        menus = Menus(self,core)
        menuFile = wx.Menu()
        ne = menuFile.Append(-1,"&New \tCTRL+n")
        menuFile.AppendSeparator()
        op = menuFile.Append(wx.ID_OPEN, "&Open \tCTRL+o")
        sa = menuFile.Append(wx.ID_SAVE, "&Save \tCTRL+s")
        ss = menuFile.Append(wx.ID_SAVEAS, "&Save as")
        menuFile.AppendSeparator()
        qu = menuFile.Append(wx.ID_EXIT, "&Quit\tCTRL+q")
        self.Bind(wx.EVT_MENU, menus.OnNew,ne)
        self.Bind(wx.EVT_MENU, menus.OnOpen,op)
        self.Bind(wx.EVT_MENU, menus.OnSave,sa)
        self.Bind(wx.EVT_MENU, menus.OnSaveAs,ss)
        self.Bind(wx.EVT_MENU, menus.OnExit,qu)
        
        menuImport = wx.Menu()
        mi1=menuImport.Append(-1,'Modflow')
        self.Bind(wx.EVT_MENU, menus.OnImportModflow,mi1)

        menuRun = wx.Menu()
        mW = menuRun.Append(-1,'Write')
        mR = menuRun.Append(-1,'Run')
        self.Bind(wx.EVT_MENU, menus.OnWrite, mW)
        self.Bind(wx.EVT_MENU, menus.OnRun, mR)
        
        self.menuAddin = wx.Menu()
        self.addin.initMenus()
        
        self.menuBar = wx.MenuBar()
        self.menuBar.Append(menuFile, "&File")
        self.menuBar.Append(menuImport, "&Import")
        self.menuBar.Append(menuRun, "&Run")
        self.menuBar.Append(self.menuAddin, "&Addins")
        self.SetMenuBar(self.menuBar)
        wx.EVT_CLOSE(self,menus.OnExit)


