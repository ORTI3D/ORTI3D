
class genericDialog:
    def __init__(self, gui, title, data, helpstring=''):
        self.gui = gui
        if gui.gtyp=='wx': 
            self.dialg = wxGenericDialog(gui, title, data, helpstring='')
        elif gui.gtyp=='qt':
            self.dialg = qtGenericDialog(gui, title, data, helpstring='')
            
    def getValues(self):
        if self.gui.gtyp=='wx':
            retour = self.dialg.ShowModal()
            if retour == wx.ID_OK: return self.dialg.getValues()
            else : return None
        elif self.gui.gtyp=='qt':
            retour = self.dialg.exec_()
            if retour == 1: 
                val = self.dialg.getValues()
                return val
            else : return None