"""
few options linked to modifSPH2023
"""
from defSPH import *

modifSPH2023 = trySPH_enableTrySPH = True
modifSPH2023_saveList = trySPH_savelist = []
modifSPH2023_debugLevel = trySPH_debugLevel = 0
modifSPH2023_hideVerticalHeaderInBackupTable = trySPH_hideVerticalHeaderInBackupTable = True
modifSPH2023_useInspect = trySPH_inspect = False
modifSPH2023_cleanSaveDirAtStart = True
# this option to avoid backups when there is no change does work in memory mode only
modifSPH2023_isThereBackupWhenNoChange = False

modifSPH2023_replayScheme = sph_possibleReplayScheme.checkedRestoreFinalValues_uncheckedReplayHistory
modifSPH2023_backupFormat = sph_possibleBackupFormat.pickle
#modifSPH2023_backupFormat = sph_possibleBackupFormat.memory

if modifSPH2023_backupFormat == sph_possibleBackupFormat.pickle:
    modifSPH2023_usePickle = trySPH_usePickle = True
    modifSPH2023_saveDir = trySPH_saveDir = 'Save'
else:
    modifSPH2023_usePickle = trySPH_usePickle = False
