##########################################
# WARNING: NO CHANGE NEEDED IN THIS FILE !
##########################################
# THIS FILE SUMMARIZES DEFINITIONS NEEDED FOR Stéphane Paulin-Henriksson (SPH) MODIFICATIONS IN THE ORTI3D PROGRAM
# THE CONFIGURATION CAN BE MADE IN THE FILE 'configSPH.py'
##########################################

import datetime as dt
import os
import enum
#from .configSPH import *
import pickle as pk

#useless imports that were commented
#from PyQt5 import QtCore
#from PyQt5 import QtWidgets
#from PyQt6.QtWidgets import QtCore, QtWidgets, QtGui
#import sys

#note pyqtgraph seems not adapted here and PyQt/matplotlib seems better
#import pyqtgraph as pg
#from pyqtgraph.Qt import QtCore, QtWidgets, QtGui

#the 'pathlib' import works only for python versions >= 3.5. ==> I dont use it because I'm not sure of this condition
# import pathlib

#############################
# definitions of enumerations
#############################
class sph_possibleReplayScheme(enum.Enum):
    """
    to indicate what is made when replaying
    it is admitted that when a backup is checked, final values are restored
    but when a backup is unchecked, it is possible to define different schemes
    """
    checkedRestoreFinalValues_uncheckedDoNothing = 1
    checkedRestoreFinalValues_uncheckedRestoreInitialValues = 2
    checkedRestoreFinalValues_uncheckedReplayHistory = 3

class sph_possibleBackupFormat(enum.Enum):
    """
    possible formats of a backup
    in 'pickle' mode, the backups are written in binary-pickle dump fmt in the file 'modifSPH2023_saveDir\fileName'
    in 'memory' mode, the backups are saved in the RAM, in a 'backupList'
    in tests, the 'pickle' mode was often preferred because more supposed to be able to deal with the graphics
    but it was not optimized because unfinished and there is a painful drawback: currently a pickle file is written for each backup.
    This means potentially 100s of files.
    finally the graphics are ignored, so now the 'memory' mode seems more convenient (and that's why pickle was unfinished)
    but pickle is still convenient for debugging
    """
    memory = 1
    pickle = 2

class sph_backupType(enum.Enum):
    """
    possible types of a backup: TBC...
    - topBar:
    - leftColumn:
    """
    topBar = 1
    leftColumn = 2

# add 240212: backupType kw
def sph_addSave(initialValues=[],finalValues=[],
                currentLine=None,modelGroup=None,
                mode=sph_possibleBackupFormat.pickle,fileName=None,
                backupList=None,overWrite=False,backupType=sph_backupType.leftColumn
                ):#preliminaryIdeaOfModel=None,currentGroup=None,
    """
    add a backup (in the backup list or among the pickle files), with initial values and final values, before and after the action respectively

    Parameters
    ----------
        initialValues : list of numbers
            values before the action
        finalValues : list of numbers
            values after the action
        currentLine : str
            name of the line that was modified by the action. For instance 'bas.2' or 'geoch.1'
        modelGroup : str
            name of the model group that was modified by the action. In the current state, may be 'Modflow series',
            'Modflow USG' or 'Min3p'
        mode : sph_possibleBackupFormat
            can be sph_possibleBackupFormat.pickle or sph_possibleBackupFormat.memory
            in 'pickle' mode the values are saved in binary-pickle dump fmt in the file 'modifSPH2023_saveDir/fileName'
            and saveList is useless
            in 'memory' mode the values are saved in the RAM, in 'backupList' (encapsulated fmt) and 'fileName' is useless
            in tests, the 'pickle' mode was often preferred because more supposed to be able to deal with the graphics
            finally the 'memory' mode seems more convenient
        overWrite : boolean
            ignored if not in 'pickle' mode
            if True, the pickle files are overwritten, if False pickle the file name is completed with '_'s
        fileName : None or 'str'
            ignored if not in 'pickle' mode
            if None the pickle file name is created otherwise it may be given (then overWrite and saveDir parameters are useless)
        backupList : list (of backups)
            ignored if not in 'memory' mode
            In the memory case, this new backup is appended to the list
        backupType: sph_backupType
            can be sph_backupType.topBar or sph_backupType.leftColumn... TBC...

    Returns
    -------
        nothing
    """
    from configSPH import modifSPH2023_isThereBackupWhenNoChange
    # add 240212: backupType kw
    currentBackup = sph_encapsulate(currentLine=currentLine,modelGroup=modelGroup,
                                    initialValues=initialValues,finalValues=finalValues,
                                    backupType=backupType
                                    )#currentGroup=currentGroup,currentModel=currentModel,

    # add 240219
    isChanged = False
    if not modifSPH2023_isThereBackupWhenNoChange:
        for count,value in enumerate(initialValues):
            if value != finalValues[count]:
                isChanged = True
                break

    if isChanged or modifSPH2023_isThereBackupWhenNoChange:
        if mode==sph_possibleBackupFormat.memory:#that's a memory save
            backupList.insert(0,(currentBackup))
        else:#that's a pickle save
            from configSPH import modifSPH2023_saveDir
            if not isinstance(fileName,str):
                # windows seems to have sometimes issues dealing with '.', so for simplicity I dont use them
                fileName = sph_buildPartFileName(data=currentBackup)#+'.pk'
                directory = os.path.join(os.getcwd(),'lib','ilibq',modifSPH2023_saveDir)
                if modifSPH2023_saveDir != '':# note the use of the pathlib module suppose python>=3.5, so I dont use this here
                    if not os.path.isdir(directory):
                        try:
                            os.makedirs(directory)
                        except:
                            # TODO: deal with this 'what the hell' case
                            pass
                    fileName = os.path.join(directory,fileName)
                if not overWrite:
                    while os.path.isfile(fileName):#+'.pk'):
                        fileName += '_'
            with open(fileName,'wb') as f:
                pk.dump(currentBackup,f,pk.HIGHEST_PROTOCOL)

# add 240212: backupType kw
def sph_encapsulate(modelGroup=None,currentLine=None,
                    initialValues=None,finalValues=None,
                    backupType=None,
                    ):#,currentModel=None):
    """
    give a convenient format to a backup

    Parameters
    ----------
        modelGroup : str
            like for sph_addSave above
        currentLine : str
            like for sph_addSave above
        initialValues : list of values
             like for sph_addSave above
        finalValues : list of values
            like for sph_addSave above
        backupType: sph_backupType
            can be sph_backupType.topBar or sph_backupType.leftColumn... TBC...

    Returns
    -------
        a dictionary representing the backup
    """
    # add 240212: backupType kw
    current = {'modelGroup':modelGroup,
               'currentLine':currentLine,
               'initialValues': initialValues,
               'finalValues': finalValues,
               'date': dt.datetime.now().strftime('%y%m%d'),
               'time': dt.datetime.now().strftime('%H%M%S'),
               'lastAction': 'debug: not implemented for now...',
               'backupType': backupType
               }#'currentModel':currentModel,

    return current

def sph_buildPartFileName(date=None,time=None,currentLine=None,data=None,modelGroup=None):#,currentGroup=None,currentModel=None):
    """
    build the name of a backup with a convenient scheme

    Parameters
    ----------
        date : None or str
            if date is a string it is used in the built name, otherwise it is considered to be in data['date'] that must exist
            reminder the default format for the date is '%y%m%d'
        time : None or str
            if date is a string it is used in the built name, otherwise it is considered to be in data['date'] that must exist
            reminder the default format is '%H%M%S'
        data : dict or None
            if any other input is not a string, data must exist with the proper fields, otherwise it is ignored
        currentLine : str
            like for sph_addSave above
        modelGroup : str
            like for sph_addSave above

    Returns
    -------
        str
    """
    if not isinstance(date, str):
        date = data['date']
    if not isinstance(time, str):
        time = data['time'].replace(':','')
    if not isinstance(currentLine,str):
        # note Windows seems sometimes annoyed with the '.' character
        # for simplicity I avoid it
        currentLine = data['currentLine'].replace('.','')
    if not isinstance(modelGroup, str):
        #for simplicity and coherence with my uses, I avoid the ' ' character and put everything in camel case
        temp = data['modelGroup'].split(' ')
        modelGroup = temp[0]
        try:
            for s in temp[1:]:
                modelGroup += s.capitalize()
        except:
            pass

    return 'orti3D_'+date+'_'+time+'_'+modelGroup+'_'+currentLine

def sph_loadSave(fileName=None,mode=sph_possibleBackupFormat.pickle,backupList=None,id=None):
    """
    load a backup from a format (currently memory mode or pickle mode)

    Parameters
    ----------
        fileName : str or None
            in the pickle mode this is the file name of the backup (and must be a str)
            in memory mode it is ignored
        mode : enum sph_possibleBackupFormat
            currently either 'pickle' or 'memory', see the 'sph_possibleBackupFormat' enum
        backupList : list of backups
            ignored in the 'pickle' mode
            TBC...
        id : int
            ignored in the 'pickle' mode
            TBC...

    Returns
    -------
        a backup
    """
    from configSPH import modifSPH2023_saveDir, modifSPH2023_debugLevel
    if modifSPH2023_debugLevel>10:
        print('debug 240212f:',mode,mode == sph_possibleBackupFormat.memory,mode == sph_possibleBackupFormat.pickle)

    # TODO: ADD WARNING
    # TODO: THE pickle LOADING PROCESS IS NOT COMPLETELY SECURED ACCORDING TO INTERNET COMMENTS !
    if mode == sph_possibleBackupFormat.memory:# then it is 'memory' mode
        #modif 240209
        if not isinstance(id,int):
            # there is a pb
            # TODO: deal with this case. for now never the case
            pass
        currentSave = backupList[id]
    elif mode == sph_possibleBackupFormat.pickle:# then it is 'pickle' mode
        with open(os.path.join(os.getcwd(),'lib','ilibq',modifSPH2023_saveDir,fileName), 'rb') as f:
            currentSave = pk.load(f)
    else:
        print('big trouble !')
        exit(0)
    return currentSave

def sph_makeBackupListForPickleMode():
    """
    list the backup directory in the 'pickle' mode (see sph_addSave above) and return the backup list
    rk: file names are assumed to be "<blabla>_<date>_<time>_<modelGroup>_<currentLine without the '.'>", this to please Windows that doesnt like the '.'
    rk2: the currentLine is assumed to be "<blabla>.<one digit>"

    Returns
    -------
        list of backups
    """
    from configSPH import modifSPH2023_saveDir
    #note the file list is anti-chronologically sorted !
    fileList = sorted(os.listdir(os.path.join(os.getcwd(), 'lib', 'ilibq', modifSPH2023_saveDir)),reverse=True)
    backupList = []
    for f in fileList:
        temp = f.split('_')
        date = temp[1]
        time = temp[2]

        # warning: this will only work for model families with less than 10 instances !!!
        # warning: because then the line is always <name>.<1 digit>
        currentLine = temp[4][:-1]+'.'+temp[4][-1]

        modelGroup = temp[3]
        backupList.append(sph_backup(modelGroup=modelGroup,
                                     currentLine=currentLine,
                                     fileName=f,
                                     date=date, time=time))
    return backupList

class sph_backup():
    """
    a convenient format to describe a backup.
    retrospectively this could have been a simple dictionary
    """
    def __init__(self,date=None,time=None,fileName=None,
                 lastAction=None,currentLine=None,modelGroup=None,
                 initValues=None,finalValues=None):
        """
        Parameters
        ----------
            date : str
                the default format is '%y%m%d' but can also be 'now'
            time : str
                the default format is '%H%M%S' but can also be 'now'
            fileName : str
            lastAction : str
            currentLine : str
            modelGroup : str
            initValues : list
            finalValues : list
        """
        if date == 'now':
            self.date = dt.datetime.now().strftime('%y%m%d')
        else:
            self.date = date

        if date == 'now' or time == 'now':
            self.time = dt.datetime.now().strftime('%H%M%S')
        else:
             self.time = time

        self.fileName = fileName
        self.lastAction = lastAction
        self.currentLine = currentLine
        self.modelGroup = modelGroup
        self.initValues = initValues
        self.finalValues = finalValues

    def __str__(self):
        return 'time='+self.time+' ; modelGroup='+self.modelGroup+' ; currentLine='+self.currentLine+\
               ' ; initialValues='+str(self.initValues)+\
               ' ; finalValues='+str(self.finalValues)

class sph_loadingMode(enum.Enum):
    """
    possible schemes of a load:
    - initial: initial values are loaded (before modif)
    - final: final values are loaded (after modif)
    - noLoading: nothing done
    originally, when a backup is checked in the replay list, I was loading final values (equiv loadingMode = 2 here)
    and when a backup was unchecked, I was doing nothing (equiv loadingMode = 3 here)
    but on 240202 OA found it was better, when unchecked, to load initial values (equiv to loadingMode = 1 here)
    """
    initial = 1
    final = 2
    noLoading = 3



def sph_tracefunc_file(frame, event, arg):
    # sometimes convenient for debugging but note this makes an open/close each time a function is called !
    # this is absolutely not optimized !
    with open('sph240212.log','at') as f:
        if event == "call":
            f.write(">>> call function "+str(frame.f_code.co_name)+"\n")
        elif event == "return":
            f.write("<<< exit function "+str(frame.f_code.co_name)+"\n")

def sph_tracefunc_print(frame, event, arg):
    if event == "call":
        print(">>> call function",str(frame.f_code.co_name))
    elif event == "return":
        print("<<< exit function ",str(frame.f_code.co_name))

def sph_constructBackupList(backupList=[]):
    """
    to build the backup list, depending on the format (currently 'pickle' or 'memory')
    note the list is sorted anti-chronologically
    in memory mode this is already the case of the list in memory, since

    """
    from configSPH import modifSPH2023_backupFormat, modifSPH2023_saveList, modifSPH2023_debugLevel
    if backupList=='dev':# or isinstance(backupList,bool):
    # this is a fast and convenient debug case that could be remove now
        backupList = [sph_backup(fileName='toto',date='now'),
                      sph_backup(fileName='tata',date='now'),
                      sph_backup(fileName='tutu',date='now')]
    elif modifSPH2023_backupFormat == sph_possibleBackupFormat.pickle:
        backupList = sph_makeBackupListForPickleMode()
    elif modifSPH2023_backupFormat == sph_possibleBackupFormat.memory:
        backupList = []
        for backup in modifSPH2023_saveList:
            backupList.append(sph_backup(modelGroup=backup['modelGroup'],
                                         currentLine=backup['currentLine'],
                                         initialValues=backup['initialValues'],
                                         finalValues=backup['finalValues'],
                                         date=backup['date'], time=backup['time']))

    if modifSPH2023_debugLevel > 20:
        print('debug 240209b:',backupList,type(backupList))

    return backupList

def sph_manageBackupDirectory():
    from configSPH import modifSPH2023_saveDir, modifSPH2023_cleanSaveDirAtStart, modifSPH2023_backupFormat
    backupDirectory = os.path.join(os.getcwd(), 'lib', 'ilibq', modifSPH2023_saveDir)
    if os.path.isfile(backupDirectory):
        # TODO: pb ! what do we do ?
        exit(0)
    if not os.path.isdir(backupDirectory):
        os.mkdir(backupDirectory)
    if modifSPH2023_cleanSaveDirAtStart and modifSPH2023_backupFormat == sph_possibleBackupFormat.pickle:
        # once more the 'pathlib' module works only for a recent python and I am not sure of this condition,
        # that's why I use 'os.remove'
        for tuple in os.walk(backupDirectory):
            # 'os.walk' is a convenient generator that cleanly separates files in a tuple
            # I use on purpose this syntax with 'for...pass' that create tuple
            # to get a proper tuple instead of a 'not subscriptable' error message
            # note this typical python trick is not stylish. maybe there is a more pythonistic way...
            pass
        for file in tuple[2]:
                os.remove(os.path.join(backupDirectory, file))
