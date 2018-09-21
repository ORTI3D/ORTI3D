from distutils.core import setup
import py2exe

from distutils.filelist import findall
import os.path
import matplotlib

## run with python setup2exe.py py2exe and include mrsvc90.dll in ilibq
######################################################
#
#  creation de la liste de donnees de l'application
#
######################################################

# liste des fichiers de donnees necessaires a l'application
all_data_files = []

# recuperation des fichiers de donnees de matplotlib (beaucoup sont inutiles, on peut les virer a la main)
matplotlibdatadir = matplotlib.get_data_path()
matplotlibdata = findall(matplotlibdatadir)
for f in matplotlibdata:
    dirname = os.path.join('matplotlibdata', f[len(matplotlibdatadir)+1:])
    all_data_files.append((os.path.split(dirname)[0], [f]))

# recuperation des fichiers de donnees de ipht3d 
utilsdata = findall('../utils')
for f in utilsdata:
    dirname = os.path.join('../utils', f[len(matplotlibdatadir)+1:])
    all_data_files.append((os.path.split(dirname)[0], [f]))
bindata = findall('../bin')
for f in bindata:
    dirname = os.path.join('../bin', f[len(matplotlibdatadir)+1:])
    all_data_files.append((os.path.split(dirname)[0], [f]))
docdata = findall('../doc')
for f in bindata:
    dirname = os.path.join('../doc', f[len(matplotlibdatadir)+1:])
    all_data_files.append((os.path.split(dirname)[0], [f]))


######################################################
#
#               execution de py2exe
#
######################################################
#includes = ['numpy','scipy','numpy.core','wx','wx.*',
#    'traits','traits.ui','traits.ui.wx','traits.ui.wx.*',
#    'pyface.*','pyface.ui.*','pyface.ui.wx.*','pyface.ui.wx.action.*',
#    'pyface.ui.wx.timer.*','pyface.ui.wx.wizard.*','pyface.ui.wx.workbench.*',
#    'enable','enable.*','enable.drawing','enable.image',
#    'enable.image_title','enable.traits','enable.wx_backend']
 
setup(
    # nom du fichier d'execution
    windows=["orti3dv2.pyw"],
    options={
    'py2exe': {
        # packages inclus
        'packages' : ['matplotlib', 'pytz', 'scipy','pyface','mayavi'],
        #'includes' : includes,
        # packages exclus
        'excludes' : ['_gtkagg','tcl',  #'numpy',
            'pywin.debugger','pywin.debugger.dbgcon','pywin.dialogs',
            'curses','email','distutil','readline','setuptools'], #,'logging'
        # dll exclus (sinon faut les mettre dans le dossier de base)
    'dll_excludes' : ['libgdk_pixbuf-2.0-0.dll', 'libgobject-2.0-0.dll',
            'libgdk-win32-2.0-0.dll',
            'api-ms-win-core-heap-l2-1-0.dll',
            'api-ms-win-core-heap-l1-2-0.dll',
            'api-ms-win-core-io-l1-1-1.dll',
            'api-ms-win-core-libraryloader-l1-2-0.dll',
            'api-ms-win-core-libraryloader-l1-2-1.dll',
            'api-ms-win-core-sysinfo-l1-2-1.dll',
            'api-ms-win-core-synch-l1-2-0.dll',
            'api-ms-win-core-errorhandling-l1-1-1.dll',
            'api-ms-win-core-com-l1-1-1.dll',
            'api-ms-win-core-memory-l1-1-2.dll',
            'api-ms-win-core-version-l1-1-0.dll',
            'api-ms-win-core-version-l1-1-1.dll',
            'api-ms-win-core-localization-l1-2-1.dll',
            'api-ms-win-security-base-l1-2-0.dll',
            'api-ms-win-eventing-provider-l1-1-0.dll',
            'api-ms-win-core-string-l1-1-0.dll',
            'api-ms-win-core-string-l2-1-0.dll',
            'api-ms-win-core-registry-l1-1-0.dll',
            'api-ms-win-core-profile-l1-1-0.dll',
            'api-ms-win-core-processthreads-l1-1-2.dll',
            'api-ms-win-core-handle-l1-1-0.dll',
            'api-ms-win-core-file-l1-2-1.dll',

            ]
        }
    },
    data_files=matplotlib.get_py2exe_datafiles() #all_data_files
)

