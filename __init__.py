# To initialize the plugin to be installed by PyMOL
from __future__ import absolute_import
from __future__ import print_function

# To provide an entry point to PyMOL's API
import os
from pymol import cmd
from pymol.Qt import QtWidgets
from pymol.Qt.utils import loadUi
from pymol.Qt.utils import getSaveFileNameWithExt

# To load the UI file into our dialog
from pymol.Qt.utils import loadUi

#This is a series of functions for mathematical calculations to read into PDB files and output a .pml text to import into the PyMOL



import math
from calc_interaction import *
def __init_plugin__(app=None):
    '''
    Add an entry to the PyMOL "Plugin" menu
    '''
    from pymol.plugins import addmenuitemqt
    addmenuitemqt('PyMIV', run_plugin_gui)

# To create global reference of the dialog variables
dialog = None

# to give filename of the UI file
uifile = os.path.join(os.path.dirname(__file__), 'pymolGUI.ui')

# To load the UI dialog
form = loadUi(uifile, dialog)


def run_plugin_gui():
    '''
    Open the custom dialog
    '''
    global dialog

    if dialog is None:
        dialog = make_dialog()

    dialog.show()


def make_dialog():
    '''
    Make the dialog window
    '''
    
    # To create a new UI window
    dialog = QtWidgets.QDialog()

    # To populate the Window from our .ui file
    uifile = os.path.join(os.path.dirname(__file__), 'pymolGUI.ui')
    form = loadUi(uifile, dialog)
   

    def disulfideFinderButton():
        '''
        Run the actions in the dialog window
        '''
        # retreive PDB file data
        pdb_file = form.lineEdit.text()
        print('User Entered Filename:', pdb_file)
        calc_disulfide(pdb_file)

        # To debug code
        
    def wcAndNonWCButton():
        '''
        Run the actions in the dialog window
        '''
        # retreive PDB file data
        pdb_file = form.lineEdit.text()
        print('User Entered Filename:', pdb_file)
        (pdb_file)
        calc_WC_and_NonWC(pdb_file)

    def alphaHeliceButtonButton():
        '''
        Run the actions in the dialog window
        '''
        # retreive PDB file data
        pdb_file = form.lineEdit.text()
        print('User Entered Filename:', pdb_file)
        (pdb_file)
        alpha_helice(pdb_file)

    # To connect clicking buttons to a value, text or command
    form.disulfideFinder.clicked.connect(disulfideFinderButton)
    #form.calculateMW.clicked.connect(calculateMWButton)
    form.wcAndNonWC.clicked.connect(wcAndNonWCButton)
    form.hydrogenBond.clicked.connect(alphaHeliceButtonButton)

    return dialog