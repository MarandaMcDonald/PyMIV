# To initialize the plugin to be installed by PyMOL
from __future__ import absolute_import
from __future__ import print_function

# To provide an entry point to PyMOL's API
import os
import math
#must do .module import
from .calc_miv import *
# pylint: disable=wrong-import-order
from pymol import cmd
from pymol.Qt import QtWidgets
from pymol.Qt.utils import loadUi
from pymol.Qt.utils import *
from pymol.plugins import addmenuitemqt
import pymol._gui
Qt = QtCore.Qt
QFileDialog = QtWidgets.QFileDialog
getOpenFileNames = QFileDialog.getOpenFileNames

############################################################
###################  Initialize Plugin #####################
############################################################

# pylint: disable=unused-argument
def __init_plugin__(app=None):
    '''
    Add an entry to the PyMOL "Plugin" menu
    '''
    addmenuitemqt('PyMIV', run_plugin_gui)

# To create global reference of the dialog variables
# pylint: disable=invalid-name
dialog = None

# to give filename of the UI file
uifile = os.path.join(os.path.dirname(__file__), 'PyMIV_GUI.ui')

# To load the UI dialog
form = loadUi(uifile, dialog)


def run_plugin_gui():
    '''
    Open the custom dialog
    '''
    # pylint: disable=global-statement
    global dialog

    if dialog is None:
        dialog = make_dialog()

    dialog.show()

def make_dialog():
    '''
    Make the dialog window
    '''

    # To create a new UI window
    # pylint: disable= redefined-outer-name
    dialog = QtWidgets.QDialog()

    # To populate the Window from our .ui file
    uifile = os.path.join(os.path.dirname(__file__), 'PyMIV_GUI.ui')
    form = loadUi(uifile, dialog)

    # To create buttons in the GUI

    def browse_filename():
        '''
        Broswe files
        '''
        filename = getOpenFileNames(
            dialog)
        clean_filename=(clean_file_path(str(filename)))
        form.lineEdit.setText(clean_filename)


    def disulfide_finder_button():
        '''
        Run the actions in the Disulfide Finder button
        '''
        # retreive PDB file data
        pdb_file = form.lineEdit.text()
        # Error Code
        if pdb_file == "":
            print("Please input a valid .pdb file name")
        else:
            print('User Entered Filename:', pdb_file)
            calc_disulfide(pdb_file)
            cmd.run("PML_Files/disulfide_bonds.pml")
            print('Yellow = Cysteine Sulfur Atoms')

    def wc_nwc_button():
        '''
        Run the actions in the WC vs Non-WC button
        '''
        # retreive PDB file data
        pdb_file = form.lineEdit.text()
        # Error Code
        if pdb_file == "":
            print("Please input a valid .pdb file name")
        else:
            print('User Entered Filename:', pdb_file)
            calc_wc_nwc(pdb_file)
            cmd.run("PML_Files/get_bonds.pml")
            print('Yellow = WC\nRed=Non-WC')

    def alpha_helix_button():
        '''
        Run the actions in the Alpha Helix Button
        '''
        # retreive PDB file data
        pdb_file = form.lineEdit.text()
        if pdb_file == "":
            print("Please input a valid .pdb file name")
        else:
            print('User Entered Filename:', pdb_file)
            alpha_helice(pdb_file)
            cmd.run("PML_Files/helix_bonds.pml")

    def calc_mw_button():
        '''
        Run the actions in the calculate peptide mw button
        '''
        # retreive PDB file data
        pdb_file = form.lineEdit.text()
        # Error Code
        if pdb_file == "":
            print("Please input a valid .pdb file name")
        else:
            print('User Entered Filename:', pdb_file)
            calc_peptide_mw(pdb_file)
            cmd.load(pdb_file)

    def end_to_end_button():
        '''
        Run the actions in the end to end distance button
        '''
        # retreive PDB file data
        pdb_file = form.lineEdit.text()
        # Error Code
        if pdb_file == "":
            print("Please input a valid .pdb file name")
        else:
            print('User Entered Filename:', pdb_file)
            end_to_end_dist(pdb_file)
            cmd.run("PML_Files/end_to_end.pml")

    # To connect clicking buttons to a value, text or command
    form.browse.clicked.connect(browse_filename)
    form.done.clicked.connect(dialog.close)
    form.disulfideFinder.clicked.connect(disulfide_finder_button)
    form.calculateMW.clicked.connect(calc_mw_button)
    form.wcAndNonWC.clicked.connect(wc_nwc_button)
    form.hydrogenBond.clicked.connect(alpha_helix_button)
    form.endToEndDistance.clicked.connect(end_to_end_button)

    return dialog
