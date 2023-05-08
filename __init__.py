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
    This function will add an entry into the PyMOL architecture

    **Parameters**

    app: *filet*
        The application to be loaded for plugin initialization
        
    **Returns**

        Menu item `PyMIV` under the `Plugins` menu item bar
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
    This function will open the custom dialog window

    **Parameters**

    None
        
    **Returns**

        Dialog window `PyMIV`
    '''
    # pylint: disable=global-statement
    global dialog

    if dialog is None:
        dialog = make_dialog()

    dialog.show()

def make_dialog():
    '''
    This function will make the dialog window if there is no window open

    **Parameters**

    None
        
    **Returns**

        Dialog window tasks `PyMIV`
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
        This function will allow selection of a valid PDB file when `Broswe` is clicked

        **Parameters**

        None
            
        **Returns**

            None
        '''
        filename = getOpenFileNames(
            dialog)
        clean_filename=(clean_file_path(str(filename)))
        form.lineEdit.setText(clean_filename)


    def disulfide_finder_button():
        '''
        This function will run the disulfide_finder function and 
        load `disulfide_bonds.pml when `Disulfide Finder` is clicked

        **Parameters**

        None
            
        **Returns**

            None
        '''
        # retreive PDB file data
        pdb_file = form.lineEdit.text()
        # Error Code
        if pdb_file == "":
            print("Please input a valid .pdb file name")
        else:
            print('User Entered Filename:', pdb_file)
            calc_disulfide(pdb_file)
            cmd.run("disulfide_bonds.pml")
            print('Yellow = Cysteine Sulfur Atoms')

    def wc_nwc_button():
        '''
        This function will run the calc_wc_nwc function and 
        load `get_bonds.pml when `WC vs Non-WC` is clicked

        **Parameters**

        None
            
        **Returns**

            None
        '''
        pdb_file = form.lineEdit.text()
        # Error Code
        if pdb_file == "":
            print("Please input a valid .pdb file name")
        else:
            print('User Entered Filename:', pdb_file)
            calc_wc_nwc(pdb_file)
            cmd.run("get_bonds.pml")
            print('Yellow = WC\nRed=Non-WC')

    def alpha_helix_button():
        '''
        This function will run alpha_helice function and 
        load `helix_bonds.pml when `Alpha Helix` is clicked

        **Parameters**

        None
            
        **Returns**

            None
        '''
        # retreive PDB file data
        pdb_file = form.lineEdit.text()
        if pdb_file == "":
            print("Please input a valid .pdb file name")
        else:
            print('User Entered Filename:', pdb_file)
            alpha_helice(pdb_file)
            cmd.run("helix_bonds.pml")

    def calc_mw_button():
        '''
        This function will run calc_peptide_mw function and 
        when `Alpha Helix` is clicked

        **Parameters**

        None
            
        **Returns**

            None
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
        This function will run end_to_end_dist function and 
        load `end_to_end.pml when `End to End Distance` is clicked

        **Parameters**

        None
            
        **Returns**

            None
        '''
        # retreive PDB file data
        pdb_file = form.lineEdit.text()
        # Error Code
        if pdb_file == "":
            print("Please input a valid .pdb file name")
        else:
            print('User Entered Filename:', pdb_file)
            end_to_end_dist(pdb_file)
            cmd.run("end_to_end.pml")

    def output_fasta_button():
        '''
        This function will run end_to_end_dist function and 
        load `end_to_end.pml when `End to End Distance` is clicked

        **Parameters**

        None
            
        **Returns**

            None
        '''
        # retreive PDB file data
        pdb_file = form.lineEdit.text()
        # Error Code
        if pdb_file == "":
            print("Please input a valid .pdb file name")
        else:
            print('User Entered Filename:', pdb_file)
            output_fasta(pdb_file)

    # To connect clicking buttons to a value, text or command
    form.browse.clicked.connect(browse_filename)
    form.done.clicked.connect(dialog.close)
    form.disulfideFinder.clicked.connect(disulfide_finder_button)
    form.calculateMW.clicked.connect(calc_mw_button)
    form.wcAndNonWC.clicked.connect(wc_nwc_button)
    form.hydrogenBond.clicked.connect(alpha_helix_button)
    form.endToEndDistance.clicked.connect(end_to_end_button)
    form.fasta.clicked.connect(output_fasta_button)

    return dialog
