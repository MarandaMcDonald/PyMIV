# This initializes the plugin to be installed by PyMOL
from __future__ import absolute_import
from __future__ import print_function

import os


def __init_plugin__(app=None):
    '''
    Add an entry to the PyMOL "Plugin" menu
    '''
    from pymol.plugins import addmenuitemqt
    addmenuitemqt('Demo "Render" Plugin', run_plugin_gui)

# global reference to avoid garbage collection of our dialog
dialog = None


def run_plugin_gui():
    '''
    Open our custom dialog
    '''
    global dialog

    if dialog is None:
        dialog = make_dialog()

    dialog.show()

# filename of the UI file
uifile = os.path.join(os.path.dirname(__file__), 'input.ui')

# load the UI file into our dialog
from pymol.Qt.utils import loadUi
form = loadUi(uifile, dialog)

def make_dialog():
    # entry point to PyMOL's API
    from pymol import cmd

    # pymol.Qt provides the PyQt5 interface, but may support PyQt4
    # and/or PySide as well
    from pymol.Qt import QtWidgets
    from pymol.Qt.utils import loadUi
    from pymol.Qt.utils import getSaveFileNameWithExt

    # create a new Window
    dialog = QtWidgets.QDialog()

    # populate the Window from our *.ui file which was created with the Qt Designer
    uifile = os.path.join(os.path.dirname(__file__), 'input.ui')
    form = loadUi(uifile, dialog)

    
    return dialog