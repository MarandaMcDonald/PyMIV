
# PyMIV
PyMIV: A Pymol Plugin as a Molecular Interactions Viewer (MIV) of disulfide bonds, WC/Non-WC bonds, protein secondary structure etc.

## Prerequisite
* [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/)
OR
* [homebrew] (https://brew.sh/)

* [PyMOL-open-source-vewier] (git@github.com:schrodinger/pymol-open-source.git)

* [PyMOL-open-source-vewier] (https://anaconda.org/conda-forge/pymol-open-source)

* [wxPython] (https://formulae.brew.sh/formula/wxpython)


### Adding PyMIV to PyMOL

In PyMOL window, go to `Plugin` -> `Plugin manager` -> `Install new plugin`, choose `PyMIV` zipped file and click `OK` on the next step. You will then see a pop-up message `Plugin "PyMIV" has been installed`.


### Using PyMIV


#### Open PyMOL


#### Load Object
Load the object to be analyzed in PyMOL, e.g. `1fdl.pdb` provided in `PDB_Files/` folder, by typing `load PDB_Files/1fdl.pdb` in pymol command line or dragging the file into PyMOL window.

Be sure to use PyMOL Command Terminal to cd into the directory containing the PDB files of interest


#### Run PyMIV
Run PyMIV through `Plugin` -> `PyMIV`


#### Analysis in PyMIV
- Change the `PyMOL selection/ object` entry to the name of your target object, e.g. `1fdl.pdb`



### PyMIV Tabs

Pyshifts has four tabs and one `Exit` button. 

* The first tab `Disulfide Finder` includes the command for finding disulfide bodns in a polypeptide structure 

* The second tab `WC/Non-Wc` performs measurements to find nucleotide interactions in nucleic acids that abide watson-cick(e.g. A-T and G-C) interactions and non-watson-crick interactions(e.g.g G-G and U-C), displays the results in the PyMOL Viewer

* The third tab `Calculate MW` calculates the molecular weight of the given molecule

* The fourth tab `Hydrogen Bond` can display hydrogen bonding in protein structure


## Author Notice

* Created by Maranda McDonald and published on May 7 2023


