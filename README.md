
# PyMIV
PyMIV: A Pymol Plugin as a Molecular Interactions Viewer (MIV) of disulfide bonds, WC/Non-WC bonds, protein secondary structure etc.

## Prerequisites
* [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/)
OR
* [homebrew] (https://brew.sh/)

* [PyMOL-open-source-vewier] (git@github.com:schrodinger/pymol-open-source.git)

* [PyMOL-open-source-vewier] (https://anaconda.org/conda-forge/pymol-open-source)

* [wxPython] (https://formulae.brew.sh/formula/wxpython)




### Adding PyMIV to PyMOL

In PyMOL window, go to `Plugin` -> `Plugin manager` -> `Install new plugin`, choose `PyMIV.zip` and click `OK` on the next step. You will then see a pop-up message `Plugin "PyMIV" has been installed`.


## Using PyMIV


#### Open PyMOL

Run PyMIV through `Plugin` -> `PyMIV`


### PyMIV Tabs

PyMIV has six tabs and one `Done` button. 

* The first tab `Browse` is to select a valid PDB format file anywhere on the user's computer. Access peptide and nucleic acid 
structures in PDB format at RCSB Protein Data Bank (https://www.rcsb.org/)

* The second tab `Disulfide Finder` includes the command for finding disulfide bodns in a polypeptide structure. This works for peptides

* The third tab `End to End Distance` will calculate the end to end distance of peptide and load the structure into the viewer. A measurement distance between the C-alpha carbon of the first and last amino acid will be displayed. This works for peptides

* The fourth tab `WC vs Non-WC` performs measurements to find nucleotide interactions in nucleic acids that abide watson-cick(e.g. A-T and G-C) interactions and non-watson-crick interactions(e.g.g G-G and U-C), displays the results in the PyMOL Viewer. This works for both DNA and RNA

* The fifth tab `Alpha Helix` will generate a FASTA sequence of a peptide and display with an `H` under amino acid resiues that are in a alpha helix secondary structure. The PDB structure will be loaded into the viewer and the alpha helices will be displayed in pink

* The sixth tab `Calculate MW` calculates the molecular weight of the given molecule. This works for peptides


## Author Notice

* Created by Maranda McDonald and published on May 7 2023

* Insipitation taken from `speleo3 ` `pymol2-demo-plugin` from the PyMOLWiki (https://pymolwiki.org/index.php/Plugins_Tutorial) 


