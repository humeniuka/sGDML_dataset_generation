#!/usr/bin/env python
"""
functions for submitting QChem jobs to the queue,
requires the script 'run_qchem.sh'
"""
import numpy as np
import os
import subprocess
from collections import OrderedDict

from sgdml_dataset_generation.readers.qchem import QChemOutputFile

def run_qchem(atoms, script='grad.in', directory=".", nprocs=1, mem="6Gb"):
    """
    run QChem input script and read results from output file

    Parameters
    ----------
    atoms:  ase.atoms.Atoms
      molecular geometry
    directory: str
      directory where the calculation should be performed

    Optional
    --------
    nprocs : int
      number of processors
    mem    : str
      allocated memory (e.g. '6Gb', '100Mb')

    Returns
    -------
    results : QChemOutputFile
      dictionary with results of calculation, which can be accessed by the following keys
        'geometry (au)'             -  cartesian geometry
        'energies (au)'             -  dictionary with excitation energies
        'gradients (au)'            -  dictionary with gradients for each state
        'derivative couplings (au)' -  dictionary with NAC vectors for combinations of states
    """
    # create directory if it does not exist already
    os.system("mkdir -p %s" % directory)
    output = script.replace(".in", ".out")
    os.system(f"cp {script} {directory}/{script}")
    # update geometry
    geom_lines = []

    for atom in atoms:
        l = "%2s    %+12.10f   %+12.10f   %+12.10f \n" % (atom.symbol, atom.x, atom.y, atom.z)
        geom_lines.append(l)

    with open(qchem_file) as fh:
        lines = fh.readlines()

    # excise old geometry block
    start = None
    end = None
    for i,l in enumerate(lines):
        if "$molecule" in l:
            start = i
        elif "$end" in l and not (start is None):
            end = i
            break
    # replace old geometry by new one, the charge and multiplicity
    # in line start+1 is left as is.
    lines = lines[0:start+2] + geom_lines + lines[end:]

    with open(qchem_file, "w") as fh:
        for l in lines:
            fh.write(l)

    # calculate electronic structure
    #print "running QChem ..."
    # submit calculation to the cluster
    ret  = os.system(f"cd {directory}; run_qchem.sh --wait {script} {nprocs} {mem}")
    if ret != 0:
        # Since the temporary files from each image are deleted, it is very difficult to
        # figure out why a calculation failed. Therefore the content of the log-file
        # is printed if the calculation failed.
        print(f" ****** content of log-file {directory}/{output} ****** " % directory)
        os.system(f"cat {directory}/{output}")
        print(f" ****** end of log-file ****** ")
        print(f"Return status = {ret}, error in QChem calculation, see error messages above !")

    with open(output) as f:
        results = QChemOutputFile(f)

    return results

def get_calculator(name):
    """
    retrieve function for calculating electronic structure
    """
    if name == "qchem":
        return run_qchem
    else:
        raise ValueError("Unknown calculator '%s'" % name)
