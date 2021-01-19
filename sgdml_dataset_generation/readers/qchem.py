#!/usr/bin/env python
# -*- coding: utf-8 -*-

# # Imports
import numpy as np
from collections import OrderedDict
import logging

# # Logging
logger = logging.getLogger(__name__)
logging.basicConfig(format="[%(module)-12s] %(message)s", level=logging.INFO)

from sgdml_dataset_generation import units

class QChemOutputFile(object):
    """
    parse output file produced by the quantum chemistry program QChem.

    Parameters
    ----------
    f    :   File
      file handle opened for reading.
      The user has to ensure the file handle is opened and closed at the end.

    The data records can be accessed by their names (see example below).

    Example
    -------

      >>> with open("qchem.out") as f:
      >>>   qchem = QChemOutputFile(f)
      >>>   # show names of all records
      >>>   print(qchem.keys())
      >>>   # element names
      >>>   print(qchem['symbols'])

    """
    def __init__(self, f):
        self.filename = f.name
        self.data = OrderedDict()
        self.data['ok'] = False
        self.data['gradients (au)'] = {}
        self.data['derivative couplings (au)'] = {}
        
        for line in f:
            parts = line.split()

            if "Standard Nuclear Orientation (Angstroms)" in line:
                f.readline()
                f.readline()
                f.readline()
                symbols = []
                geometry = []
                while True:
                    line = f.readline()
                    if "-----" in line:
                        break
                    parts = line.split()
                    symbols.append(parts[1])
                    pos_xyz = list(map(float, parts[2:5]))
                    geometry.append(pos_xyz)
                geometry = np.array(geometry)
                self._getset('symbols', symbols)
                self._getset('geometry (au)', geometry / units.bohr_to_angs)

            elif "Total energy in final basis set" in line:
                parts = line.split()
                self._getset('scf energy', float(parts[8]))
            
            elif "TDDFT Excitation Energies" in line:
                exc_energies = self._getset('excitation energies (eV)', {})
                osc_strengths = self._getset('oscillator strengths', {})
                f.readline()
                f.readline()
                while True:
                    line = f.readline()
                    parts = line.split()
                    if "-----" in line:
                        break
                    if "Excited state" in line:
                        state = int(parts[2].replace(":", ""))
                        exc_energies[state] = float(parts[7])
                    elif "Strength" in line:
                        osc_strengths[state] = float(parts[2])
                # total energies of all states
                energies = self._getset('energies (au)', {})
                energies[0] = self['scf energy']
                for state, en_ex in self['excitation energies (eV)'].items():
                    energies[state] = energies[0] + en_ex / units.hartree_to_ev
            elif "between states " in line:
                derivative_couplings = self._getset('derivative couplings (au)', {})
                parts = line.split()
                I, J = int(parts[2]), int(parts[4])
                while True:
                    line = f.readline()
                    if "DC between" in line and "with ETF" in line:
                        break
                for i in range(2):
                    f.readline()
                nacv = []
                while True:
                    line = f.readline()
                    parts = line.split()
                    if "-----" in line:
                        break
                    nacv_i = list(map(float, parts[1:4]))
                    nacv.append(nacv_i)
                nacv = np.array(nacv)
                
                derivative_couplings[(I,J)] = nacv

            elif "Gradient of SCF Energy" in line:
                grad_scf = []
                while True:
                    line = f.readline()
                    if "Max gradient" in line:
                        break
                    parts = line.split()
                    atom_indices = list(map(lambda p: int(p)-1, parts))
                    grad_xyz = []
                    for xyz in [0,1,2]:
                        line = f.readline()
                        parts = line.split()
                        dEdx = list(map(lambda p: float(p), parts[1:]))
                        grad_xyz.append( dEdx )
                    grad_xyz = list(np.array(grad_xyz).transpose())
                    grad_scf += grad_xyz
                grad_scf = np.array(grad_scf)
                gradients = self._getset('gradients (au)', {})
                gradients[0] = grad_scf

            elif "RPA" in line and "State Energy is" in line:
                state = int(line.split()[1])
                
            elif "Gradient of the state energy" in line:
                grad_ex = []
                while True:
                    line = f.readline()
                    if "Gradient time" in line:
                        break
                    parts = line.split()
                    atom_indices = list(map(lambda p: int(p)-1, parts))
                    grad_xyz = []
                    for xyz in [0,1,2]:
                        line = f.readline()
                        parts = line.split()
                        dEdx = list(map(lambda p: float(p), parts[1:]))
                        grad_xyz.append( dEdx )
                    grad_xyz = list(np.array(grad_xyz).transpose())
                    grad_ex += grad_xyz
                grad_ex = np.array(grad_ex)
                gradients = self._getset('gradients (au)', {})
                gradients[state] = grad_ex

            elif "Thank you very much for using Q-Chem." in line:
                self.data['ok'] = True
                
    def _getset(self, key, default):
        item = self.data.get(key, default)
        self.data[key] = item
        return item

    def __getitem__(self, key):
        """
        access data fields by their names

        Parameters
        ----------
        key     :   str
          name of field that should be retrieved (e.g. 'gradients (au)')

        Returns
        -------
        field   :  float, int or ndarray
          a KeyError is raised if the field is not present
        """
        return self.data[key]
    def keys(self):
        """
        list names of all fields read

        Returns
        -------
        keys  :  list of str
          field names
        """
        return self.data.keys()

