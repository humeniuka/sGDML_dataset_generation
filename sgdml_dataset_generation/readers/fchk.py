#!/usr/bin/env python
# -*- coding: utf-8 -*-
__all__ = ["FormattedCheckpointFile"]

# # Imports
import numpy as np
import scipy.linalg as sla
from collections import OrderedDict
import re
import logging

# # Local Imports
from sgdml_dataset_generation import units
from sgdml_dataset_generation.units import hbar

# # Logging
logger = logging.getLogger(__name__)
logging.basicConfig(format="[%(module)-12s] %(message)s", level=logging.INFO)

class FormattedCheckpointFile(object):
    """
    reads all fields from formatted checkpoint files produced by the quantum chemistry
    programs Gaussian 16 and QChem.

    Parameters
    ----------
    f    :   File
      file handle opened for reading a formatted checkpoint file
      The user has to ensure the file handle is opened and closed at the end.

    The fields of the checkpoint file can be accessed by their names (see example below).
    Array fields are stored as 1D numpy arrays of float (R) or integer (I) type.

    Example
    -------

      >>> with open("freq.fchk") as f:
      >>>   fchk = FormattedCheckpointFile(f)
      >>>   print(fchk["Number of atoms"])

    """
    def __init__(self, f):
        self.filename = f.name
        self.data = OrderedDict()
        # accumulate all lines belonging to the same field (whithout newlines)
        acc = ""
        dtype = None
        for line_number, line in enumerate(f.readlines()):
            # count lines starting from 1
            line_number += 1
            # The name of a field starts in the first column and with a capital letter
            if re.match(r"^[A-Z].*", line):
                if len(acc) > 0 and not dtype is None:
                    # All lines belonging to the previous field must have been read,
                    # so we convert it to a numpy array.
                    try:
                        if dtype == str:
                            self.data[field] = acc
                        else:
                            # numerical types
                            array = np.fromstring(acc, dtype=dtype, sep=" ")
                            assert len(array) == count
                            self.data[field] = array
                    except (ValueError,AssertionError) as err:
                        logger.warning(f"A problem occurred reading field `{field}` in line {line_number:10} in {f.name} .")
                        logger.warning(err)
                        self.data[field] = np.zeros(count, dtype=dtype)
                # reset accumulator
                acc = ""

                try:
                    if len(line) < 43:
                        # skip title and method
                        logger.debug(f"skipping line {line_number:10} in {f.name}: `{line.strip()}`")
                        continue
                    # First 43 columns are reserved for the field name
                    field = line[0:43].strip()
                    logger.debug(f"field `{field}` encountered")
                    # Colum 43 contains a character indicating the data type:
                    #   I -> integer
                    #   R -> real
                    type_char = line[43]
                    if type_char == "I":
                        dtype = int
                    elif type_char == "R":
                        dtype = float
                    elif type_char == "C":
                        dtype = str
                    else:
                        dtype = None
                        # skip lines without I or R data type markers
                        logger.debug(f"skipping line {line_number:10} in {f.name}: `{line.strip()}` .")
                        continue
                    # If column 47-48 contain the string "N=",  we are dealing with an array
                    # and the last integer indicates the number of elements
                    if line[47:49] == "N=":
                        count = int(line[49:])
                    else:
                        # scalar value
                        self.data[field] = dtype(line[49:])
                except Exception as err:
                    logger.error(f"An error occurred while reading line {line_number:10} in {f.name} .")
                    raise err
            else:
                acc += " " + line
        # read last field
        if len(acc)  > 0:
            self.data[field] = np.fromstring(acc, dtype=dtype, sep=" ")
            assert len(self.data[field]) == count    
    def __getitem__(self, key):
        """
        access data fields by their names

        Parameters
        ----------
        key     :   str
          name of field that should be retrieved (e.g. 'Number of atoms')

        Returns
        -------
        field   :  float, int or ndarray
          a KeyError is raised if the field is not present in the formatted checkpoint file
        """
        return self.data[key]
    def keys(self):
        """
        list names of all fields present in the formatted checkpoint file

        Returns
        -------
        keys  :  list of str
          field names
        """
        return self.data.keys()
    def harmonic_approximation(self):
        """
        extract the position, gradient and Hessian of the potential energy in cartesian coordinates 
        
        The potential is expanded to second order around the current position x0:

             E(x) = E(x0) + grad(E)^T.(x-x0) + 1/2 (x-x0)^T . hess(E) . (x-x0)

        A frequency calculation has to be present in the formatted checkpoint file.
        The frequency calculation should be performed in a separate Gaussian 16 job using the
        following route line for the ground state calculation:

             #P functional/basis  Freq NoSymm IOp(7/32=5)

        and the following route line for an excited state frequency calculation:

             #P functional/basis TD=(Nstates=2, Root=1, NAC) Freq NoSymm IOp(7/32=5)

        Returns
        -------
        pos    :  ndarray (3*nat,)
          cartesian coordinates x0
        energy :  ndarray (1,)
          total energy E(x0) of state of interest (in Hartree) 
        grad   :  ndarray (3*nat,)
          cartesian gradient dE/dx(x0)  (in Hartree/bohr)
        hess   :  ndarray (3*nat,3*nat)
          cartesian force constants d^2E/(dxdx)(x0)  (in Hartree/bohr^2)
        """
        try:
            nat = self.data["Number of atoms"]
            # total energy of state of interest
            energy = np.array(self.data["Total Energy"])
            # geometry
            pos = self.data["Current cartesian coordinates"]
            # cartesian gradient
            grad = self.data["Cartesian Gradient"]
            # Only the lower triangular part of the Hessian is stored.
            hess = np.zeros((3*nat,3*nat))
            row, col = np.tril_indices(3*nat)
            hess[row,col] = self.data["Cartesian Force Constants"]
            # Hessian is symmetric, H^T = H
            hess[col,row] = hess[row,col]
        except KeyError as err:
            logger.error(f"A required field could not be found in formatted checkpoint file {self.filename} .")
            raise err
        return pos, energy, grad, hess
    def nonadiabatic_coupling(self):
        """
        extract non-adiabatic coupling vector between ground and excited state (Root=I), if present.

        Only Gaussian 16 saves the NAC vector in the checkpoint file, while QChem writes it to the output file.

        Returns
        -------
        nac  :  ndarray (3*nat,)
          1st order derivative coupling <0|d/dx|I>
        """
        try:
            nac = self.data["Nonadiabatic coupling"]
        except KeyError as err:
            logger.error(f"The field `Nonadiabatic coupling` could not be found in the formatted checkpoint file {self.filename} .")
            raise err
        if (nac == 0.0).all():
            logger.warning(f"All components of non-adiabatic coupling vector in {self.filename} are zero.")
        return nac
    def vibrational_groundstate(self, zero_threshold=100.0):
        """
        The vibrational ground state belonging to the harmonic potential is given by

                                          1/4                 T
           psi (x) = (det(Gamma ) / pi^N)    exp{ -1/2 (x-x )  Gamma  (x-x ) }
              0                0                           0        0     0

        provided that x0 is the minimum. This function computes the width parameter matrix
        Gamma_0 from the Hessian at the minimum.

        Optional
        --------
        zero_threshold   :  float > 0
          threshold for considering normal mode frequencies as zero (in cm-1)

        Returns
        -------
        x0      : ndarray (3*nat,)
          center of Gaussian, in cartesian coordinates (bohr)
        Gamma0  : ndarray (3*nat,3*nat)
          symmetric, positive semi-definite matrix of width parameters (bohr^{-2})
        en_zpt  : float
          zero-point energy (Hartree)
        """
        x0, energy, grad, hess = self.harmonic_approximation()
        mass = self.masses()
        # diagonals of M^{1/2} and M^{-1/2}
        msq = np.sqrt(mass)
        imsq = 1.0/msq
        # mass-weighted Hessian H
        hess_mwc = np.einsum('i,ij,j->ij', imsq, hess, imsq)
        # diagonalize symmetric  H = V.diag(w).V^T
        w2,V = sla.eigh(hess_mwc)

        # vibrational energies
        w = np.sqrt(w2)
        # zero-point energy
        en_zpt = 0.5 * hbar * np.sum(w)
        
        logger.info("Normal mode frequencies (cm-1)")
        logger.info(w*units.hartree_to_wavenumbers)

        if not (w * units.hartree_to_wavenumbers > zero_threshold).all():
            logger.warning("At a minimum all frequencies should be positive, found imaginary ones.")
        
        # select non-zero vibrational modes
        non_zero = (w * units.hartree_to_wavenumbers) > zero_threshold
        # number of non singular dimensions
        num_non_zero = np.count_nonzero( non_zero )

        dim = x0.shape[0]
        logger.info(f"number of zero modes : {dim - num_non_zero}")

        # L = hbar^{-1/2} M^{1/2} V w^{1/2}
        L = hbar**(-1/2) * np.einsum('i,ij,j->ij', msq, V[:,non_zero], np.sqrt(w[non_zero]))

        # Gamma_0 = L . L^T
        Gamma_0 = np.einsum('ij,kj->ik', L, L)

        return x0, Gamma_0, en_zpt
        
    def masses(self):
        """
        atomic masses in a.u.

        Returns
        -------
        masses  :  ndarray (3*nat,)
          masses for each cartesian coordinate in multiples of electron mass
        """
        mass = self.data["Real atomic weights"] * units.amu_to_aumass
        mass = np.repeat(mass, 3)
        return mass
    def atomic_numbers(self):
        """
        atomic numbers 

        Returns
        -------
        numbers :  ndarray(nat,)
          atomic number for each atom
        """
        return self.data["Atomic numbers"]
