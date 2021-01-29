
Description
-----------
In order to train machine-learned potentials such as sGDML, one has to generate
large amounts of training data. Given a sequence of geometries, this
set of scripts distributes QChem jobs over a SLURM queue and stores
the forces and NAC vectors in the format expected by sGDML. 

Requirements
------------

 * QChem
 * SLURM queue
   
Installation
------------
.. code-block:: bash

   $ pip install -e .
   
Getting Started
---------------

In 'examples/formaldehyde', run

.. code-block:: bash

   $ sbatch --job-name='create dataset' << EOF
     #!/bin/bash
	forces_trajectory.py   geometries.xyz  qchem.in  --parallel_images=10
     EOF

After some time, this should produce the following files in extended XYZ format (using atomic units):

 * forces_0.xyz      -   forces in S0
 * forces_1.xyz      -   forces in S1
 * nacvec_0-1.xyz    -   NAC vectors between S0 and S1


Geometries at which the forces are calculated, can be sampled from the Maxwell-Boltzmann distribution
by running molecular dynamics at high temperature with the ANI1cxx force field.

.. code-block:: bash

   $ ani1cxx_dynamics.py  initial.xyz --temperature=500.0  -o geometries.xyz 

