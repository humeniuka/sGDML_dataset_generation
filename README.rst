
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

   $ sbatch --jobname='create dataset' << EOF
	forces_trajectory.py   geometries.xyz  qchem.in  --parallel_images=100
     EOF

After some time, this should produce the following files in extended XYZ format:

 * forces_0.xyz      -   forces in S0
 * forces_1.xyz      -   forces in S1
 * nacvec_0-1.xyz    -   NAC vectors between S0 and S1


   
