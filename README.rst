
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
   
.. code-block:: bash

   $ forces_trajectory.py   geometries.xyz  qchem.in  --parallel_images=100


   
