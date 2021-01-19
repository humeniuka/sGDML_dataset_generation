#!/usr/bin/env python
import numpy as np
import argparse
import itertools
import sys
import os
import os.path
from multiprocessing import Pool
import tempfile
import shutil
import ase.io
import ase.atoms

from sgdml_dataset_generation.calculators import get_calculator

# This is a wrapper function for making run_calculator compatible with Pool.map(...)
# The run_calculator() function has to be defined at the top level.
def run_calculator_map(args):
    atoms, igeo, kwds = args
    # create a temporary directory
    directory = "TMP/GEOM_%5.5d" % igeo
    if not os.path.exists(directory):
        os.makedirs(directory)
    res = run_calculator(atoms, directory=directory, igeo=igeo, **kwds)
    shutil.rmtree(directory)
    return res
    
if __name__ == "__main__":
    
    usage = """
       Usage: forces_trajectory.py  geometries.xyz  script.in

    Compute forces and coupling vectors for all geometries in the xyz-file 
    and save them together with the energies and geometries in the extended xyz-format. 

    Input Files:
       geometries.xyz   -   xyz-file with geometries, the atom order has to be the same 
                            for all geometries
       script.in        -   

    Output Files:
       forces_<I>.xyz       -   extended xyz-file with geometries, energies and forces for
                                each geometry.
       nacvec_<I>-<J>.xyz   -   extended xyz-file with geometries and non-adiabatic coupling
                                vectors for each geometry.

    Separate output files are created for each gradient (of state I) and coupling vector
    (between states I and J) that are found in the output file.
    """
    
    parser = argparse.ArgumentParser(prog="forces_trajectory.py", usage=usage)

    parser.add_argument("geometries_file", type=str, default="geometries.xyz",
                      help="XYZ file containing geometries [default: geometries.xyz]")
    parser.add_argument("script_file", type=str, default="grad.in",
                        help="QChem script controlling the calculation of gradients and couplings.")
    

    # options
    parser.add_argument("-i", "--procs_per_image", dest="procs_per_image", type=int, default=1,
                      help="Number of processors used in the calculation for a single image [default: 1]")
    parser.add_argument("-m", "--mem_per_image", dest="mem_per_image", type=str, default="6Gb",
                      help="Amount of memory to be allocated for each image calculation [default: 6Gb]")
    parser.add_argument("-p", "--parallel_images", dest="parallel_images", type=int, default=1,
                      help="How many images should be processed in parallel? Each image will be run with `procs_per_image` processors [default: 1]")
    parser.add_argument("-s", "--scratch_dir", dest="scratch_dir", type=str, default="/sscratch/${SLURM_JOBID}",
                      help="Path to scratch directory [default: /sscratch/${SLURM_JOBID}]")
    parser.add_argument("--calculator", dest="calculator", type=str, default="qchem",
                        help="Choose electronic structure program [default: qchem]")
    
    args = parser.parse_args()

    """
    try:
        os.environ["SLURM_JOBID"]
    except KeyError:
        raise KeyError("This script can only be run through a SLURM queue, environment variable SLURM_JOBID not found!")
    """

    if args.calculator == "qchem":
        if not os.path.exists(args.script_file):
            raise RuntimeError(f"ERROR: QChem input script '{args.script_file}' not found in current folder!")
        ext = os.path.splitext(args.script_file)[1]
        if ext != ".in":
            raise RuntimeError(f"ERROR: File of extension of QChem input script should be '.in' but got '{ext}'")

    molecules = ase.io.iread(args.geometries_file)
    
    global run_calculator
    run_calculator = get_calculator(args.calculator)

    pool = Pool(args.parallel_images)

    # additional keywords controlling electronic structure calculation
    kwds = {
        "script" : args.script_file,
        "nprocs" : args.procs_per_image,
        "mem": args.mem_per_image
    }

    # With python2 imap(...) returns a list of all results after the last
    # calculation has finished, with python3 the results are returned as
    # soon as they are available.
    results_list = pool.imap(run_calculator_map,
                        zip(
                            molecules,
                            itertools.count(1),
                            itertools.repeat(kwds)) )

    # In the first step the force file is overwritten if it exists
    # already, later mode is set to append.
    append = False

    # Save results to file as they keep coming in.
    for igeom,results in enumerate(results_list):
        mol = ase.atoms.Atoms(symbols=results["symbols"])
        mol.set_positions(results["geometry (au)"])
        mol.info["Units"] = "a.u."

        # a separate file for each gradient
        for state,gradient in results["gradients (au)"].items():
            #
            mol.info["Energy"] = results["energies (au)"][state]
            # forces are stored as 'momenta'
            mol.set_momenta( - gradient)

            ase.io.extxyz.write_extxyz(f"forces_{state}.xyz", [mol],
                                       columns=['symbols', 'positions', 'momenta'],
                                       write_results=False, append=append)

        # separate file for each NAC vector
        for (state1,state2),nac_vector in results["derivative couplings (au)"].items():
            #
            mol.info["Energy"] = results["energies (au)"][state]
            # NAC vectors are stored as 'momenta'
            mol.set_momenta(nac_vector)

            ase.io.extxyz.write_extxyz(f"nacvec_{state1}-{state2}.xyz", [mol],
                                       columns=['symbols', 'positions', 'momenta'],
                                       write_results=False, append=append)
            
            
        append = True

    print("FINISHED")
    
