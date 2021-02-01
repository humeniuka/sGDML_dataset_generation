#!/usr/bin/env python
"""
run molecular dynamics with ANI-2x at constant temperature to generate a Maxwell-Boltzmann distributed
ensemble of geometries
"""
import argparse
import sys
import numpy as np
import logging

import tqdm

import torch
import torchani

import ase
import ase.io
from ase.md.nvtberendsen import NVTBerendsen
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution, Stationary, ZeroRotation

# # Logging
logger = logging.getLogger(__name__)
logging.basicConfig(format="[%(module)-12s] %(message)s", level='DEBUG')

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument(
    '--version',
    action='version',
    version='%(prog)s '
        + ' [Python {}, NumPy {}'.format(
            '.'.join(map(str, sys.version_info[:3])), np.__version__
        )
        + ', PyTorch {}'.format(torch.__version__)
        + ', ASE {}'.format(ase.__version__)
        + ']',
)
parser.add_argument('geometry_input',
                    type=str, metavar='initial.xyz', help='molecular geometry in any format supported by ASE')
parser.add_argument('-o', '--output',
                    dest='samples_output',
                    type=str,
                    metavar='geometries.xyz',
                    default='geometries.xyz',
                    help='samples along ground state trajectory are saved to this file.')
parser.add_argument(
    '--cuda',
    type=int, dest='cuda', default=0, metavar='id', 
    help="select id of cuda device if more than one is available, i.e. 0 for 'cuda:0'")

parser.add_argument("--save_every",
                    dest="save_every",
                    type=int,
                    metavar='N',
                    default=100,
                    help="save every N-th dynamics step")
parser.add_argument("--samples",
                    dest="samples",
                    type=int,
                    default=10000,
                    help="number of samples, the total number of time steps equals SAMPLES x SAVE_EVERY")
parser.add_argument("--dt_fs",
                    dest="dt_fs",
                    type=float,
                    default=1.0,
                    help="duration of a single time step (in fs)")
parser.add_argument("--temperature",
                    dest="temperature",
                    type=float,
                    default=300.0,
                    help="temperature in Kelvin")

args = parser.parse_args()

logger.info(f"arguments: {args}")

torch.set_default_dtype(torch.float64)
if torch.cuda.is_available():
    device = torch.device(f"cuda:{args.cuda}")
else:
    device = torch.device('cpu')
logger.info(f"running on device {device}")
    
atoms = ase.io.read(args.geometry_input)

assert set(atoms.get_chemical_symbols()) - set(['H', 'C', 'N', 'O', 'F', 'S', 'Cl']) == set(), "ANI-2x is only parametrized for H,C,N,O,F,S,Cl"

model = torchani.models.ANI2x(periodic_table_index=True).double().to(device)
atoms.calc = torchani.ase.Calculator(atoms.get_chemical_symbols(), model)

# Set the momenta corresponding to given temperature
MaxwellBoltzmannDistribution(atoms, temperature_K=args.temperature)
# eliminate translation and rotation
Stationary(atoms)
ZeroRotation(atoms)


logger.info(f"initial temperature T= {atoms.get_temperature():.1f} K")

dyn = NVTBerendsen(atoms, args.dt_fs * ase.units.fs,
                   args.temperature, 
                   taut=0.5*1000*ase.units.fs)

# save initial geometry
ase.io.extxyz.write_extxyz(args.samples_output, [atoms],
                           columns=['symbols', 'positions'],
                           write_results=False, append=False)

logger.info(f"starting dynamics, target temperature= {args.temperature} K")
with tqdm.tqdm(total=args.samples) as progress_bar:
    for i in range(0, args.samples):
        # run for several steps, then save a sample
        dyn.run(args.save_every)
        
        ase.io.extxyz.write_extxyz(args.samples_output, [atoms],
                                   columns=['symbols', 'positions'],
                                   write_results=False, append=True)

        progress_bar.set_description(f" T= {atoms.get_temperature():.1f} K   time/fs= {(i+1)*args.dt_fs*args.save_every}")
        progress_bar.update(1)

    
