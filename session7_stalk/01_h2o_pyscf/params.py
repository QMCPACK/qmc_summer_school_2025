#!/usr/bin/env python3

from numpy import array, sin, cos, ndarray, pi

from pyscf import dft
from pyscf import gto
from pyscf.geomopt.geometric_solver import optimize
from pyscf.gto.mole import tofile

from stalk.params.util import bond_angle, mean_distances, mean_param
from stalk import ParameterStructure
from stalk.params import PesFunction


# Natural forward mapping using bond lengths and angles
def forward_natural(pos: ndarray):
    pos = pos.reshape(-1, 3)  # make sure of the shape
    # for easier comprehension, list particular atoms
    O0 = pos[0]
    H0 = pos[1]
    H1 = pos[2]

    # for redundancy, calculate mean bond lengths
    r = mean_distances([
        (O0, H0),
        (O0, H1)
    ])
    a = bond_angle(H0, O0, H1, units='rad')
    params = [r, a]
    return params
# end def


# Backward mapping: produce array of atomic positions from bond length and angle
def backward_natural(params: ndarray):
    r, a = tuple(params)
    # Transform bond angle to triangular angle
    a = (pi - params[1]) / 2
    # place atoms on the xy-plane
    O0 = [0.0, 0.0, 0.0]
    H0 = [r * cos(a), r * sin(a), 0.0]
    H1 = [-r * cos(a), r * sin(a), 0.0]
    pos = array([O0, H0, H1])
    return pos
# end def


# Auxiliary forward mapping using Cartesian support variables
def forward(pos: ndarray):
    pos = pos.reshape(-1, 3)  # make sure of the shape
    # for easier comprehension, list particular atoms
    O0 = pos[0]
    H0 = pos[1]
    H1 = pos[2]

    # Assume that all atoms lie parallel to xy-plane
    # Assume that 2 H atoms lie at the same y with O in the middle

    # Parameter #0, the horizontal offset
    x = mean_param([
        (H0 - O0)[0],  # first H on the 'right'
        -(H1 - O0)[0]  # second H on the 'left'
    ])
    y = mean_param([
        (H0 - O0)[1],  # positive value when H's are 'above' O
        (H1 - O0)[1]
    ])
    params = array([x, y])
    return params
# end def


# Backward mapping: produce array of atomic positions from x, y offsets
def backward(params: ndarray):
    x, y = tuple(params)
    O0 = [0.0, 0.0, 0.0]
    H0 = [x, y, 0.0]  # first H on the 'right'
    H1 = [-x, y, 0.0]  # second H on the 'left'
    pos = array([O0, H0, H1])
    return pos
# end def


# Common definitions of the PySCF kernel
def kernel_pyscf(structure: ParameterStructure):
    atom = []
    for el, pos in zip(structure.elem, structure.pos):
        atom.append([el, tuple(pos)])
    # end for
    mol = gto.Mole()
    mol.atom = atom
    mol.verbose = 2
    mol.basis = 'ccpvdz'
    mol.unit = 'A'
    mol.ecp = 'ccecp'
    mol.charge = 0
    mol.spin = 0
    mol.symmetry = False
    mol.cart = True
    mol.build()

    mf = dft.RKS(mol)
    return mf
# end def


# Relaxation job takes in a structure and outputs (to file) the relaxed structure
def relax_pyscf(structure: ParameterStructure, outfile='relax.xyz', xc='pbe'):
    mf = kernel_pyscf(structure=structure)
    mf.xc = xc
    mf.kernel()
    mol_eq = optimize(mf, maxsteps=100)
    # Write to external file
    tofile(mol_eq, outfile, format='xyz')
# end def


# PES job takes in a structure and returns the energy and errorbar
def pes_pyscf(structure: ParameterStructure, xc='pbe', **kwargs):
    print(f'Computing: {structure.label} ({xc})')
    mf = kernel_pyscf(structure=structure)
    mf.xc = xc
    e_scf = mf.kernel()
    return e_scf, 0.0
# end def


# The PES functions can be wrapped and passed to line-searches with
# dynamic arguments using the PesFunction class
pes_pbe = PesFunction(pes_pyscf, {'xc': 'pbe'})
pes_b3lyp = PesFunction(pes_pyscf, {'xc': 'b3lyp'})
