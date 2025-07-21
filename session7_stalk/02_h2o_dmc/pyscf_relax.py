#!/usr/bin/env python3

# This is a template files used by Nexus to generate PySCF relaxation function

from pyscf import dft

$system

### generated calculation text ###
mf = dft.RKS(mol)
mf.xc = 'pbe'
e_scf = mf.kernel()
### end generated calculation text ###

from pyscf.geomopt.geometric_solver import optimize
mol_eq = optimize(mf, maxsteps=100)
from pyscf.gto.mole import tofile
# Write to external file
tofile(mol_eq, 'relax.xyz', format='xyz')
# Write to output file 
print(mol_eq.atom_coords())
