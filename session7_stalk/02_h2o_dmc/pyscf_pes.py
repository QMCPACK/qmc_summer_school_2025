#!/usr/bin/env python3

# This is a template files used by Nexus to generate PySCF PES function

from pyscf import dft
from numpy import savetxt

$system

### generated calculation text ###
mf = dft.RKS(mol)
mf.xc = 'pbe'
e_scf = mf.kernel()
### end generated calculation text ###

savetxt('energy.dat', [[e_scf, 0.0]])