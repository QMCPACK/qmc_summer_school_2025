#!/usr/bin/env python3

from numpy import array, pi

from stalk import ParameterStructure
from stalk.io import XyzGeometry

from params import forward, backward, backward_natural, relax_pyscf, pes_b3lyp


# Let us initiate a ParameterStructure object that implements the parametric mappings
# There are other ways to do this, but here we happen to have Wikipedia and the
# natural mappings at hand, so why not use them.
params_init_natural = array([0.97, 104.0 / 180 * pi])
pos_init = backward_natural(params_init_natural)
elem = ['O'] + 2 * ['H']
structure_init = ParameterStructure(
    forward=forward,
    backward=backward,
    pos=pos_init,
    elem=elem,
    units='A',
    label='init'
)

# This looks a bit hacky, because the code usually runs on Nexus and the intermediate
# files have not been implemented for non-Nexus workflows
outfile = 'relax_b3lyp.xyz'
try:
    # If the file is found, load it
    geom = XyzGeometry({'suffix': outfile}).load('./')
except FileNotFoundError:
    # Else, use the relaxation run to produce it
    relax_pyscf(structure_init, outfile, xc='b3lyp')
    geom = XyzGeometry({'suffix': outfile}).load('./')
# end try
# Obtain the new parameters by mapping forward the XYZ coordinates
new_params = structure_init.map_forward(geom.get_pos())
structure_relax = structure_init.copy(params=new_params, label='eqm')
# Evaluate the PES at eqm to get the reference energy
pes_b3lyp.evaluate(structure_relax)

if __name__ == '__main__':
    print(structure_init)
    print(structure_relax)
# end if
