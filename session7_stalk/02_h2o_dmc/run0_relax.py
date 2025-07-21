#!/usr/bin/env python3

from numpy import array, pi

from stalk.nexus import NexusStructure

from params import forward, backward, backward_natural, relax_pyscf


interactive = __name__ == "__main__"

# Let us initiate a ParameterStructure object that implements the parametric mappings
params_init_natural = array([0.97, 104.0 / 180 * pi])
pos_init = backward_natural(params_init_natural)
elem = ['O'] + 2 * ['H']
structure_init = NexusStructure(
    forward=forward,
    backward=backward,
    pos=pos_init,
    elem=elem,
    units='A',
    label='init'
)

structure_relax = structure_init.copy()

relax_pyscf.relax(
    structure_relax,
    path='relax/',
    interactive=interactive,
)

if interactive:
    print('Initial params:')
    print(structure_init.params)
    print('Relaxed params:')
    print(structure_relax.params)
# end if
