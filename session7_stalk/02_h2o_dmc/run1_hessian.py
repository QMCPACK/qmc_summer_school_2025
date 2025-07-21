#!/usr/bin/env python

from stalk import ParameterHessian

from params import pes_pyscf
from run0_relax import structure_relax

interactive = __name__ == "__main__"

hessian = ParameterHessian(structure=structure_relax)
hessian.compute_fdiff(
    path='fdiff',
    pes=pes_pyscf,
    dp=0.001,  # Finite displacements along each parameter
    interactive=interactive,
)

if interactive:
    print('Hessian:')
    print(hessian)
# end if
