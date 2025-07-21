#!/usr/bin/env python3

from matplotlib import pyplot as plt

from stalk import LineSearchIteration

from run0_relax_b3lyp import structure_relax as structure_b3lyp
from params import pes_b3lyp
from run2_surrogate import surrogate


# Copy the surrogate structure and shift to make for more challenge
shifted_structure = surrogate.structure.copy()
shifted_structure.shift_params([-0.2, 0.2])
# Then generate line-search iteration object based on the shifted surrogate
srg_ls = LineSearchIteration(
    surrogate=surrogate,
    structure=shifted_structure,
    path='b3lyp_ls',
    pes=pes_b3lyp,
)
# Iterate 4 times
for i in range(4):
    srg_ls.propagate(i, add_sigma=True)
# end for
# Evaluate the latest eqm structure
srg_ls.pls().evaluate_eqm(add_sigma=True)

if __name__ == '__main__':
    print(srg_ls)
    print('True (B3LYP) energy and params:')
    print(structure_b3lyp.value, structure_b3lyp.params)
    srg_ls.plot_convergence(targets=structure_b3lyp.params)
    plt.show()
# end if
