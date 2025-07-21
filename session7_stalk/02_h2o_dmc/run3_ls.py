#!/usr/bin/env python3

from matplotlib import pyplot as plt

from stalk import LineSearchIteration

from params import pes_dmc
from run2_surrogate import surrogate


interactive = __name__ == "__main__"

# Copy the surrogate eqm structure and characterize QMC error
structure_qmc = surrogate.structure.copy()
var_eff = pes_dmc.get_var_eff(
    structure_qmc,
    path='dmc_var_eff',
    samples=10,
    interactive=interactive
)
pes_dmc.args['var_eff'] = var_eff

# Then generate line-search iteration object based on the shifted surrogate
structure_qmc.shift_params([-0.2, 0.2])
dmc_ls = LineSearchIteration(
    surrogate=surrogate,
    structure=structure_qmc,
    path='dmc_ls',
    pes=pes_dmc,
)
# Iterate 3 times
for i in range(3):
    dmc_ls.propagate(i, interactive=interactive)
# end for
# Evaluate the latest eqm structure
dmc_ls.pls().evaluate_eqm()

if interactive:
    print(dmc_ls)
    dmc_ls.plot_convergence()
    plt.show()
# end if
