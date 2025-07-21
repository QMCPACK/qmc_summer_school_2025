#!/usr/bin/env python

from stalk import TargetParallelLineSearch
from matplotlib import pyplot as plt

from params import pes_pbe
from run1_hessian import hessian


surrogate_file = 'surrogate.p'
surrogate = TargetParallelLineSearch(
    path='surrogate/',  # Choose path for storing auxiliary files
    fit_kind='pf3',  # cubic fit
    load=surrogate_file,  # Try to load from disk to avoid recomputation
    structure=hessian.structure,  # Use this eqm structure
    hessian=hessian,  # Use these directions
    pes=pes_pbe,  # Use this PES
    window_frac=0.25,  # Set the initial displacements: W = sqrt(Lambda) * window_frac
    M=15  # Use this many points per direction
)
# Eliminate relaxation deficiencies by bracketing true minimum
surrogate.bracket_target_biases()

epsilon_p = [0.02, 0.02]
surrogate.optimize(
    epsilon_p=epsilon_p,  # Set parameter tolerances
    fit_kind='pf3',  # Optimize for cubic fit
    M=7,  # Grid size
    N=400,  # This much resampling of the errorbars
    reoptimize=False,  # Do not reoptimize on import
    write=surrogate_file,  # Write to disk afterwards
)

if __name__ == '__main__':
    print('W_opt: ', surrogate.W_opt)
    print('sigma_opt:', surrogate.sigma_opt)
    surrogate.plot()
    surrogate.plot_error_surfaces()
    plt.show()
# end if
