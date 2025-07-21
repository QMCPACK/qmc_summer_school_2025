# Session 7: STALK examples

The following couple of examples show how STALK can be used to analyze
line-searches and relax a simple H2O geometry using DMC. 

For some more tutorial discussion, see:
https://github.com/QMCPACK/stalk/tree/master/examples/nexus

## Morse line-search analysis

For convenience, `00_morse_analysis/` contains a file for regeneration and
reproduction of selected figures appearing in the presentation. The script is
not thoroughly documented but serves as an example for those determined to
understand under-the-hood operations of the STALK code.

```
    cd 00_morse_analysis/
    python3 plot_figures.py
```

Feel free to try out different parameters for the data and plots.

## Example 1: H2O relaxation in native PySCF

The first complete example allows to study the parallel line-search workflow in 

```
    cd 01_h2o_pyscf/
    python3 run0_relax.py
    python3 run1_hessian.py
    python3 run2_surrogate.py
    python3 run3_ls_b3lyp.py
```

If you're feeling lucky, invoke push-button mode by running just
`run3_ls_b3lyp.py`. This is effectively a push-button mode where dependencies
are cascaded through imports and diagnostic printouts and plots are omitted.


### Exercises

1. Repeat runs 0 and 1 with the natural parameters, `forward_natural` and `backward_natural`. How does the Hessian compare to that of the original setup?
    * Hint: Make copies of `run0_relax.py` and `run1_hessian.py` and make necessary edits.
1. Set up yet another SCF PES (e.g. LDA) and perform line-search iteration with it, using PBE as a surrogate.
    * Hint: Make a copy of `run3_ls_b3lyp.py` and edit necessary parts there and in `params.py`. Shifting the structure and comparison to "exact" reference are not mandatory.
1. Reoptimize the surrogate model and observe statistical cost
    * In python console: `from run2_surrogate import surrogate`
    * Run `surrogate.optimize(epsilon_p=[x, x], M=y)` where x is a chosen tolerance and 3 < y < 11 is the number of grid points
    * Observe the estimated statistical cost with `surrogate.statistical_cost`.
    * Observe the target errorbars with `surrogate.sigma_opt`.
    * How does the statistical cost depend on the number of points? Or does it?

## Example 2: H2O relaxation with DMC

The second example is the same as first, but using Nexus to manage jobs and DMC
with QMCPACK to provide the stochastic PES. To run all steps:

```
    cd 02_h2o_dmc/
    python3 run0_relax.py
    python3 run1_hessian.py
    python3 run2_surrogate.py
    python3 run3_ls.py
```

To suppress plotting and interactive Nexus prompts, set `interactive=False` in
the scripts. One full iterations takes between 10-30 minutes on the VM, as it
involves a few dozen DMC workflows. 

Note that the error estimation with DMC requires a test run that characterizes
the effective DMC variance. This is done in the beginning of `run3_ls.py`.

### Exercises

1. After 3 iterations, what is the DMC eqm structure?
    * Hint: Import `from run3_ls import dmc_ls` and collect `s = dmc_ls.pls(i).structure_next` for every i, then consider the statistical aggregates of values `s.params` and their errorbars `s.params_err`.
    * Note: This analysis feature coming up in next release.
1. Does the structure meet the tolerances? If not, why could it be?
1. Is the assumption successful: DMC variance/errorbar remains constant in all displaced positions?
    * After line-search, load the DMC data in Python: `from run3_ls import dmc_ls`, then `print(dmc_ls.pls(0).ls(0).errors)`.
    * Compare to the target error `print(dmc_ls.pls(0).ls(0).sigma)`
    * What is the consequence for the line-search optimization?