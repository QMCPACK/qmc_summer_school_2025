#!/usr/bin/env python3

# The following script generates and plots selected figures appearing in the presentation.

import numpy as np
from matplotlib import pyplot as plt

from stalk import ParameterHessian, ParameterSet, TargetLineSearch
from stalk.ls import LineSearchBase
from stalk.params import PesFunction


# Define a simple morse potential
def morse(structure: ParameterSet, p=None, **kwargs):
    r = structure.params[0]
    return p[2] * ((1 - np.exp(-(r - p[0]) / p[1]))**2 - 1) + p[3], 0.0
# end def


# Request a simple line-search based on interpolated data
def get_ls(tls: TargetLineSearch, R, M=5, fit_kind='pf2', sigma=None):
    offsets = tls.figure_out_adjusted_offsets(R=R, M=M)
    values = tls.evaluate_target(offsets)
    if sigma is not None:
        errors = len(values) * [sigma]
    else:
        errors = None
    # end if
    ls = LineSearchBase(fit_kind=fit_kind, offsets=offsets, values=values, errors=errors)
    return ls
# end def


# well stiffness
h = 1.0
# eqm distance
r0 = 2.0
# Morse parameters (eqm value, stiffness, well depth, E_inf)
p0 = [r0, h, 0.5, 0.0]
pes = PesFunction(morse, {'p': p0})
s = ParameterSet([r0])
hessian = ParameterHessian(structure=s, hessian=[[h]])

tls = TargetLineSearch(
    structure=s,
    pes=pes,
    hessian=hessian,
    fit_kind='pf2',
    d=0,
    M=21,
    R=1.5,
)

R_full = tls.figure_out_adjusted_offsets(R=tls.R_max, M=101)
E_full = tls.evaluate_target(R_full)

# Plot a few energy curves
ls05 = get_ls(tls, R=0.5)
ls07 = get_ls(tls, R=0.7)
f, ax = plt.subplots()
ax.plot(R_full, E_full, 'k')
ls05.plot(ax=ax, color='tab:blue', marker='<')
ls07.plot(ax=ax, color='tab:orange', marker='>')
ax.set_ylim([-0.6, 0.2])
plt.show()

# Plot bias vs R
f, ax = plt.subplots()
Rs = np.linspace(0.0, 1.5, 51)
Rs2, Ws2, biases_pf2 = tls.compute_bias_of(R=Rs, M=7, fit_kind='pf2')
Rs3, Ws3, biases_pf3 = tls.compute_bias_of(R=Rs, M=7, fit_kind='pf3')
Rs4, Ws4, biases_pf4 = tls.compute_bias_of(R=Rs, M=7, fit_kind='pf4')
ax.set_title('Fitting bias vs R')
ax.set_xlabel('Grid extent R')
ax.set_ylabel('Fitting bias')
ax.plot(Ws2, biases_pf2, label='n=2')
ax.plot(Ws3, biases_pf3, label='n=3')
ax.plot(Ws4, biases_pf4, label='n=4')
plt.legend()
plt.show()

# Plot collection of PES fits
f, ax0 = plt.subplots()
f, ax1 = plt.subplots()
sigma = 0.02
Gs = np.random.randn(500, 5) * sigma
xgrid = np.linspace(-0.6, 0.6, 101)
x0s = []
for G in Gs:
    ls3 = get_ls(tls, R=0.5, fit_kind='pf3', M=5)
    ls3.values += G
    res = ls3.search()
    x0s.append(res.x0)
    ygrid = np.polyval(res.fit, xgrid)
    ax0.plot(xgrid, ygrid, 'k', alpha=0.05)
# end for
ls3 = get_ls(tls, R=0.5, fit_kind='pf3', M=5)
ax0.errorbar(ls3.offsets, ls3.values, 5 * [2 * sigma])
ax0.set_xlabel('R')
ax1.set_ylabel('E')
ax1.hist(abs(np.array(x0s - ls3.x0)), bins=20)
ax1.set_xlabel('Fitting error')
ax1.set_ylabel('Distribution')
plt.show()


# Plot noise vs R
f, ax = plt.subplots()
Rs_noise = np.linspace(0.4, 1.0, 12)
errors2 = []
errors3 = []
errors4 = []
sigma = 0.005
Gs = np.random.randn(500, 5)
for R in Rs_noise:
    ls2 = get_ls(tls, R=R, fit_kind='pf2', M=5, sigma=sigma)
    ls3 = get_ls(tls, R=R, fit_kind='pf3', M=5, sigma=sigma)
    ls4 = get_ls(tls, R=R, fit_kind='pf4', M=5, sigma=sigma)
    errors2.append(tls.compute_errorbar(grid=ls2, fit_kind='pf2', Gs=Gs)[0])
    errors3.append(tls.compute_errorbar(grid=ls3, fit_kind='pf3', Gs=Gs)[0])
    errors4.append(tls.compute_errorbar(grid=ls4, fit_kind='pf4', Gs=Gs)[0])
# end for
ax.plot(Rs_noise, errors2, label='n=2')
ax.plot(Rs_noise, errors3, label='n=3')
ax.plot(Rs_noise, errors4, label='n=4')
ax.set_title('Fitting uncertainty vs R')
ax.set_ylabel('Fitting uncertainty')
ax.set_xlabel('Grid extent R')
plt.legend()
plt.show()


# Plot total error surface
f, ax = plt.subplots()
tls.optimize(
    fit_kind='pf3',
    M=5,
    N=500,
    epsilon=0.04
)
tls.plot_error_surface(ax=ax)
f, ax = plt.subplots()
tls.optimize(
    epsilon=0.08,
)
tls.plot_error_surface(ax=ax)
plt.show()

# plot statistical cost, W_opt, sigma_opt vs tolerance
f, (ax0, ax1, ax2) = plt.subplots(3, 1, sharex=True)
f.subplots_adjust(hspace=0.02)
costs = []
sigmas = []
windows = []
epsilons = np.linspace(0.01, 0.1, 10)
for epsilon in epsilons:
    tls.optimize(epsilon=epsilon)
    costs.append(tls.statistical_cost())
    sigmas.append(tls.sigma_opt)
    windows.append(tls.W_opt)
# end for
ax0.plot(epsilons, costs)
ax1.plot(epsilons, sigmas)
ax2.plot(epsilons, windows)
ax0.set_ylabel('Stat. cost')
ax1.set_ylabel(r'$W_opt$')
ax2.set_ylabel(r'$sigma_opt$')
ax2.set_xlabel('Error tolerance')
ax0.set_yticks([])
ax1.set_yticks([])
ax2.set_yticks([])
plt.show()
