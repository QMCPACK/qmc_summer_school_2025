#!/usr/bin/env python3

import numpy as np
from matplotlib import pyplot as plt

from stalk import LineSearchIteration, ParameterHessian, ParameterSet
from stalk.params import PesFunction


def hessian_pes(structure: ParameterSet, H=None, **kwargs):
    p = structure.params
    if H is None:
        # If Hessian is not provided, revert to identity
        H = np.diag(len(p) * [1])
    # end if
    return 0.5 * p.T @ H @ p, 0.0
# end def


H0 = None
H1 = np.array([
    [0.7, -1.0],
    [-1.0, 2.0]
])

p_init = [2, 2.5]
pes = PesFunction(hessian_pes, {'H': H1})

# Start from same position
p0 = ParameterSet(p_init)
p1 = ParameterSet(p_init)
# Hessian
h0 = ParameterHessian(structure=p0, hessian=H0)
h1 = ParameterHessian(structure=p1, hessian=H1)
# Parallel line-search iterations
lsi0 = LineSearchIteration(path='lsi0', structure=p0, hessian=h0, pes=pes, fit_kind='pf2', R=3.0, M=7)
lsi1 = LineSearchIteration(path='lsi1', structure=p1, hessian=h1, pes=pes, fit_kind='pf2', R=3.0, M=7)

# Iterate both line-searches
for i in range(10):
    lsi0.propagate(i=i)
    lsi1.propagate(i=i)
# end for

print(lsi0)
print(lsi1)


# Construct the PES on a grid
xs = np.linspace(-4, 4, 51)
ys = np.linspace(-4, 4, 51)
try:
    # Try to save and load the data to avoid recomputations.
    pes_data = np.loadtxt('pes.dat')
except FileNotFoundError:
    pes_data = []
    for x in xs:
        pes_row = []
        for y in ys:
            e = hessian_pes(ParameterSet([x, y]), H=H1)[0]
            pes_row.append(e)
        # end for
        pes_data.append(pes_row)
    # end for
    pes_data = np.array(pes_data)
    np.savetxt('pes.dat', pes_data)
# end try


X, Y = np.meshgrid(xs, ys)
Lambda, U = np.linalg.eig(H1)

f, ax = plt.subplots(figsize=(4, 4))

levels = np.linspace(0.0, np.max(pes_data)**0.5, 21)**2
ax.contourf(X, Y, pes_data, levels=levels, alpha=0.2, cmap='cividis')
ax.contour(X, Y, pes_data, levels=levels, alpha=0.5)
ax.plot(*p_init, 'kd')

s = 2.0
ax.arrow(0, 0, s, 0, color='k', head_width=0.2, linestyle=':')
ax.arrow(0, 0, 0, s, color='k', head_width=0.2, linestyle=':', label='Parameter directions')

ax.arrow(0, 0, s * U.T[0, 0], s * U.T[1, 0], color='tab:red', head_width=0.2, linestyle=':')
ax.arrow(0, 0, s * U.T[0, 1], s * U.T[1, 1], color='tab:red', head_width=0.2, linestyle=':', label='Conjugate directions')

for i in range(len(lsi0) - 1):
    p0 = lsi0.pls(i).structure.params
    p0_next = lsi0.pls(i).structure_next.params
    p1 = lsi1.pls(i).structure.params
    p1_next = lsi1.pls(i).structure_next.params
    ax.arrow(*p0, *(p0_next - p0), color='k', linestyle='-', head_width=0.1)
    ax.arrow(*p1, *(p1_next - p1), color='tab:red', linestyle='-', head_width=0.1)
# end for

ax.set_xticks([])
ax.set_yticks([])

plt.legend()
plt.show()
