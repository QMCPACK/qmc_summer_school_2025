#!/usr/bin/env python3

from numpy import array, sin, cos, ndarray, pi

from nexus import obj, job, settings, generate_pyscf, Structure
from nexus import generate_physical_system, generate_qmcpack, generate_convert4qmc

from stalk.util import EffectiveVariance
from stalk.nexus import NexusGeometry, NexusPes, QmcPes
from stalk.params.util import bond_angle, mean_distances, mean_param
from stalk.io import XyzGeometry, FilesLoader


qmcpseudos = ['O.ccECP.xml', 'H.ccECP.xml']
cores = 4
nx_settings = obj(
    sleep=3,
    pseudo_dir='./pseudos',
    runs='',
    results='',
    status_only=0,
    generate_only=0,
    machine=f'ws{cores}',
)

# Make sure to init Nexus only once
if len(settings) == 0:
    settings(**nx_settings)
# end if

# Configure
presub = ''
qmcapp = 'qmcpack'
pyscfjob = obj(app='python3', serial=True)
optjob = obj(app=qmcapp, cores=cores, ppn=cores, presub=presub)
dmcjob = obj(app=qmcapp, cores=cores, ppn=cores, presub=presub)


# Natural forward mapping using bond lengths and angles
def forward_natural(pos: ndarray):
    pos = pos.reshape(-1, 3)  # make sure of the shape
    # for easier comprehension, list particular atoms
    O0 = pos[0]
    H0 = pos[1]
    H1 = pos[2]

    # for redundancy, calculate mean bond lengths
    r = mean_distances([
        (O0, H0),
        (O0, H1)
    ])
    a = bond_angle(H0, O0, H1, units='rad')
    params = [r, a]
    return params
# end def


# Backward mapping: produce array of atomic positions from bond length and angle
def backward_natural(params: ndarray):
    r, a = tuple(params)
    # Transform bond angle to triangular angle
    a = (pi - params[1]) / 2
    # place atoms on the xy-plane
    O0 = [0.0, 0.0, 0.0]
    H0 = [r * cos(a), r * sin(a), 0.0]
    H1 = [-r * cos(a), r * sin(a), 0.0]
    pos = array([O0, H0, H1])
    return pos
# end def


# Auxiliary forward mapping using Cartesian support variables
def forward(pos: ndarray):
    pos = pos.reshape(-1, 3)  # make sure of the shape
    # for easier comprehension, list particular atoms
    O0 = pos[0]
    H0 = pos[1]
    H1 = pos[2]

    # Assume that all atoms lie parallel to xy-plane
    # Assume that 2 H atoms lie at the same y with O in the middle

    # Parameter #0, the horizontal offset
    x = mean_param([
        (H0 - O0)[0],  # first H on the 'right'
        -(H1 - O0)[0]  # second H on the 'left'
    ])
    y = mean_param([
        (H0 - O0)[1],  # positive value when H's are 'above' O
        (H1 - O0)[1]
    ])
    params = array([x, y])
    return params
# end def


# Backward mapping: produce array of atomic positions from x, y offsets
def backward(params: ndarray):
    x, y = tuple(params)
    O0 = [0.0, 0.0, 0.0]
    H0 = [x, y, 0.0]  # first H on the 'right'
    H1 = [-x, y, 0.0]  # second H on the 'left'
    pos = array([O0, H0, H1])
    return pos
# end def


# Define common SCF Mole arguments to keep consistent between relaxation and PES.
scf_mole_args = obj(
    spin=0,
    verbose=2,
    ecp='ccecp',
    basis='ccpvdz',
    symmetry=False,
)


# Nexus generator for SCF relaxation workflow
def scf_relax_job(structure: Structure, path, **kwargs):
    system = generate_physical_system(
        structure=structure,
        O=6,
        H=1,
    )
    relax = generate_pyscf(
        template='./pyscf_relax.py',
        system=system,
        identifier='relax',
        job=job(**pyscfjob),
        path=path,
        mole=scf_mole_args
    )
    return [relax]
# end def


# Nexus generator for SCF PES workflow
def scf_pes_job(structure: Structure, path, **kwargs):
    system = generate_physical_system(
        structure=structure,
        O=6,
        H=1,
    )
    scf = generate_pyscf(
        template='./pyscf_pes.py',
        system=system,
        identifier='scf',
        job=job(**pyscfjob),
        path=path,
        mole=scf_mole_args
    )
    return [scf]
# end def


# Nexus generator for DMC PES workflow
def dmc_pes_job(
    structure: Structure,
    path,
    sigma=None,
    samples=10,
    var_eff=None,
    **kwargs
):
    # Estimate the relative number of samples needed
    if isinstance(var_eff, EffectiveVariance):
        dmcsteps = var_eff.get_samples(sigma)
    else:
        dmcsteps = samples
    # end if

    system = generate_physical_system(
        structure=structure,
        O=6,
        H=1,
    )
    # Generate orbitals with PySCF
    scf = generate_pyscf(
        template='./pyscf_pes.py',
        system=system,
        identifier='scf',
        job=job(**pyscfjob),
        path=path + 'scf',
        mole=obj(
            verbose=4,
            ecp='ccecp',
            basis='ccecp-ccpvqz',  # Use larger basis to promote QMC performance
            symmetry=False,
        ),
        # Save orbitals for QMC
        save_qmc=True,
    )
    c4q = generate_convert4qmc(
        identifier='c4q',
        path=path + 'scf',
        job=job(cores=1),
        dependencies=(scf, 'orbitals'),
    )
    # Optimize J1 + J2
    opt = generate_qmcpack(
        system=system,
        path=path + 'opt',
        job=job(**optjob),
        dependencies=[(c4q, 'orbitals')],
        cycles=6,
        identifier='opt',
        qmc='opt',
        input_type='basic',
        pseudos=qmcpseudos,
        J2=True,
        J1_size=6,
        J1_rcut=6.0,
        J2_size=8,
        J2_rcut=8.0,
        minmethod='oneshift',
        blocks=200,
        substeps=2,
        samples=100000,
        minwalkers=0.3,
    )
    # Generate DMC job, where the number of steps per block is set
    # dynamically to meet the requested errorbar
    dmc = generate_qmcpack(
        system=system,
        path=path + 'dmc',
        job=job(**dmcjob),
        dependencies=[(c4q, 'orbitals'), (opt, 'jastrow')],
        steps=dmcsteps,
        identifier='dmc',
        qmc='dmc',
        input_type='basic',
        pseudos=qmcpseudos,
        jastrows=[],
        walkers_per_rank=128,
        blocks=200,
        timestep=0.01,
        ntimesteps=1,
    )
    # Store the relative samples for printout
    dmc.samples = dmcsteps
    return [scf, c4q, opt, dmc]
# end def


# Finally, wrap the Nexus job generators as defined above with appropriate loader arguments
# to be used as finalized relaxation/PES recipes.
relax_pyscf = NexusGeometry(
    scf_relax_job,
    # pyscf_relax.py is configured to output relaxed geometry in relax.xyz
    loader=XyzGeometry({'suffix': 'relax.xyz'})
)
pes_pyscf = NexusPes(
    scf_pes_job,
    # pyscf_pes.py is configured to output SCF energy in energy.dat
    loader=FilesLoader({'suffix': 'energy.dat'})
)
pes_dmc = NexusPes(
    dmc_pes_job,
    # Nexus QmcpackAnalyzer returns DMC energy for the first time-step after walker
    # generation, so at index->1
    loader=QmcPes({'suffix': '/dmc/dmc.in.xml', 'qmc_idx': 1})
)
