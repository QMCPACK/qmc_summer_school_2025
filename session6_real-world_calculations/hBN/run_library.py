#!/usr/bin/env python
# nexus imports
from nexus import Job, obj
from nexus import settings
from nexus import linear, loop, vmc, dmc
from qmcpack_input import spindensity
# general settings for nexus
settings(
    pseudo_dir    = './pseudos',
    status_only   = 0,                    # only show status of runs
    generate_only = 0,                    # only make input files
    sleep         = 3,                    # check on runs every 3 secondsa
    machine       = 'perlmutter',         # Perlmutter NERSC machine
    account       = '<account_name>',     # User account name    
    )

def get_dft_settings(kgrid    = None, 
                    ecutwfc   = None,
                    pseudos   = None,
                    start_mag = None,
                    hubbard   = None,
                    tot_magnetization = None):
    

    if settings.machine == 'perlmutter':
        qe_modules = '' # Modules used to build QE
        qe_bin     = '' # QE build directory    
        dft_job = Job(cores  = 4,
                    threads= 1,
                    hours  = 12,
                    app    = qe_bin+'/pw.x',
                    constraints = 'cpu', # Default
                    presub= qe_modules)
        conv_job = Job(cores=1,
                        hours=1,
                        app='/pw2qmcpack.x',
                        constraints = 'cpu',
                        presub=qe_modules)
    else:
        print('Error: Unknown computer for DFT, using {}'.format(settings.machine))
        exit()

    qe_shared = obj(
        job          = dft_job,
        input_type   = 'generic',
        ecutwfc      = ecutwfc,          # DFT planewave energy cutoff
        input_DFT    = 'PBE',            # DFT functional
        conv_thr     = 1e-8,             # SCF convergence threshold
        wf_collect   = True,             # write orbitals
        pseudos      = pseudos,          # QE Pseudopotentials
        start_mag    = start_mag,        # Starting magnetization
        hubbard      = hubbard,          # Hubbard U parameters
        occupations  = 'smearing',       # Occupation scheme
        smearing     = 'gauss',          # Smearing type
        degauss      = 0.001,            # Smearing widt
        tot_magnetization = tot_magnetization
    )

    scf_shared = obj(
        nosym        = False,            # use symmetry
        identifier   = 'scf',            # identifier/file prefix
        calculation  = 'scf',            # perform scf calculation
        kgrid        = kgrid,            # Converged DFT k-grid
        **qe_shared
    )

    nscf_shared = obj(
        nosym        = True,             # don't use symmetry
        identifier   = 'nscf',           # identifier/file prefix
        calculation  = 'nscf',           # perform nscf calculation
        diagonalization = 'cg',          # Diagonalization method
        **qe_shared
    )    

    conv_shared = obj(
            identifier   = 'conv',           # identifier/file prefix
            job          = conv_job,         # Job object for PW2QMCPACK
            write_psir   = False,            # Don't write psir
    )

    return scf_shared, nscf_shared, conv_shared

def get_qmc_settings(system      = None,
                    hybrid_rcut = None,
                    hybrid_lmax = None, 
                    meshfactor  = 1.0,
                    pseudos     = None):
    
    if settings.machine == 'perlmutter':
        qmcpack_modules = ''    # Modules used to build QMCPACK
        qmcpack_bin     = ''    # QMCPACK build directory
        qmcpack_exec    = qmcpack_bin+'/qmcpack_complex'
        opt_job = Job(nodes   = 12,
                    threads = 16,
                    hours   = 12,
                    constraint = 'gpu',
                    app     = qmcpack_exec,
                    presub  = qmcpack_modules)
        dmc_job = Job(nodes   = 24,
                    threads = 16,
                    hours   = 12,
                    constraint = 'gpu',
                    app     = qmcpack_exec,
                    presub  = qmcpack_modules)        
    else:
        print('Error: Unknown computer for QMC, using {}'.format(settings.machine))
        exit()

    system.structure.change_units('B')
    rwigner = system.structure.rwigner()

    qmc_settings = obj(
        system          = system,        # PhysicalSystem object containing structural info
        input_type      = 'basic',       # Simple input format for QMCPACK
        pseudos         = pseudos,       # Pseudopotential files for QMC
        driver          = 'batched',     # Use batched driver in QMCPACK
        hybrid_rcut     = hybrid_rcut,   # Cutoff radius for hybrid orbital representation 
        hybrid_lmax     = hybrid_lmax,   # Max angular momentum for hybrid orbitals
        meshfactor      = meshfactor,    # Controls fineness of real-space (Spline) grid
        lr_handler      = 'ewald',       # Use Ewald summation for long-range interactions
        lr_dim_cutoff   = 30,            # Cutoff for long-range Ewald sums
        spin_polarized  = True,          # Enable spin-polarized calculations
    )

    opt_parameters = obj(
        num_varmin_j2   = 12,     # Number of variance minimization iterations for 2-body Jastrow
        num_emin_j2     = 8,      # Number of energy minimization iterations for 2-body Jastrow
        num_emin_j3     = 6,      # Number of energy minimization iterations for 3-body Jastrow
        j2_init         = "rpa",  # Initialize 2-body Jastrow with Random Phase Approximation
        num_j1_jastrows = 10,     # Number of 1-body Jastrow parameters to optimize
        num_j2_jastrows = 10,     # Number of 2-body Jastrow parameters to optimize
        num_j3_jastrows = 3,      # Number of 3-body Jastrow parameters to optimize
        j3_rcut         = 4.0 if rwigner > 4.0 else rwigner,  # 3-body Jastrow cutoff radius (min of 4.0 or Wigner radius)
        timestep        = 1.0     # VMC timestep for optimization
    )

    opt_settings = obj(
        job             = opt_job,
        twistnum        = 0,     # Twist index
    )
    opt_settings = opt_settings.set(qmc_settings)

    # Variance minimization settings
    varmin = linear(
        energy               = 0.0,                # Weight for energy minimization (0 = pure variance min)
        unreweightedvariance = 1.0,                # Weight for unreweighted variance minimization
        reweightedvariance   = 0.0,                # Weight for reweighted variance minimization
        minwalkers           = 1e-4,               # Lower bound of the effective walker weight
        shift_i              = 0.05,               # (OneShiftOnly Optimizer) Direct stabilizer shift
        shift_s              = 1.0,                # (OneShiftOnly Optimizer) Stabilizer shift based on overlap matrix
        warmupsteps          = 200,                # Number of steps before measurements begin
        blocks               = 100,                # Number of statistical measurement blocks
        steps                = 1,                  # Steps per block
        timestep             = 1.0,                # VMC timestep 
        minmethod            = "OneShiftOnly",     # Minimization algorithm to use
        substeps             = 10,                 # Number of MC steps between parameter updates
    )    

    # Energy minimization settings
    emin = varmin.copy() # Copy from varmin
    emin.minwalkers             = 0.5  # Use larger minwalkers, since varmin provides a better starting point
    emin.energy                 = 0.95 # Mixed cost function 0.95 energy / 0.05 variance
    emin.unreweightedvariance   = 0.0
    emin.reweightedvariance     = 0.05
    emin.shift_i                = 0.01 # Reduced shift_i, since we are closer to the minimum

    j2_settings     = obj(
        calculations = [loop(max=opt_parameters.num_varmin_j2, qmc=varmin), 
                        loop(max=opt_parameters.num_emin_j2,   qmc=emin)],
        J1_size = opt_parameters.num_j1_jastrows, 
        J2_size = opt_parameters.num_j2_jastrows, 
        J1_rcut = rwigner, 
        J2_rcut = rwigner, 
        J2_init = opt_parameters.j2_init,
        **opt_settings
    )

    j3_settings     = obj(
        calculations = [loop(max=opt_parameters.num_emin_j3, qmc=emin)],
        J3=True,
        J3_isize = opt_parameters.num_j3_jastrows,
        J3_esize = opt_parameters.num_j3_jastrows,
        J3_rcut  = opt_parameters.j3_rcut,
        **opt_settings
    )    

    dmc_parameters = obj(
        vmcdt                = 0.3,     # VMC timestep in atomic units
        vmcwarmup            = 25,      # Number of VMC blocks to equilibrate
        vmcblocks            = 100,     # Number of VMC measurement blocks
        vmcsubsteps          = 4,       # VMC steps between measurements
        dmc_eq_dt            = 0.02,    # DMC equilibration timestep
        dmc_eq_blocks        = 100,     # Number of DMC equilibration blocks
        dmcdt                = 0.005,   # DMC production timestep
        dmcblocks            = 500,     # Number of DMC production blocks
        dmcwarmup            = 100,     # Number of DMC blocks to equilibrate
        dmcsteps             = 10,      # Steps per DMC block
        vmc_walkers_per_rank = 240,     # Number of VMC walkers per MPI rank
        dmc_walkers_per_rank = 240,     # Number of DMC walkers per MPI rank
        nonlocalmoves        = False,   # Use T-moves for non-local pseudopotentials
    )
    vmc_dmc = obj(
        warmupsteps = dmc_parameters.vmcwarmup,
        blocks      = dmc_parameters.vmcblocks,
        steps       = 1,
        timestep    = dmc_parameters.vmcdt,
        substeps    = dmc_parameters.vmcsubsteps,
        walkers_per_rank = dmc_parameters.vmc_walkers_per_rank
    )
    dmc_eq  = obj(
        warmupsteps = dmc_parameters.dmcwarmup,
        blocks      = dmc_parameters.dmc_eq_blocks,
        steps       = dmc_parameters.dmcsteps,
        timestep    = dmc_parameters.dmc_eq_dt,
        walkers_per_rank = dmc_parameters.dmc_walkers_per_rank,
        nonlocalmoves = dmc_parameters.nonlocalmoves, 
    )
    dmc_stat = obj(
        warmupsteps = dmc_parameters.dmcwarmup,
        blocks      = dmc_parameters.dmcblocks,
        steps       = dmc_parameters.dmcsteps,
        timestep    = dmc_parameters.dmcdt,
        walkers_per_rank = dmc_parameters.dmc_walkers_per_rank,
        nonlocalmoves = dmc_parameters.nonlocalmoves, 
    )

    dmc_settings = obj(
        job           = dmc_job,
        calculations  = [vmc(**vmc_dmc), dmc(**dmc_eq), dmc(**dmc_stat)],
        estimators    = [spindensity(dr=3*[0.3])],
        **qmc_settings    
    )

    return j2_settings, j3_settings, dmc_settings