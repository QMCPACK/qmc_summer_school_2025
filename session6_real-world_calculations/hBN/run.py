#!/usr/bin/env python
# user library imports
from run_library import get_dft_settings, get_qmc_settings
# nexus imports
from nexus import run_project, read_structure, obj
from nexus import generate_physical_system
from nexus import generate_pwscf
from nexus import generate_pw2qmcpack
from nexus import generate_qmcpack

# structure files and interlayer separations in Angstroms
structures = {3.0:  'structures/hBN_d_3000.xsf',
            2.5:  'structures/hBN_d_2500.xsf',                
            3.25: 'structures/hBN_d_3250.xsf',
            3.5:  'structures/hBN_d_3500.xsf',
            4.0:  'structures/hBN_d_4000.xsf',
            4.5:  'structures/hBN_d_4500.xsf',
            5.0:  'structures/hBN_d_5000.xsf',
            'mono' : 'structures/hBN_mono.xsf'}
interlayer_separations  = list(structures.keys())

# Supercell tiling vectors and respective kgrids
tiling_vectors          = [(2,2,1), (3,3,1), (4,4,1)]
tiling_kgrids           = {(2,2,1):(4,4,1), 
                           (3,3,1):(2,2,1), 
                           (4,4,1):(2,2,1)}

# DFT and QMC settings shared across all calculations
system_shared = obj(
    B        = 3,        # Boron PP valency
    N        = 5,        # Nitrogen PP valency
    net_spin = 0         # Net spin of the system
)

dft_shared = obj(
    kgrid    = (8,8,1),  # K-point grid for DFT calculations
    ecutwfc  = 400,      # Plane-wave cutoff energy in Rydberg
    pseudos  = 'B.ccECP.upf N.ccECP.upf'.split()  # DFT PP files for Boron and Nitrogen
)

qmc_shared = obj(
    # Hybrid representation cutoff radius for Boron and Nitrogen in atomic units (a.u.)
    hybrid_rcut  = obj(B=1.1, N=1.1),  
    # Maximum angular momentum for hybrid representation for Boron and Nitrogen
    hybrid_lmax  = obj(B=5, N=5),     
    # Blip-spline Mesh factor for QMC calculations
    meshfactor   = 0.5,               
    # Pseudopotential files for QMC calculations
    pseudos      = 'B.ccECP.xml  N.ccECP.xml'.split()  
)

# SCF, NSCF and PW2QMCPACK settings
scf_shared, nscf_shared, conv_shared = get_dft_settings(**dft_shared)

# Binding energy workflow start 
for d in interlayer_separations:
    # Convert interlayer separation to an int for file naming
    if isinstance(d, (int, float)):
        d_name = int(d*1000)
    else:
        d_name = d

    scf_path = 'scf_{}'.format(d_name)
    
    # Generate the primitive cell system
    prim_system = generate_physical_system(
        structure = read_structure(structures[d]),
        **system_shared
        )
    # SCF calculation
    scf_run = generate_pwscf(
        system = prim_system,
        path = scf_path,
        **scf_shared
        )
    for t in tiling_vectors:
        # Directory for the NSCF calculation
        nscf_path = 'nscf_{}_{}'.format(d_name, t[0])
        
        # Generate the supercell system
        tiled_system = generate_physical_system(
            structure = read_structure(structures[d]),
            tiling   = t,
            kgrid    = tiling_kgrids[t],
            **system_shared
        )
        # NSCF calculation
        nscf_run = generate_pwscf(
            system = tiled_system,
            path = nscf_path,
            **nscf_shared
        )        
        # PW2QMCPACK conversion calculation
        conv_run = generate_pw2qmcpack(
            path         = nscf_path,    # Use the same path as the NSCF calculation
            dependencies = (nscf_run, 'orbitals'),
            **conv_shared
        )        

        dmc_path = 'dmc_{}_{}'.format(d_name, t[0])
        
        # Optimize jastrows using the first structure listed in interlayer_separations
        # In this example, this is d == 3.0 since dictionary keys are always ordered in Python 3.7+
        j2_settings, j3_settings, dmc_settings  = get_qmc_settings(system = tiled_system, 
                                                                        **qmc_shared)
        if d == interlayer_separations[0]: 
            j2_path = 'j2_{}_{}'.format(d_name, t[0])
            j3_path = 'j3_{}_{}'.format(d_name, t[0])
            # J2, J3 optimizations and DMC calculation settings
            # Here each "settings" object is specific to the system size
            
            
            # J2 optimization calculation
            j2_run = generate_qmcpack(path = j2_path,
                                    dependencies = (conv_run, 'orbitals'),
                                    **j2_settings)
            # J3 optimization calculation
            j3_run = generate_qmcpack(path = j3_path,
                                    dependencies = [(j2_run, 'jastrow'), (conv_run, 'orbitals')],
                                    **j3_settings)
        # DMC calculation
        dmc_run = generate_qmcpack(path = dmc_path,
                                    dependencies = [(j3_run, 'jastrow'),(conv_run, 'orbitals')],
                                    **dmc_settings)
run_project()
