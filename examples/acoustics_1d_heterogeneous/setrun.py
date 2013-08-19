# setrun file for 1D acoustics with a material interface

def setrun(claw_pkg='classic'):
    from clawpack.clawutil import data

    # 1D general data object
    rundata = data.ClawRunData(claw_pkg, 1)    # 1 = number of space dimensions

    # Problem-specific data
    probdata = rundata.new_UserData(name='probdata', fname='setprob.data')

    probdata.add_param('ic',    1, 'Initial condition type')
    probdata.add_param('beta', 5., 'Gaussian hump width parameter')
    probdata.add_param('rhol', 1., 'Density left of interface')
    probdata.add_param('cl',   1., 'Sound speed left of interface')
    probdata.add_param('rhor', 4., 'Density right of interface')
    probdata.add_param('cr',  0.5, 'Sound speed right of interface')


    # General Clawpack control data
    clawdata = rundata.clawdata    # Initialized when rundata is created

    # Size of computational domain
    clawdata.lower = -5.
    clawdata.upper =  5.

    # Grid dimensions
    clawdata.num_cells = 500

    # Size of system
    clawdata.num_eqn = 2
    clawdata.num_aux = 2
    clawdata.capa_index = 0    # No capacity function for this problem

    # Output control
    clawdata.output_style = 1
    clawdata.num_output_times = 20
    clawdata.tfinal = 10.

    # Time stepping
    clawdata.dt_initial = 1.
    clawdata.cfl_desired = 0.9
    clawdata.cfl_max = 1.
    clawdata.steps_max = 500
    clawdata.dt_variable = True

    # Details of the numerical method
    clawdata.order = 2
    clawdata.verbosity = 0
    clawdata.source_split = 0    # No source term
    clawdata.use_fwaves = False

    # Waves and limiters
    clawdata.num_waves = 2
    clawdata.limiter = [4, 4]

    # Boundary conditions
    # Zero-order extrapolation at both ends
    clawdata.bc_lower = [1]
    clawdata.bc_upper = [1]

    return rundata


if __name__ == '__main__':
    # Set up run-time parameters and write all data files.
    import sys
    if len(sys.argv) == 2:
        rundata = setrun(sys.argv[1])
    else:
	rundata = setrun()

    rundata.write()
