# setrun file for 3D variable coefficient acoustics example

def setrun(claw_pkg='classic'):
    from clawpack.clawutil import data

    # General data object
    rundata = data.ClawRunData(claw_pkg, 3)    # 3 = number of dimensions

    # Problem-specific data
    probdata = rundata.new_UserData(name='probdata', fname='setprob.data')

    probdata.add_param('z1', 2., 'Impedance in lower part of domain')
    probdata.add_param('c1', 2., 'Sound speed in lower part of domain')
    probdata.add_param('z2', 1., 'Impedance in upper part of domain')
    probdata.add_param('c2', 1., 'Sound speed in upper part of domain')


    # General Clawpack control data
    clawdata = rundata.clawdata    # Initialized when rundata is created

    # Size of computational domain
    clawdata.lower = [-1., -1., -1.]
    clawdata.upper = [ 1.,  1.,  1.]

    # Dimensions of grid
    clawdata.num_cells = [50, 50, 50]

    # Size of system
    clawdata.num_eqn = 4
    clawdata.num_aux = 2
    clawdata.capa_index = 0    # No capacity function

    # Output control
    clawdata.output_style = 1
    clawdata.num_output_times = 6
    clawdata.tfinal = 1.2

    # Time stepping
    clawdata.dt_initial = 0.01
    clawdata.cfl_desired = 0.9
    clawdata.cfl_max = 1.
    clawdata.steps_max = 500
    clawdata.dt_variable = 1

    # Details of the numerical method
    clawdata.order = 2
    clawdata.transverse_waves = 22
    clawdata.verbosity = 1
    clawdata.source_split = 0

    # Waves and limiting
    clawdata.num_waves = 2
    clawdata.limiter = [3]*clawdata.num_waves

    # Boundary conditions
    clawdata.bc_lower = [1, 1, 1]
    clawdata.bc_upper = [1, 1, 0]

    return rundata


if __name__ == '__main__':
    # Set up run-time parameters and write all data files.
    import sys
    if len(sys.argv) == 2:
        rundata = setrun(sys.argv[1])
    else:
	rundata = setrun()

    rundata.write()
