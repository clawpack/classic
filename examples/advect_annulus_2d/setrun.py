# setrun file for 2D advection in an annulus

def setrun(claw_pkg='classic'):
    from clawpack.clawutil import clawdata

    # 2D general data object
    rundata = clawdata.ClawRunData(claw_pkg, 2)    # 2 = number of dimensions

    # Problem-specific data
    probdata = rundata.new_UserData(name='probdata', fname='setprob.data')

    probdata.add_param('A1',     1., 'Amplitude of first Gaussian hump')
    probdata.add_param('beta1', 40., 'Decay rate of first Gaussian hump')
    probdata.add_param('x1',   -0.5, 'X coordinate of center of first Gaussian hump')
    probdata.add_param('y1',     0., 'Y coordinate of center of first Gaussian hump')
    probdata.add_param('A2',    -1., 'Amplitude of second Gaussian hump')
    probdata.add_param('beta2', 40., 'Decay rate of second Gaussian hump')
    probdata.add_param('x2',    0.5, 'X coordinate of center of second Gaussian hump')
    probdata.add_param('y2',     0., 'Y coordinate of center of second Gaussian hump')


    # General Clawpack control data
    clawdata = rundata.clawdata    # Initialized when rundata is created

    # Size of computational domain
    clawdata.lower = [0.2, 0.]
    clawdata.upper = [1., 6.28318530718]

    # Grid dimensions
    clawdata.num_cells = [40, 120]    # 40 radial, 120 azimuthal

    # Size of system
    clawdata.num_eqn = 1
    clawdata.num_aux = 3
    clawdata.capa_index = 3

    # Output control
    clawdata.output_style = 1
    clawdata.num_output_times = 25
    clawdata.tfinal = 2.5

    # Time stepping
    clawdata.dt_initial = 0.1
    clawdata.cfl_desired = 0.9
    clawdata.cfl_max = 1.
    clawdata.steps_max = 1000
    clawdata.dt_variable = 1

    # Details of the numerical method
    clawdata.order = 2
    clawdata.transverse_waves = 2
    clawdata.verbosity = 0
    clawdata.source_split = 0

    # Waves and limiting
    clawdata.num_waves = 1
    clawdata.limiter = [3]

    # Boundary conditions
    # Zero-order extrapolation in radial direction, periodic in azimuthal
    clawdata.bc_lower = [1, 2]
    clawdata.bc_upper = [1, 2]

    return rundata


if __name__ == '__main__':
    # Set up run-time parameters and write all data files.
    import sys
    if len(sys.argv) == 2:
        rundata = setrun(sys.argv[1])
    else:
	rundata = setrun()

    rundata.write()
