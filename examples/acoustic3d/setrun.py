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
    clawdata.xlower = -1.
    clawdata.xupper =  1.
    clawdata.ylower = -1.
    clawdata.yupper =  1.
    clawdata.zlower = -1.
    clawdata.zupper =  1.

    # Dimensions of grid
    clawdata.mx = 20
    clawdata.my = 20
    clawdata.mz = 20

    # Size of system
    clawdata.meqn = 4
    clawdata.maux = 2
    clawdata.mcapa = 0    # No capacity function

    # Output control
    clawdata.outstyle = 1
    clawdata.nout = 6
    clawdata.tfinal = 1.2

    # Time stepping
    clawdata.dt_initial = 0.01
    clawdata.cfl_desired = 0.9
    clawdata.cfl_max = 1.
    clawdata.max_steps = 500
    clawdata.dt_variable = 1

    # Details of the numerical method
    clawdata.order = 2
    clawdata.order_trans = 22
    clawdata.verbosity = 1
    clawdata.src_split = 0

    # Waves and limiting
    clawdata.mwaves = 2
    clawdata.mthlim = [3]*clawdata.mwaves

    # Boundary conditions
    clawdata.mthbc_xlower = 1    # Zero-order extrapolation
    clawdata.mthbc_xupper = 1
    clawdata.mthbc_ylower = 1
    clawdata.mthbc_yupper = 1
    clawdata.mthbc_zlower = 1
    clawdata.mthbc_zupper = 0    # Custom BC defined in this problem's bc3 routine

    return rundata


if __name__ == '__main__':
    # Set up run-time parameters and write all data files.
    import sys
    if len(sys.argv) == 2:
        rundata = setrun(sys.argv[1])
    else:
	rundata = setrun()

    rundata.write()
