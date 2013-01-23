program advection1d

    use solver_module
    use solution_module
    use rp1_advection, only: u, rp1

    use clawdata_module

    implicit none
    
    type(solution_type) :: solution
    type(solver_type) :: solver
    type(clawdata_type) :: clawdata

    real(kind=8) :: beta

    ! Read in advective speed and initial condition data
    call open_data_file(13,'./setprob.data')
    read(13,*) u
    read(13,*) beta
    close(13)

    ! Read input data
    call new(clawdata,'./claw.data')

    ! Initialize solution
    call new(solution,clawdata)
    solution%q(1,:) = exp(-beta * (solution%centers - 0.3d0)**2)

    ! Initialize solver
    call new(solver,clawdata)
    solver%rp1 => rp1

    ! Simple single step test
    if (clawdata%output_t0) then
        call output_solution(solution, 0.d0, 0, .false., './_output/')
    endif
    call evolve_to_time(clawdata%t_final,solution,solver)
    call output_solution(solution, clawdata%t_final, 1, .false., './_output/')

end program advection1d
